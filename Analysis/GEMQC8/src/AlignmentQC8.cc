#include "Analysis/GEMQC8/interface/AlignmentQC8.h"

using namespace std;
using namespace edm;

AlignmentQC8::AlignmentQC8(const edm::ParameterSet& cfg): GEMBaseValidation(cfg)
{
  time_t rawTime;
  time(&rawTime);
  printf("Begin of AlignmentQC8::AlignmentQC8() at %s\n", asctime(localtime(&rawTime)));
  isMC = cfg.getParameter<bool>("isMC");
  InputTagToken_RH = consumes<GEMRecHitCollection>(cfg.getParameter<edm::InputTag>("recHitsInputLabel"));
  InputTagToken_TR = consumes<vector<reco::Track>>(cfg.getParameter<edm::InputTag>("tracksInputLabel"));
  InputTagToken_TS = consumes<vector<TrajectorySeed>>(cfg.getParameter<edm::InputTag>("seedInputLabel"));
  InputTagToken_TJ = consumes<vector<Trajectory>>(cfg.getParameter<edm::InputTag>("trajInputLabel"));
  InputTagToken_TI = consumes<vector<int>>(cfg.getParameter<edm::InputTag>("chNoInputLabel"));
  InputTagToken_TT = consumes<vector<unsigned int>>(cfg.getParameter<edm::InputTag>("seedTypeInputLabel"));
  InputTagToken_DG = consumes<GEMDigiCollection>(cfg.getParameter<edm::InputTag>("gemDigiLabel"));
  edm::ParameterSet serviceParameters = cfg.getParameter<edm::ParameterSet>("ServiceParameters");
  theService = new MuonServiceProxy(serviceParameters);
  minCLS = cfg.getParameter<double>("minClusterSize");
  maxCLS = cfg.getParameter<double>("maxClusterSize");
  maxRes = cfg.getParameter<double>("maxResidual");
  SuperChamType = cfg.getParameter<vector<string>>("SuperChamberType");
  vecChamType = cfg.getParameter<vector<double>>("SuperChamberSeedingLayers");
  edm::ParameterSet smootherPSet = cfg.getParameter<edm::ParameterSet>("MuonSmootherParameters");
  theSmoother = new CosmicMuonSmoother(smootherPSet, theService);
  theUpdator = new KFUpdator();
  time(&rawTime);

  edm::Service<TFileService> fs;

  // Histograms declaration

  goodVStriggeredEvts = fs->make<TH1D>("goodVStriggeredEvts","Events with track vs triggered events",2,0,2);
  string name = "";
  for(int i = 0; i < 15; i++) // SuperChamber#
  {
    for(int j = 0; j < 8; j++) // iEta
    {
      name = Form("h_resX_eta_%d_%d",i+1,j+1);
      h_resX_eta[i][j] = fs->make<TH1D>(name.c_str(),name.c_str(),200,-3.0,3.0);
    }
  }
  for(int i = 0; i<3; i++) // column number
  {
    for(int j=0;j<8;j++) // iEta
    {
      name = Form("h_PxPz_col_eta_%d_%d",i+1,j+1);
      h_PxPz_col_eta[i][j] = fs->make<TH1D>(name.c_str(),name.c_str(),500,-0.05,0.05);
    }
  }

  // Tree branches declaration

  tree = fs->make<TTree>("tree", "Tree for QC8");
  tree->Branch("run",&run,"run/I");
  tree->Branch("lumi",&lumi,"lumi/I");
  tree->Branch("ev",&nev,"ev/I");
  tree->Branch("trajPhi",&trajPhi,"trajPhi/F");
  tree->Branch("trajTheta",&trajTheta,"trajTheta/F");
  tree->Branch("trajX",&trajX,"trajX/F");
  tree->Branch("trajY",&trajY,"trajY/F");
  tree->Branch("trajZ",&trajZ,"trajZ/F");
  tree->Branch("trajPx",&trajPx,"trajPx/F");
  tree->Branch("trajPy",&trajPy,"trajPy/F");
  tree->Branch("trajPz",&trajPz,"trajPz/F");
  tree->Branch("nRecHitsTraj",&nRecHitsTraj,"nRecHitsTraj/I");
  tree->Branch("chi2Traj",&chi2Traj,"chi2Traj/F");
  tree->Branch("ndofTraj",&ndofTraj,"ndofTraj/I");
  tree->Branch("testTrajHitX",&testTrajHitX,"testTrajHitX[30]/F");
  tree->Branch("testTrajHitY",&testTrajHitY,"testTrajHitY[30]/F");
  tree->Branch("testTrajHitZ",&testTrajHitZ,"testTrajHitZ[30]/F");
  tree->Branch("testTrajHitXerr",&testTrajHitXerr,"testTrajHitXerr[30]/F");
  tree->Branch("testTrajHitYerr",&testTrajHitYerr,"testTrajHitYerr[30]/F");
  tree->Branch("confTestHitX",&confTestHitX,"confTestHitX[30]/F");
  tree->Branch("confTestHitY",&confTestHitY,"confTestHitY[30]/F");
  tree->Branch("confTestHitZ",&confTestHitZ,"confTestHitZ[30]/F");
  tree->Branch("confTestHitXerr",&confTestHitXerr,"confTestHitXerr[30]/F");
  tree->Branch("confTestHitYerr",&confTestHitYerr,"confTestHitYerr[30]/F");
  tree->Branch("confTestHitClSize",&confTestHitClSize,"confTestHitClSize[30]/I");
  tree->Branch("confTestHitiEta",&confTestHitiEta,"confTestHitiEta[30]/I");

  printf("End of AlignmentQC8::AlignmentQC8() at %s\n", asctime(localtime(&rawTime)));
}

void AlignmentQC8::bookHistograms(DQMStore::IBooker & ibooker, edm::Run const & Run, edm::EventSetup const & iSetup ) {
  GEMGeometry_ = initGeometry(iSetup);
  if ( GEMGeometry_ == nullptr) return ;

  const std::vector<const GEMSuperChamber*>& superChambers_ = GEMGeometry_->superChambers();
  for (auto sch : superChambers_)
  {
    int n_lay = sch->nChambers();
    for (int l=0;l<n_lay;l++)
    {
      gemChambers.push_back(*sch->chamber(l+1));
    }
  }
  n_ch = gemChambers.size();
}

int AlignmentQC8::findIndex(GEMDetId id_) {
  return id_.chamberId().chamber()+id_.chamberId().layer()-2;
}

const GEMGeometry* AlignmentQC8::initGeometry(edm::EventSetup const & iSetup) {
  const GEMGeometry* GEMGeometry_ = nullptr;
  try {
    edm::ESHandle<GEMGeometry> hGeom;
    iSetup.get<MuonGeometryRecord>().get(hGeom);
    GEMGeometry_ = &*hGeom;
  }
  catch( edm::eventsetup::NoProxyException<GEMGeometry>& e) {
    edm::LogError("MuonGEMBaseValidation") << "+++ Error : GEM geometry is unavailable on event loop. +++\n";
    return nullptr;
  }
  return GEMGeometry_;
}

AlignmentQC8::~AlignmentQC8() {
  time_t rawTime;
  time(&rawTime);
  printf("End of AlignmentQC8::AlignmentQC8() at %s\n", asctime(localtime(&rawTime)));
}

void AlignmentQC8::analyze(const edm::Event& e, const edm::EventSetup& iSetup){

  run = e.id().run();
  lumi = e.id().luminosityBlock();
  nev = e.id().event();

  goodVStriggeredEvts->Fill(0);

  for(int i=0;i<30;i++) //loop over the layers of each superchambers
  {
    testTrajHitX[i] = testTrajHitY[i] = testTrajHitZ[i] = -999.9;
    testTrajHitXerr[i] = testTrajHitYerr[i] = -999.9;
    confTestHitX[i] = confTestHitY[i] = confTestHitZ[i] = -999.9;
    confTestHitXerr[i] = confTestHitYerr[i] = -999.9;
    confTestHitClSize[i] = confTestHitiEta[i] = -1;
  }

  trajPhi = trajTheta = -999.9;
  trajX = trajY = trajZ = -999.9;
  trajPx = trajPy = trajPz = -999.9;
  chi2Traj = -999.9;
  ndofTraj = nRecHitsTraj = 0;

  theService->update(iSetup);

  // recHits

  edm::Handle<GEMRecHitCollection> gemRecHits;
  e.getByToken( this->InputTagToken_RH, gemRecHits);
  if (!gemRecHits.isValid())
  {
    edm::LogError("AlignmentQC8") << "Cannot get strips by Token RecHits Token.\n";
    return ;
  }

  edm::Handle<std::vector<int>> idxChTraj;
  e.getByToken( this->InputTagToken_TI, idxChTraj);

  edm::Handle<std::vector<TrajectorySeed>> seedGCM;
  e.getByToken( this->InputTagToken_TS, seedGCM);

  edm::Handle<std::vector<Trajectory>> trajGCM;
  e.getByToken( this->InputTagToken_TJ, trajGCM);

  edm::Handle<vector<reco::Track>> trackCollection;
  e.getByToken( this->InputTagToken_TR, trackCollection);

  edm::Handle<std::vector<unsigned int>> seedTypes;
  e.getByToken( this->InputTagToken_TT, seedTypes);

  if(trackCollection->size() == 0) return;

  // Get the propagators

  edm::ESHandle<Propagator> propagatorAlong;
  edm::ESHandle<Propagator> propagatorOpposite;
  iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny",propagatorAlong);
  iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny",propagatorOpposite);

  // Getting the track info

  TrajectorySeed bestSeed = (*trajGCM->begin()).seed();
  Trajectory bestTraj = *trajGCM->begin();
  reco::Track bestTrack = *trackCollection->begin();

  const FreeTrajectoryState* ftsAtVtx = bestTraj.geometricalInnermostState().freeState();

  GlobalPoint trackPCA = ftsAtVtx->position();
  GlobalVector gvecTrack = ftsAtVtx->momentum();

  trajTheta = gvecTrack.theta();
  trajPhi = gvecTrack.phi();
  trajX = trackPCA.x();
  trajY = trackPCA.y();
  trajZ = trackPCA.z();
  trajPx = gvecTrack.x();
  trajPy = gvecTrack.y();
  trajPz = gvecTrack.z();
  nRecHitsTraj = size(bestTraj.recHits());
  chi2Traj = bestTraj.chiSquared();
  ndofTraj = bestTraj.ndof();

  FreeTrajectoryState startPoint(bestTraj.lastMeasurement().updatedState().globalParameters().position(),bestTraj.lastMeasurement().updatedState().globalParameters().momentum(),bestTrack.charge(),&*theService->magneticField());
  startPoint.setCurvilinearError(bestTraj.lastMeasurement().updatedState().curvilinearError());

  for(int c=0; c<n_ch;c++)
  {
    GEMChamber ch = gemChambers[c];

    int index = findIndex(ch.id());

    // Find the ieta partition ( -> mRoll )

    int mRoll = -1;
    int n_eta = ch.nEtaPartitions();

    TrajectoryStateOnSurface tsosTestProp;
    TrajectoryStateOnSurface tsosProp;

    for (int ieta=0; ieta<n_eta and mRoll==-1; ieta++)
    {
      const BoundPlane& bpeta = GEMGeometry_->idToDet(ch.etaPartition(ieta+1)->id())->surface();
      tsosTestProp = propagatorAlong->propagate(startPoint, bpeta);
      if (!tsosTestProp.isValid()) tsosTestProp = propagatorOpposite->propagate(startPoint, bpeta);
      if (!tsosTestProp.isValid()) continue;
      if (bpeta.bounds().inside(tsosTestProp.localPosition()))
      {
        mRoll = ieta+1;
        tsosProp = tsosTestProp;
      }
    }

    if (mRoll == -1) continue; // The track was not passing through that chamber... (another column, etc...)

    Global3DPoint gtrp = tsosProp.globalPosition();

    testTrajHitX[index] = gtrp.x();
    testTrajHitY[index] = gtrp.y();
    testTrajHitZ[index] = gtrp.z();

    double xx = tsosProp.curvilinearError().matrix()(3,3);
    double yy = tsosProp.curvilinearError().matrix()(4,4);
    double xy = tsosProp.curvilinearError().matrix()(4,3);
    testTrajHitXerr[index] = sqrt(0.5*(xx+yy-sqrt((xx-yy)*(xx-yy)+4*xy*xy)));
    testTrajHitYerr[index] = sqrt(0.5*(xx+yy+sqrt((xx-yy)*(xx-yy)+4*xy*xy)));

    // Check if there's a matching recHit in the test chamber (confRecHit) - confirming recHit. Otherwise put all the other hits outside the range in non asssociated

    GEMRecHit confRecHit;

    double maxR = 99.9;

    for ( GEMRecHitCollection::const_iterator hit = gemRecHits->begin(); hit != gemRecHits->end(); ++hit )
    {
      // cluster size selection
      if ((*hit).clusterSize()<minCLS or (*hit).clusterSize()>maxCLS) continue;

      GEMDetId hitID((*hit).rawId());
      int hitCh = findIndex(hitID);

      if (hitCh == index) // hit in test chamber
      {
        GlobalPoint hitGP = GEMGeometry_->idToDet((*hit).gemId())->surface().toGlobal(hit->localPosition());
        int hitiEta = hitID.roll();

        if (fabs(hitGP.x() - gtrp.x()) > maxRes or abs(hitiEta - mRoll) > 1) continue;

        // Associated hits? Choosing the closest one
        double deltaR = (hitGP - gtrp).mag();
        if (deltaR < maxR)
        {
          confRecHit = *hit;
          maxR = deltaR;
        }
      }
    }

    if (confRecHit.rawId()>0)
    {
      GEMDetId confHitID(confRecHit.rawId());
      int confHitiEta = confHitID.roll();
      confTestHitiEta[index] = confHitiEta;

      Global3DPoint confHitGP = GEMGeometry_->idToDet(confRecHit.gemId())->surface().toGlobal(confRecHit.localPosition());

      confTestHitX[index] = confHitGP.x();
      confTestHitY[index] = confHitGP.y();
      confTestHitZ[index] = confHitGP.z();

      confTestHitClSize[index] = confRecHit.clusterSize();
      confTestHitXerr[index] = confRecHit.localPositionError().xx();
      confTestHitYerr[index] = confRecHit.localPositionError().yy();

      if (fabs(trajPy)<0.03)
      {
        h_resX_eta[int(index/2)][confHitiEta]->Fill(confHitGP.x()-gtrp.x());
        h_PxPz_col_eta[int(index/10)][confHitiEta]->Fill(trajPx/trajPz);
      }
    }
  }
  tree->Fill();
  goodVStriggeredEvts->Fill(1);
}

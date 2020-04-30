#include "Analysis/GEMQC8/interface/ValidationQC8.h"

using namespace std;
using namespace edm;

ValidationQC8::ValidationQC8(const edm::ParameterSet& cfg): GEMBaseValidation(cfg)
{
  time_t rawTime;
  time(&rawTime);
  printf("Begin of ValidationQC8::ValidationQC8() at %s\n", asctime(localtime(&rawTime)));
  isMC = cfg.getParameter<bool>("isMC");
  InputTagToken_RH = consumes<GEMRecHitCollection>(cfg.getParameter<edm::InputTag>("recHitsInputLabel"));
  InputTagToken_TR = consumes<vector<reco::Track>>(cfg.getParameter<edm::InputTag>("tracksInputLabel"));
  InputTagToken_TS = consumes<vector<TrajectorySeed>>(cfg.getParameter<edm::InputTag>("seedInputLabel"));
  InputTagToken_TJ = consumes<vector<Trajectory>>(cfg.getParameter<edm::InputTag>("trajInputLabel"));
  InputTagToken_TI = consumes<vector<int>>(cfg.getParameter<edm::InputTag>("chNoInputLabel"));
  InputTagToken_TT = consumes<vector<unsigned int>>(cfg.getParameter<edm::InputTag>("seedTypeInputLabel"));
  InputTagToken_DG = consumes<GEMDigiCollection>(cfg.getParameter<edm::InputTag>("gemDigiLabel"));
  if ( isMC ) InputTagToken_US = consumes<edm::HepMCProduct>(cfg.getParameter<edm::InputTag>("genVtx"));
  edm::ParameterSet serviceParameters = cfg.getParameter<edm::ParameterSet>("ServiceParameters");
  theService = new MuonServiceProxy(serviceParameters);
  minCLS = cfg.getParameter<double>("minClusterSize");
  maxCLS = cfg.getParameter<double>("maxClusterSize");
  maxRes = cfg.getParameter<double>("maxResidual");
  SuperChamType = cfg.getParameter<vector<string>>("SuperChamberType");
  vecChamType = cfg.getParameter<vector<double>>("SuperChamberSeedingLayers");
  TripEventsPerCh = cfg.getParameter<vector<string>>("tripEvents");
  edm::ParameterSet smootherPSet = cfg.getParameter<edm::ParameterSet>("MuonSmootherParameters");
  theSmoother = new CosmicMuonSmoother(smootherPSet, theService);
  theUpdator = new KFUpdator();
  time(&rawTime);

  edm::Service<TFileService> fs;

  // Histograms declaration

  goodVStriggeredEvts = fs->make<TH1D>("goodVStriggeredEvts","Events with track vs triggered events",2,0,2);
  hitsVFATnum = fs->make<TH3D>("hitsVFATnum","confirmedHits per VFAT (numerator of efficiency)",3,0,3,8,0,8,30,0,30);
  hitsVFATdenom = fs->make<TH3D>("hitsVFATdenom","trajHits per VFAT (denominator of efficiency)",3,0,3,8,0,8,30,0,30);
  digiStrips = fs->make<TH3D>("digiStrips","digi per strip",384,0,384,8,0,8,30,0,30);
  digisPerEvtPerCh = fs->make<TH2D>("digisPerEvtPerCh","digis per event per chamber",30,0,30,20,0,20);
  recHits3D = fs->make<TH3D>("recHits3D","recHits 3D map",200,-100,100,156,-61,95,83,-12,154); // volume defined by the scintillators
  recHits2DPerLayer = fs->make<TH3D>("recHits2DPerLayer","recHits per layer",2000,-100,100,8,0,8,10,0,10);
  associatedHits2DPerLayer = fs->make<TH3D>("associatedHits2DPerLayer","Associated recHits to track per layer",2000,-100,100,8,0,8,10,0,10);
  nonAssociatedHits2DPerLayer = fs->make<TH3D>("nonAssociatedHits2DPerLayer","Non associated recHits to track per layer",2000,-100,100,8,0,8,10,0,10);
  num2DPerLayer = fs->make<TH3D>("num2DPerLayer","Numerator per layer",600,-100,100,8,0,8,10,0,10);
  denom2DPerLayer = fs->make<TH3D>("denom2DPerLayer","Denominator per layer",600,-100,100,8,0,8,10,0,10);
  recHitsPerEvt = fs->make<TH1D>("recHitsPerEvt","recHits per event",1000,0,1000);
  nonAssRecHitsPerEvt = fs->make<TH2D>("nonAssRecHitsPerEvt","Non associated recHits per event",30,0,30,100,0,100);
  clusterSize = fs->make<TH3D>("clusterSize","clusterSize per chamber per eta partition",30,0,30,24,0,24,20,0,20);
  associatedHitsClusterSize = fs->make<TH3D>("associatedHitsClusterSize","clusterSize of associated hits per chamber per eta partition",30,0,30,24,0,24,20,0,20);
  nonAssociatedHitsClusterSize = fs->make<TH3D>("nonAssociatedHitsClusterSize","clusterSize of non associated hits per chamber per eta partition",30,0,30,24,0,24,20,0,20);
  residualPhi = fs->make<TH1D>("residualPhi","residualPhi",400,-5,5);
  residualEta = fs->make<TH1D>("residualEta","residualEta",400,-20,20);
  recHitsPerTrack = fs->make<TH1D>("recHitsPerTrack","recHits per reconstructed track",15,0,15);

  if(isMC)
  {
    genMuAngX = fs->make<TH1D>("genMuAngX","genAngX (XZ plane)",1000,-1,1);
    genMuAngY = fs->make<TH1D>("genMuAngY","genAngY (YZ plane)",1000,-0.8,0.8);
  }

  // Tree branches declaration

  tree = fs->make<TTree>("tree", "Tree for QC8");
  tree->Branch("run",&run,"run/I");
  tree->Branch("lumi",&lumi,"lumi/I");
  tree->Branch("ev",&nev,"ev/I");
  tree->Branch("nTraj",&nTraj,"nTraj/I");
  tree->Branch("trajPhi",&trajPhi,"trajPhi[30]/F");
  tree->Branch("trajTheta",&trajTheta,"trajTheta[30]/F");
  tree->Branch("trajX",&trajX,"trajX[30]/F");
  tree->Branch("trajY",&trajY,"trajY[30]/F");
  tree->Branch("trajZ",&trajZ,"trajZ[30]/F");
  tree->Branch("trajPx",&trajPx,"trajPx[30]/F");
  tree->Branch("trajPy",&trajPy,"trajPy[30]/F");
  tree->Branch("trajPz",&trajPz,"trajPz[30]/F");
  tree->Branch("nRecHitsTraj",&nRecHitsTraj,"nRecHitsTraj[30]/I");
  tree->Branch("chi2Traj",&chi2Traj,"chi2Traj[30]/F");
  tree->Branch("ndofTraj",&ndofTraj,"ndofTraj[30]/I");
  tree->Branch("testTrajHitX",&testTrajHitX,"testTrajHitX[30]/F");
  tree->Branch("testTrajHitY",&testTrajHitY,"testTrajHitY[30]/F");
  tree->Branch("testTrajHitZ",&testTrajHitZ,"testTrajHitZ[30]/F");
  tree->Branch("testTrajHitXerr",&testTrajHitXerr,"testTrajHitXerr[30]/F");
  tree->Branch("testTrajHitYerr",&testTrajHitYerr,"testTrajHitYerr[30]/F");
  tree->Branch("testTrajHitZerr",&testTrajHitZerr,"testTrajHitZerr[30]/F");
  tree->Branch("confTestHitX",&confTestHitX,"confTestHitX[30]/F");
  tree->Branch("confTestHitY",&confTestHitY,"confTestHitY[30]/F");
  tree->Branch("confTestHitZ",&confTestHitZ,"confTestHitZ[30]/F");
  tree->Branch("confTestHitXerr",&confTestHitXerr,"confTestHitXerr[30]/F");
  tree->Branch("confTestHitYerr",&confTestHitYerr,"confTestHitYerr[30]/F");
  tree->Branch("confTestHitClSize",&confTestHitClSize,"confTestHitClSize[30]/I");
  tree->Branch("confTestHitiPhi",&confTestHitiPhi,"confTestHitiPhi[30]/I");
  tree->Branch("confTestHitiEta",&confTestHitiEta,"confTestHitiEta[30]/I");


  if (isMC)
  {
    // Tree for gen events
    genTree = fs->make<TTree>("genTree", "gen info for QC8");
    genTree->Branch("genMuPx",&genMuPx,"genMuPx[30]/F");
    genTree->Branch("genMuPy",&genMuPy,"genMuPy[30]/F");
    genTree->Branch("genMuPz",&genMuPz,"genMuPz[30]/F");
    genTree->Branch("genMuPt",&genMuPt,"genMuPt[30]/F");
    genTree->Branch("genMuPhi",&genMuPhi,"genMuPhi[30]/F");
    genTree->Branch("genMuTheta",&genMuTheta,"genMuTheta[30]/F");
    genTree->Branch("genMuX",&genMuX,"genMuX[30]/F");
    genTree->Branch("genMuY",&genMuY,"genMuY[30]/F");
    genTree->Branch("genMuZ",&genMuZ,"genMuZ[30]/F");
    genTree->Branch("genAngX",&genAngX,"genAngX[30]/F");
    genTree->Branch("genAngY",&genAngY,"genAngY[30]/F");
  }

  printf("End of ValidationQC8::ValidationQC8() at %s\n", asctime(localtime(&rawTime)));
}

void ValidationQC8::bookHistograms(DQMStore::IBooker & ibooker, edm::Run const & Run, edm::EventSetup const & iSetup ) {
  time_t rawTime;
  time(&rawTime);
  printf("Begin of ValidationQC8::bookHistograms() at %s\n", asctime(localtime(&rawTime)));
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
  time(&rawTime);

  printf("End of ValidationQC8::bookHistograms() at %s\n", asctime(localtime(&rawTime)));
}

int ValidationQC8::findIndex(GEMDetId id_)
{
  return id_.chamberId().chamber()+id_.chamberId().layer()-2;
}

int ValidationQC8::findiPhi(float x, float a, float b) {
  float step = fabs(b-a)/3.0;
  if ( x < (min(a,b)+step) ) { return 1;}
  else if ( x < (min(a,b)+2.0*step) ) { return 2;}
  else { return 3;}
}

const GEMGeometry* ValidationQC8::initGeometry(edm::EventSetup const & iSetup) {
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

int g_nEvt = 0;
int g_nNumTrajHit = 0;
int g_nNumMatched = 0;

ValidationQC8::~ValidationQC8() {
  printf("Number of events : %i\n", g_nEvt);
  printf("Number of trajHits : %i (g_nNumTrajHit)\n", g_nNumTrajHit);
  printf("Number of matching trajHits : %i (g_nNumMatched)\n", g_nNumMatched);
  if(g_nNumTrajHit>0) printf("eff : %f\n", double(g_nNumMatched)/g_nNumTrajHit);
}

int g_nNumTest = 0;

void ValidationQC8::analyze(const edm::Event& e, const edm::EventSetup& iSetup){

  g_nEvt++;

  run = e.id().run();
  lumi = e.id().luminosityBlock();
  nev = e.id().event();

  goodVStriggeredEvts->Fill(0);

  // Variables initialization

  nrecHit = 0;
  nTraj = 0;

  for (int i=0; i<30; i++)
  {
    // Not in tree

    nDigisPerCh[i] = 0;
    nOfNonAssHits[i] = 0;

    // In tree

    trajPhi[i] = trajTheta[i] = -999.9;
    trajX[i] = trajY[i] = trajZ[i] = -999.9;
    trajPx[i] = trajPy[i] = trajPz[i] = -999.9;
    chi2Traj[i] = -999.9;
    ndofTraj[i] = nRecHitsTraj[i] = 0;
    testTrajHitX[i] = testTrajHitY[i] = testTrajHitZ[i] = -999.9;
    testTrajHitXerr[i] = testTrajHitYerr[i] = testTrajHitZerr[i] = -999.9;
    confTestHitX[i] = confTestHitY[i] = confTestHitZ[i] = -999.9;
    confTestHitXerr[i] = confTestHitYerr[i] = -999.9;
    confTestHitiEta[i] = confTestHitiPhi[i] = confTestHitClSize[i] = -1;

    // In genTree

    if (isMC)
    {
      genMuPx[i] = -999.9;
      genMuPy[i] = -999.9;
      genMuPz[i] = -999.9;
      genMuPt[i] = -999.9;
      genMuTheta[i] = -999.9;
      genMuPhi[i] = -999.9;
      genMuX[i] = -999.9;
      genMuY[i] = -999.9;
      genMuZ[i] = -999.9;
      genAngX[i] = -999.9;
      genAngY[i] = -999.9;
    }
  }

  // Get the events when a chamber was tripping

  string delimiter = "";
  string line = "";
  string interval = "";
  int ch, beginEvt, endEvt;
  vector<int> beginTripEvt[30];
  vector<int> endTripEvt[30];
  for (unsigned int i = 0; i < TripEventsPerCh.size(); i++)
  {
    line = TripEventsPerCh[i];

    delimiter = ",";
    ch = stoi(line.substr(0, line.find(delimiter)));
    line.erase(0, line.find(delimiter) + delimiter.length());

    int numberOfIntervals = count(line.begin(), line.end(), ',') + 1; // intervals are number of separators + 1
    for (int badInterv = 0; badInterv < numberOfIntervals; badInterv++)
    {
      delimiter = ",";
      interval = line.substr(0, line.find(delimiter));
      delimiter = "-";
      beginEvt = stoi(interval.substr(0, interval.find(delimiter)));
      interval.erase(0, interval.find(delimiter) + delimiter.length());
      endEvt = stoi(interval);
      if (beginEvt < endEvt)
      {
        beginTripEvt[ch].push_back(beginEvt);
        endTripEvt[ch].push_back(endEvt);
      }
      if (endEvt < beginEvt)
      {
        beginTripEvt[ch].push_back(endEvt);
        endTripEvt[ch].push_back(beginEvt);
      }
      delimiter = ",";
      line.erase(0, line.find(delimiter) + delimiter.length());
    }
  }

  theService->update(iSetup);

  // digis

  if (!isMC)
  {
    edm::Handle<GEMDigiCollection> digis;
    e.getByToken( this->InputTagToken_DG, digis);
    int nNumDigi = 0;
    for (GEMDigiCollection::DigiRangeIterator gemdgIt = digis->begin(); gemdgIt != digis->end(); ++gemdgIt)
    {
      nNumDigi++;
      const GEMDetId& gemId = (*gemdgIt).first;
      int chIdRecHit = findIndex(gemId);
      const GEMDigiCollection::Range& range = (*gemdgIt).second;
      for ( auto digi = range.first; digi != range.second; ++digi )
      {
        digiStrips->Fill(digi->strip(),gemId.roll()-1,chIdRecHit); // Strip#=[0,383] -> OK , Eta#=[1,8] -> -1
        nDigisPerCh[chIdRecHit]++;
      }
    }
  }

  for (int i=0; i<30; i++)
  {
    digisPerEvtPerCh->Fill(i,nDigisPerCh[i]);
  }

  // recHits

  edm::Handle<GEMRecHitCollection> gemRecHits;
  e.getByToken(this->InputTagToken_RH, gemRecHits);
  if (!gemRecHits.isValid())
  {
    edm::LogError("ValidationQC8") << "Cannot get strips by Token RecHits Token.\n";
    return ;
  }

  for ( GEMRecHitCollection::const_iterator rechit = gemRecHits->begin(); rechit != gemRecHits->end(); ++rechit )
  {
    // calculation of chamber id
    GEMDetId hitID((*rechit).rawId());
    int chIdRecHit = findIndex(hitID);

    // cluster size plot per VFAT
    int rhEta = hitID.roll();

    GEMChamber ch = gemChambers[0];
    GEMChamber chHit = gemChambers[0];

    for(int c=0; c<n_ch;c++)
    {
      ch = gemChambers[c];
      if (findIndex(ch.id()) == chIdRecHit)
      {
        chHit = ch;
        break;
      }
    }

    int n_strip = chHit.etaPartition(rhEta)->nstrips();
    double min_x = chHit.etaPartition(rhEta)->centreOfStrip(0).x();
    double max_x = chHit.etaPartition(rhEta)->centreOfStrip(n_strip-1).x();

    Local3DPoint rhLP = rechit->localPosition();
    int rhPhi = findiPhi(rhLP.x(), min_x, max_x);

    clusterSize->Fill(chIdRecHit,8-rhEta+8*(rhPhi-1),(*rechit).clusterSize());

    // cluster size selection
    if ((*rechit).clusterSize()<minCLS) continue;
    if ((*rechit).clusterSize()>maxCLS) continue;

    GlobalPoint recHitGP = GEMGeometry_->idToDet((*rechit).gemId())->surface().toGlobal(rechit->localPosition());
    recHits3D->Fill(recHitGP.x(),recHitGP.y(),recHitGP.z());

    recHits2DPerLayer->Fill(recHitGP.x(),hitID.roll()-1,chIdRecHit%10);

    nrecHit++;
  }

  recHitsPerEvt->Fill(nrecHit);

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

  if ( idxChTraj->size() == 0 ) return;

  // Get the propagators

  edm::ESHandle<Propagator> propagatorAlong;
  edm::ESHandle<Propagator> propagatorOpposite;
  iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny",propagatorAlong);
  iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny",propagatorOpposite);

  // Efficiency calculation per chamber from here on!

  int countTC = 0;

  for (auto ch : gemChambers)
  {
    countTC += 1;

    // Select trajectory correspondent to test chamber excluded from fit

    std::vector<int>::const_iterator it1 = idxChTraj->begin();
    std::vector<TrajectorySeed>::const_iterator it2 = seedGCM->begin();
    std::vector<Trajectory>::const_iterator it3 = trajGCM->begin();
    std::vector<reco::Track>::const_iterator it4 = trackCollection->begin();

    TrajectorySeed bestSeed;
    Trajectory bestTraj;
    reco::Track bestTrack;

    for ( ; it1 != idxChTraj->end() ; ) {
      if ( *it1 == countTC ) {
        bestTraj = *it3;
        bestSeed = (*it3).seed();
        bestTrack = *it4;
        break;
      }
      it1++;
      it2++;
      it3++;
      it4++;
    }

    if ( it1 == idxChTraj->end() ) continue;

    const FreeTrajectoryState* ftsAtVtx = bestTraj.geometricalInnermostState().freeState();

    GlobalPoint trackPCA = ftsAtVtx->position();
    GlobalVector gvecTrack = ftsAtVtx->momentum();

    nTraj++;

    int index = findIndex(ch.id());

    trajTheta[index] = gvecTrack.theta();
    trajPhi[index] = gvecTrack.phi();
    trajX[index] = trackPCA.x();
    trajY[index] = trackPCA.y();
    trajZ[index] = trackPCA.z();
    trajPx[index] = gvecTrack.x();
    trajPy[index] = gvecTrack.y();
    trajPz[index] = gvecTrack.z();
    nRecHitsTraj[index] = size(bestTraj.recHits());
    recHitsPerTrack->Fill(size(bestTraj.recHits()));
    chi2Traj[index] = bestTraj.chiSquared();
    ndofTraj[index] = bestTraj.ndof();

    if (isMC)
    {
      HepMC::GenParticle *genMuon = NULL;

      edm::Handle<edm::HepMCProduct> genVtx;
      e.getByToken( this->InputTagToken_US, genVtx);
      genMuon = genVtx->GetEvent()->barcode_to_particle(1);

      genMuPx[index] = float(genMuon->momentum().x());
      genMuPy[index] = float(genMuon->momentum().y());
      genMuPz[index] = float(genMuon->momentum().z());
      genMuPt[index] = float(genMuon->momentum().perp());
      genMuPhi[index] = float(genMuon->momentum().phi());
      genMuTheta[index] = float(genMuon->momentum().theta());
      genAngX[index] = atan(genMuPx[index]/genMuPz[index]);
      genAngY[index] = atan(genMuPy[index]/genMuPz[index]);

      genMuAngX->Fill(genAngX[index]);
      genMuAngY->Fill(genAngY[index]);

      float dUnitGen = 0.1;

      genMuX[index] = float(dUnitGen * genMuon->production_vertex()->position().x());
      genMuY[index] = float(dUnitGen * genMuon->production_vertex()->position().y());
      genMuZ[index] = float(dUnitGen * genMuon->production_vertex()->position().z());

      genTree->Fill();
    }

    FreeTrajectoryState startPoint(bestTraj.lastMeasurement().updatedState().globalParameters().position(),bestTraj.lastMeasurement().updatedState().globalParameters().momentum(),bestTrack.charge(),&*theService->magneticField());
    startPoint.setCurvilinearError(bestTraj.lastMeasurement().updatedState().curvilinearError());

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
    Local3DPoint tlp = tsosProp.localPosition();

    // Find the iphi partition

    int n_strip = ch.etaPartition(mRoll)->nstrips();
    double min_x = ch.etaPartition(mRoll)->centreOfStrip(0).x();
    double max_x = ch.etaPartition(mRoll)->centreOfStrip(n_strip-1).x();

    int iPhi = findiPhi(tlp.x(), min_x, max_x);

    // For testing the edge eta partition on the top and bottom layers only vertical seeds are allowed!

    int topSeedIndex = findIndex(bestTraj.firstMeasurement().recHit()->detUnit()->geographicalId().rawId());
    int bottomSeedIndex = findIndex(bestTraj.lastMeasurement().recHit()->detUnit()->geographicalId().rawId());

    if ((vecChamType[index]==1 or vecChamType[index]==2) and
        (mRoll==1 or mRoll==8) and
        (fabs(atan(gvecTrack.y()/gvecTrack.z()))>0.15)) continue;

    int seedsDiffCol = abs(int(topSeedIndex/10.0) - int(bottomSeedIndex/10.0));

    if ((tlp.x()<(min_x + 1.5) or tlp.x()>(max_x - 1.5)) and seedsDiffCol!=0) continue;

    bool validEvent = true;
    for (unsigned int i = 0; i < beginTripEvt[index].size(); i++)
    {
      if (beginTripEvt[index].at(i) <= nev && nev <= endTripEvt[index].at(i))
      {
        validEvent = false;
      }
    }

    if (validEvent)
    {
      testTrajHitX[index] = gtrp.x();
      testTrajHitY[index] = gtrp.y();
      testTrajHitZ[index] = gtrp.z();

      double a = tsosProp.curvilinearError().matrix()(3,3);
      double b = tsosProp.curvilinearError().matrix()(4,4);
      double c = tsosProp.curvilinearError().matrix()(4,3);
      testTrajHitXerr[index] = sqrt(0.5*(a+b-sqrt((a-b)*(a-b)+4*c*c)));
      testTrajHitYerr[index] = sqrt(0.5*(a+b+sqrt((a-b)*(a-b)+4*c*c)));
      testTrajHitZerr[index] = 0.0;

      g_nNumTrajHit++;

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
          int hitiPhi = findiPhi(hit->localPosition().x(), ch.etaPartition(hitiEta)->centreOfStrip(0).x(), ch.etaPartition(hitiEta)->centreOfStrip(n_strip-1).x());

          // Non associated hits? Fill and continue

          if (fabs(hitGP.x() - gtrp.x()) > maxRes or fabs(hitiEta - mRoll) > 1)
          {
            if (fabs(hitiEta - mRoll) > 0)
            {
              nonAssociatedHitsClusterSize->Fill(hitCh,8-hitiEta+8*(hitiPhi-1),(*hit).clusterSize());
              nonAssociatedHits2DPerLayer->Fill(hitGP.x(),hitiEta-1,index%10);
              nOfNonAssHits[index]++;
            }
            continue;
          }

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
        // All confRecHit details

        GEMDetId confHitID(confRecHit.rawId());
        int confHitiEta = confHitID.roll();
        confTestHitiEta[index] = confHitiEta;

        int n_strip = ch.etaPartition(confHitiEta)->nstrips();
        double min_x = ch.etaPartition(confHitiEta)->centreOfStrip(0).x();
        double max_x = ch.etaPartition(confHitiEta)->centreOfStrip(n_strip-1).x();

        int confHitiPhi = findiPhi(confRecHit.localPosition().x(), min_x, max_x);
        confTestHitiPhi[index] = confHitiPhi;

        Global3DPoint confHitGP = GEMGeometry_->idToDet(confRecHit.gemId())->surface().toGlobal(confRecHit.localPosition());

        confTestHitX[index] = confHitGP.x();
        confTestHitY[index] = confHitGP.y();
        confTestHitZ[index] = confHitGP.z();

        associatedHitsClusterSize->Fill(index,8-confHitiEta+8*(confHitiPhi-1),confRecHit.clusterSize());
        confTestHitClSize[index] = confRecHit.clusterSize();
        confTestHitXerr[index] = confRecHit.localPositionError().xx();
        confTestHitYerr[index] = confRecHit.localPositionError().yy();

        hitsVFATdenom->Fill(confHitiPhi-1,confHitiEta-1,index);
        denom2DPerLayer->Fill(confHitGP.x(),confHitiEta-1,index%10);

        hitsVFATnum->Fill(confHitiPhi-1,confHitiEta-1,index);
        num2DPerLayer->Fill(confHitGP.x(),confHitiEta-1,index%10);

        g_nNumMatched++;

        residualPhi->Fill(confHitGP.x()-gtrp.x());
        residualEta->Fill(confHitGP.y()-gtrp.y());

        associatedHits2DPerLayer->Fill(confHitGP.x(),confHitiEta-1,index%10);
      }
      if (!(confRecHit.rawId()>0))
      {
        hitsVFATdenom->Fill(iPhi-1,mRoll-1,index);
        denom2DPerLayer->Fill(gtrp.x(),mRoll-1,index%10);
      }
      nonAssRecHitsPerEvt->Fill(index,nOfNonAssHits[index]);
    }
  }

  g_nNumTest++;

  tree->Fill();
  goodVStriggeredEvts->Fill(1);
}

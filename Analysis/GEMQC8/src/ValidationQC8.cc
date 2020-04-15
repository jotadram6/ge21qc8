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
    chi2Traj[i] = 0.0;
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
      int chIdRecHit = gemId.chamberId().chamber()+gemId.chamberId().layer()-2;
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
  	int chIdRecHit = hitID.chamberId().chamber() + hitID.chamberId().layer() - 2;

  	// cluster size plot per VFAT
    int rhEta = hitID.roll();

    GEMChamber ch = gemChambers[0];
    GEMChamber chHit = gemChambers[0];

  	for(int c=0; c<n_ch;c++)
    {
      ch = gemChambers[c];
      if (ch.id().chamber()+ch.id().layer()-2 == chIdRecHit)
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

  int countTC = 0;

  for (auto tch : gemChambers)
  {
    countTC += 1;

    // Create collection of recHits in the test chamber

    MuonTransientTrackingRecHit::MuonRecHitContainer testRecHits;

    for (auto etaPart : tch.etaPartitions())
    {
      GEMDetId etaPartID = etaPart->id();
      GEMRecHitCollection::range range = gemRecHits->get(etaPartID);
      for (GEMRecHitCollection::const_iterator rechit = range.first; rechit!=range.second; ++rechit)
      {
        const GeomDet* geomDet(etaPart);
        if ((*rechit).clusterSize()<minCLS) continue;
        if ((*rechit).clusterSize()>maxCLS) continue;
        testRecHits.push_back(MuonTransientTrackingRecHit::specificBuild(geomDet,&*rechit));
      }
    }

    // Select trajectory correspondent to test chamber excluded from fit

    std::vector<int>::const_iterator it1 = idxChTraj->begin();
    std::vector<TrajectorySeed>::const_iterator it2 = seedGCM->begin();
    std::vector<Trajectory>::const_iterator it3 = trajGCM->begin();
    std::vector<reco::Track>::const_iterator it4 = trackCollection->begin();
    std::vector<unsigned int>::const_iterator it5 = seedTypes->begin();

    TrajectorySeed bestSeed;
    Trajectory bestTraj;
    reco::Track bestTrack;
    unsigned int unTypeSeed = 0;

    for ( ; it1 != idxChTraj->end() ; ) {
      if ( *it1 == countTC ) {
        bestTraj = *it3;
        bestSeed = (*it3).seed();
        bestTrack = *it4;
        unTypeSeed = *it5;
        break;
      }
      it1++;
      it2++;
      it3++;
      it4++;
      it5++;
    }

    if ( it1 == idxChTraj->end() ) continue;

    const FreeTrajectoryState* ftsAtVtx = bestTraj.geometricalInnermostState().freeState();

    GlobalPoint trackPCA = ftsAtVtx->position();
    GlobalVector gvecTrack = ftsAtVtx->momentum();

    PTrajectoryStateOnDet ptsd1(bestSeed.startingState());
    DetId did(ptsd1.detId());
    const BoundPlane& bp = theService->trackingGeometry()->idToDet(did)->surface();
    TrajectoryStateOnSurface tsosCurrent = trajectoryStateTransform::transientState(ptsd1,&bp,&*theService->magneticField());

    nTraj++;

    int chIndex = findIndex(tch.id());

    trajTheta[chIndex] = gvecTrack.theta();
    trajPhi[chIndex] = gvecTrack.phi();
    trajX[chIndex] = trackPCA.x();
    trajY[chIndex] = trackPCA.y();
    trajZ[chIndex] = trackPCA.z();
    trajPx[chIndex] = gvecTrack.x();
    trajPy[chIndex] = gvecTrack.y();
    trajPz[chIndex] = gvecTrack.z();
    nRecHitsTraj[chIndex] = size(bestTraj.recHits());
    recHitsPerTrack->Fill(size(bestTraj.recHits()));
    chi2Traj[chIndex] = bestTraj.chiSquared();
    ndofTraj[chIndex] = bestTraj.ndof();

    if (isMC)
    {
      HepMC::GenParticle *genMuon = NULL;

      edm::Handle<edm::HepMCProduct> genVtx;
      e.getByToken( this->InputTagToken_US, genVtx);
      genMuon = genVtx->GetEvent()->barcode_to_particle(1);

      genMuPx[chIndex] = float(genMuon->momentum().x());
      genMuPy[chIndex] = float(genMuon->momentum().y());
      genMuPz[chIndex] = float(genMuon->momentum().z());
      genMuPt[chIndex] = float(genMuon->momentum().perp());
      genMuPhi[chIndex] = float(genMuon->momentum().phi());
      genMuTheta[chIndex] = float(genMuon->momentum().theta());
      genAngX[chIndex] = atan(genMuPx[chIndex]/genMuPz[chIndex]);
      genAngY[chIndex] = atan(genMuPy[chIndex]/genMuPz[chIndex]);

      genMuAngX->Fill(genAngX[chIndex]);
      genMuAngY->Fill(genAngY[chIndex]);

      float dUnitGen = 0.1;

      genMuX[chIndex] = float(dUnitGen * genMuon->production_vertex()->position().x());
      genMuY[chIndex] = float(dUnitGen * genMuon->production_vertex()->position().y());
      genMuZ[chIndex] = float(dUnitGen * genMuon->production_vertex()->position().z());

      genTree->Fill();
    }

    // Extrapolation to all the chambers, test chamber selected for efficiency calculation

    for(int c=0; c<n_ch;c++)
    {
      GEMChamber ch = gemChambers[c];
      const BoundPlane& bpch = GEMGeometry_->idToDet(ch.id())->surface();
      tsosCurrent = theService->propagator("SteppingHelixPropagatorAny")->propagate(tsosCurrent, bpch);
      if (!tsosCurrent.isValid()) continue;

      if (ch==tch)
      {
        double propPointX = (bestTraj.firstMeasurement().updatedState().globalPosition().x()-bestTraj.lastMeasurement().updatedState().globalPosition().x())/(bestTraj.firstMeasurement().updatedState().globalPosition().z()-bestTraj.lastMeasurement().updatedState().globalPosition().z())*(tsosCurrent.freeTrajectoryState()->position().z()-bestTraj.lastMeasurement().updatedState().globalPosition().z())+bestTraj.lastMeasurement().updatedState().globalPosition().x();
        double propPointY = (bestTraj.firstMeasurement().updatedState().globalPosition().y()-bestTraj.lastMeasurement().updatedState().globalPosition().y())/(bestTraj.firstMeasurement().updatedState().globalPosition().z()-bestTraj.lastMeasurement().updatedState().globalPosition().z())*(tsosCurrent.freeTrajectoryState()->position().z()-bestTraj.lastMeasurement().updatedState().globalPosition().z())+bestTraj.lastMeasurement().updatedState().globalPosition().y();

        Global3DPoint gtrp = Global3DPoint(propPointX,propPointY,tsosCurrent.freeTrajectoryState()->position().z());

        Local3DPoint tlp = bpch.toLocal(gtrp);
        if (!bpch.bounds().inside(tlp)) continue;

        // Find the ieta partition ( -> mRoll )

        int n_roll = ch.nEtaPartitions();
        double minDeltaY = 50.;
        int mRoll = -1;
        for (int r=0; r<n_roll; r++)
        {
          const BoundPlane& bproll = GEMGeometry_->idToDet(ch.etaPartition(r+1)->id())->surface();
          Local3DPoint rtlp = bproll.toLocal(gtrp);
          if (minDeltaY > fabs(rtlp.y()))
          {
            minDeltaY = fabs(rtlp.y());
            mRoll = r+1;
          }
        }

        if (mRoll == -1)
        {
          cout << "no mRoll" << endl;
          continue;
        }

        // Define region 'inside' the ieta of the chamber

        int n_strip = ch.etaPartition(mRoll)->nstrips();
        double min_x = ch.etaPartition(mRoll)->centreOfStrip(0).x();
        double max_x = ch.etaPartition(mRoll)->centreOfStrip(n_strip-1).x();

        if ( (tlp.x()>(min_x)) & (tlp.x() < (max_x)) )
        {
          // For testing the edge eta partition on the top and bottom layers only vertical seeds are allowed!

          if ( ( vecChamType[ countTC - 1 ] == 2 || vecChamType[ countTC - 1 ] == 1 ) &&
              ( mRoll == 1 || mRoll == 8 ) &&
              ( unTypeSeed & QC8FLAG_SEEDINFO_MASK_REFVERTROLL18 ) == 0 ) continue;

          uint32_t unDiffCol = ( unTypeSeed & QC8FLAG_SEEDINFO_MASK_DIFFCOL ) >> QC8FLAG_SEEDINFO_SHIFT_DIFFCOL;

          if ( ! ( (tlp.x()>(min_x + 1.5)) & (tlp.x() < (max_x - 1.5)) ) )
          {
            if ( unDiffCol != 0 )
            {
              continue;
            }
            else if ( ( vecChamType[ countTC - 1 ] == 2 || vecChamType[ countTC - 1 ] == 1 ) )
            {
              continue;
            }
          }

          int index = findIndex(ch.id());
          int iPhi = findiPhi(tlp.x(), min_x, max_x);

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

            g_nNumTrajHit++;

            // Check if there's a matching recHit in the test chamber (tmpRecHit)

            double maxR = 99.9;
            shared_ptr<MuonTransientTrackingRecHit> tmpRecHit;

            for (auto hit : testRecHits)
            {
              GEMDetId hitID(hit->rawId());
              if (hitID.chamberId() != ch.id()) continue;

              GlobalPoint hitGP = hit->globalPosition();

              if (fabs(hitGP.x() - gtrp.x()) > maxRes) continue;
              if (fabs(hitID.roll() - mRoll) > 1) continue;

              // Choosing the closest one

              double deltaR = (hitGP - gtrp).mag();
              if (deltaR < maxR)
              {
                tmpRecHit = hit;
                maxR = deltaR;
              }
            }

            if(tmpRecHit)
            {
              double a = 0.0, b = 0.0, c = 0.0;

              Global3DPoint tempHitGP = tmpRecHit->globalPosition();

              testTrajHitX[index] = (bestTraj.firstMeasurement().updatedState().globalPosition().x()-bestTraj.lastMeasurement().updatedState().globalPosition().x())/(bestTraj.firstMeasurement().updatedState().globalPosition().z()-bestTraj.lastMeasurement().updatedState().globalPosition().z())*(tempHitGP.z()-bestTraj.lastMeasurement().updatedState().globalPosition().z())+bestTraj.lastMeasurement().updatedState().globalPosition().x();
              testTrajHitY[index] = (bestTraj.firstMeasurement().updatedState().globalPosition().y()-bestTraj.lastMeasurement().updatedState().globalPosition().y())/(bestTraj.firstMeasurement().updatedState().globalPosition().z()-bestTraj.lastMeasurement().updatedState().globalPosition().z())*(tempHitGP.z()-bestTraj.lastMeasurement().updatedState().globalPosition().z())+bestTraj.lastMeasurement().updatedState().globalPosition().y();
              testTrajHitZ[index] = tempHitGP.z();

              a = bestTraj.firstMeasurement().updatedState().curvilinearError().matrix()(3,3);
              b = bestTraj.firstMeasurement().updatedState().curvilinearError().matrix()(4,4);
              c = bestTraj.firstMeasurement().updatedState().curvilinearError().matrix()(4,3);
              testTrajHitXerr[index] = sqrt(0.5*(a+b-sqrt((a-b)*(a-b)+4*c*c)));
              testTrajHitYerr[index] = sqrt(0.5*(a+b+sqrt((a-b)*(a-b)+4*c*c)));
              testTrajHitZerr[index] = 0.0;

              confTestHitX[index] = tempHitGP.x();
              confTestHitY[index] = tempHitGP.y();
              confTestHitZ[index] = tempHitGP.z();

              n_roll = ch.nEtaPartitions();
        	    minDeltaY = 50.;
              int hitRoll = -1;
              for (int r=0; r<n_roll; r++)
              {
                const BoundPlane& bproll = GEMGeometry_->idToDet(ch.etaPartition(r+1)->id())->surface();
                Local3DPoint hittlp = bproll.toLocal(tempHitGP);
                if (minDeltaY > fabs(hittlp.y()))
                {
                  minDeltaY = fabs(hittlp.y());
                  hitRoll = r+1;
                }
              }

              int n_strip = ch.etaPartition(hitRoll)->nstrips();
          	  double min_x = ch.etaPartition(hitRoll)->centreOfStrip(0).x();
          	  double max_x = ch.etaPartition(hitRoll)->centreOfStrip(n_strip-1).x();

          	  Local3DPoint tempHitLP = tmpRecHit->localPosition();

          	  int hitiPhi = findiPhi(tempHitLP.x(), min_x, max_x);

              hitsVFATdenom->Fill(hitiPhi-1,hitRoll-1,index);
              denom2DPerLayer->Fill(tempHitGP.x(),hitRoll-1,index%10);

              hitsVFATnum->Fill(hitiPhi-1,hitRoll-1,index);
              num2DPerLayer->Fill(tempHitGP.x(),hitRoll-1,index%10);

              g_nNumMatched++;

              residualPhi->Fill(tempHitGP.x()-gtrp.x());
              residualEta->Fill(tempHitGP.y()-gtrp.y());

              GEMDetId confirmedHitID((*tmpRecHit).rawId());
              int chConfHit = confirmedHitID.chamber()+confirmedHitID.layer()-2;
              int etaConfHit = confirmedHitID.roll()-1;

              confTestHitiPhi[chConfHit] = hitiPhi;
              confTestHitiEta[chConfHit] = hitRoll;

              associatedHits2DPerLayer->Fill(tempHitGP.x(),etaConfHit,chConfHit%10);

              // to find the (non) associated hits cluster size, need to match tempHit with corresponding recHit, since tempHit (being a trackHit) has no info about it

              for ( GEMRecHitCollection::const_iterator rechit = gemRecHits->begin(); rechit != gemRecHits->end(); ++rechit )
              {
                GEMDetId recHitID((*rechit).rawId());
                int recHitCh = recHitID.chamber()+recHitID.layer()-2;

                if (recHitCh==chConfHit)
                {
                  GlobalPoint rechitGP = GEMGeometry_->idToDet((*rechit).gemId())->surface().toGlobal(rechit->localPosition());

                  if (fabs(rechitGP.x()-tempHitGP.x())<0.01 && fabs(rechitGP.y()-tempHitGP.y())<0.01)
                  {
                    associatedHitsClusterSize->Fill(recHitCh,8-hitRoll+8*(hitiPhi-1),(*rechit).clusterSize()); // once matched, ieta and iphi are the same and we can use the ones above
                    confTestHitClSize[recHitCh] = (*rechit).clusterSize();
                    confTestHitXerr[recHitCh] = (*rechit).localPositionError().xx();
                    confTestHitYerr[recHitCh] = (*rechit).localPositionError().yy();
                  }
                  else if (fabs(rechitGP.x()-tempHitGP.x())>maxRes && fabs(rechitGP.y()-tempHitGP.y())>1.0)
                  {
                    nonAssociatedHitsClusterSize->Fill(recHitCh,8-hitRoll+8*(hitiPhi-1),(*rechit).clusterSize()); // once matched, ieta and iphi are the same and we can use the ones above
                    nonAssociatedHits2DPerLayer->Fill(rechitGP.x(),hitRoll-1,recHitCh%10);
                    nOfNonAssHits[recHitCh]++;
                  }
                }
              }
            }
            if(!tmpRecHit)
            {
              hitsVFATdenom->Fill(iPhi-1,mRoll-1,index);
              denom2DPerLayer->Fill(gtrp.x(),mRoll-1,index%10);
            }
            nonAssRecHitsPerEvt->Fill(index,nOfNonAssHits[index]);
          }
        }
      }
    }
  }

  g_nNumTest++;

  tree->Fill();
  goodVStriggeredEvts->Fill(1);
}

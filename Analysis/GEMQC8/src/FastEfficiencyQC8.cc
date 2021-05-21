
#include "Analysis/GEMQC8/interface/FastEfficiencyQC8.h"

using namespace std;
using namespace edm;

FastEfficiencyQC8::FastEfficiencyQC8(const edm::ParameterSet& cfg): GEMBaseValidation(cfg)
{
    time_t rawTime;
    time(&rawTime);
    printf("Begin of FastEfficiencyQC8::FastEfficiencyQC8() at %s\n", asctime(localtime(&rawTime)));
    InputTagToken_RH = consumes<GEMRecHitCollection>(cfg.getParameter<edm::InputTag>("recHitsInputLabel"));
    edm::ParameterSet serviceParameters = cfg.getParameter<edm::ParameterSet>("ServiceParameters");
    theService = new MuonServiceProxy(serviceParameters);
    minCLS = cfg.getParameter<double>("minClusterSize");
    maxCLS = cfg.getParameter<double>("maxClusterSize");
    TripEventsPerCh = cfg.getParameter<vector<string>>("tripEvents");
    theUpdator = new KFUpdator();
    time(&rawTime);

    edm::Service<TFileService> fs;

    // Histograms declaration

    recHits3D = fs->make<TH3D>("recHits3D","recHits 3D map",200,-100,100,156,-61,95,83,-12,154); // volume defined by the scintillators
    recHits2DPerLayer = fs->make<TH3D>("recHits2DPerLayer","recHits per layer",400,-100,100,8,0,8,10,0,10);
    clusterSize = fs->make<TH3D>("clusterSize","clusterSize per chamber per eta partition",30,0,30,8,0,8,20,0,20);
    occupancyIfConfirmedHits = fs->make<TH3D>("occupancyIfConfirmedHits","occupancy per strip if hits found in both the chambers",384,0,384,8,0,8,30,0,30);
    DxCorrespondingRecHits = fs->make<TH1D>("DxCorrespondingRecHits","Delta x of corresponding recHits",600,-30,30);
    DiEtaCorrespondingRecHits = fs->make<TH1D>("DiEtaCorrespondingRecHits","Delta iEta of corresponding recHits",16,-8,8);

    // nRecHitsPerEventPerChamber: evolution in time
    nRecHitsPerEvtPerCh = fs->make<TH3D>("nRecHitsPerEvtPerCh","recHits per ieta per ch vs event (packages of 500 evts = 5 sec)",48000,0,24000000,30,0,30,8,0,8);

    // Numerator and denominator
    numerator = fs->make<TH1D>("numerator","numerator",30,0,30);
    denominator = fs->make<TH1D>("denominator","denominator",30,0,30);

    // Numerator and denominator: evolution in time
    numeratorPerEvt = fs->make<TH2D>("numeratorPerEvt","numerator per ch vs event (packages of 6k evts = 1 min)",4000,0,24000000,30,-0.5,29.5);
    denominatorPerEvt = fs->make<TH2D>("denominatorPerEvt","denominator per ch vs event (packages of 6k evts = 1 min)",4000,0,24000000,30,-0.5,29.5);

    // Tree branches declaration

    tree = fs->make<TTree>("tree", "Tree for QC8");
    tree->Branch("run",&run,"run/I");
    tree->Branch("lumi",&lumi,"lumi/I");
    tree->Branch("ev",&nev,"ev/I");

    printf("End of FastEfficiencyQC8::FastEfficiencyQC8() at %s\n", asctime(localtime(&rawTime)));
}

void FastEfficiencyQC8::bookHistograms(DQMStore::IBooker & ibooker, edm::Run const & Run, edm::EventSetup const & iSetup )
{
    time_t rawTime;
    time(&rawTime);
    printf("Begin of FastEfficiencyQC8::bookHistograms() at %s\n", asctime(localtime(&rawTime)));
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

    printf("End of FastEfficiencyQC8::bookHistograms() at %s\n", asctime(localtime(&rawTime)));
}

const GEMGeometry* FastEfficiencyQC8::initGeometry(edm::EventSetup const & iSetup)
{
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

FastEfficiencyQC8::~FastEfficiencyQC8() {
    printf("\n\nFast Efficiency Successfully Done!\n\n");
}

void FastEfficiencyQC8::analyze(const edm::Event& e, const edm::EventSetup& iSetup)
{
    run = e.id().run();
    lumi = e.id().luminosityBlock();
    nev = e.id().event();

    theService->update(iSetup);

    edm::Handle<GEMRecHitCollection> gemRecHits;
    e.getByToken(this->InputTagToken_RH, gemRecHits);
    if (!gemRecHits.isValid())
    {
        edm::LogError("gemcrValidation") << "Cannot get strips by Token RecHits Token.\n";
        return ;
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

    // Array: have a valid event
    bool validEvent[30];

    for (int ch=0; ch<30; ch++)
    {
        validEvent[ch] = true;

        for (unsigned int i = 0; i < beginTripEvt[ch].size(); i++)
        {
            if (beginTripEvt[ch].at(i) <= nev && nev <= endTripEvt[ch].at(i))
            {
                if (ch%2 == 0)
                {
                    validEvent[ch] = false;
                    validEvent[ch+1] = false;
                }
                if (ch%2 == 1)
                {
                    validEvent[ch] = false;
                    validEvent[ch-1] = false;
                }
            }
        }
    }

    // Array of vectors: recHits positions per chamber
    vector<float> xRecHitRef[30];
    vector<int> iEtaRecHitRef[30];
    vector<float> xRecHitTest[30];
    vector<int> iEtaRecHitTest[30];
    vector<int> FirstStripRecHitTest[30];
    vector<int> ClusterSizeRecHitTest[30];

    // recHits loop
    for (GEMRecHitCollection::const_iterator rechit = gemRecHits->begin(); rechit != gemRecHits->end(); ++rechit)
    {
        // calculation of chamber id
        GEMDetId hitID((*rechit).rawId());
        int chIdRecHit = hitID.chamberId().chamber() + hitID.chamberId().layer() - 2;

        if (validEvent[chIdRecHit])
        {
            // cluster size plot and selection
            clusterSize->Fill(chIdRecHit,hitID.roll()-1,(*rechit).clusterSize());
            if ((*rechit).clusterSize()<minCLS) continue;
            if ((*rechit).clusterSize()>maxCLS) continue;

            // recHits plots
            GlobalPoint recHitGP = GEMGeometry_->idToDet((*rechit).gemId())->surface().toGlobal(rechit->localPosition());
            recHits3D->Fill(recHitGP.x(),recHitGP.y(),recHitGP.z());
            recHits2DPerLayer->Fill(recHitGP.x(),hitID.roll()-1,chIdRecHit % 10);
            nRecHitsPerEvtPerCh->Fill(nev,chIdRecHit,hitID.roll()-1);

            // Filling the details of the test recHit per chamber
            xRecHitTest[chIdRecHit].push_back(recHitGP.x());
            iEtaRecHitTest[chIdRecHit].push_back(hitID.roll());
            FirstStripRecHitTest[chIdRecHit].push_back(rechit->firstClusterStrip());
            ClusterSizeRecHitTest[chIdRecHit].push_back(rechit->clusterSize());

            // Reference hits only if in fiducial area

            // Find which chamber is the one in which we had the hit
            int index=-1;
            for (int c=0 ; c<n_ch ; c++)
            {
                if ((gemChambers[c].id().chamber() + gemChambers[c].id().layer() - 2) == chIdRecHit) index = c;
            }
            GEMChamber chHit = gemChambers[index];

            // Define region 'inside' the ieta of the chamber
            int n_strip = chHit.etaPartition(hitID.roll())->nstrips();
            double min_x = chHit.etaPartition(hitID.roll())->centreOfStrip(0).x();
            double max_x = chHit.etaPartition(hitID.roll())->centreOfStrip(n_strip-1).x();

            LocalPoint recHitLP = rechit->localPosition();
            if (((min_x + 4.5) < recHitLP.x()) and (recHitLP.x() < (max_x - 4.5)))
            {
                xRecHitRef[chIdRecHit].push_back(recHitGP.x());
                iEtaRecHitRef[chIdRecHit].push_back(hitID.roll());
            }
        }
    }

    float DxRecHits = 999.9;
    int DiEtaRecHits = 999;
    int confIndex = 999;
    bool numOK = false;

    // Efficiency check: numerator and denominator

    // Reference: even chambers, Under test: odd chambers
    for (int ch=0; ch<=28; ch+=2)
    {
        if (xRecHitRef[ch].size()>0)
        {
            denominator->Fill(ch+1); // Denominator of corresponding odd chamber +1
            denominatorPerEvt->Fill(nev,ch+1);
        }
        else continue;

        numOK = false;

        for (unsigned int i_ref = 0; i_ref < xRecHitTest[ch].size(); i_ref++)
        {
            if (xRecHitTest[ch+1].size()==0) continue;

            confIndex = 0;
            DxRecHits = xRecHitTest[ch+1].at(0) - xRecHitTest[ch].at(i_ref);
            DiEtaRecHits = iEtaRecHitTest[ch+1].at(0) - iEtaRecHitTest[ch].at(i_ref);

            for (unsigned int i_test = 1; i_test < xRecHitTest[ch+1].size(); i_test++)
            {
                if ( ( fabs(iEtaRecHitTest[ch+1].at(i_test) - iEtaRecHitTest[ch].at(i_ref)) <= fabs(DiEtaRecHits+1) ) and ( abs(xRecHitTest[ch+1].at(i_test) - xRecHitTest[ch].at(i_ref)) < abs(DxRecHits) ) )
                {
                    confIndex = i_test;
                    DxRecHits = xRecHitTest[ch+1].at(i_test) - xRecHitTest[ch].at(i_ref);
                    DiEtaRecHits = iEtaRecHitTest[ch+1].at(i_test) - iEtaRecHitTest[ch].at(i_ref);
                }
            }

            DxCorrespondingRecHits->Fill(DxRecHits);
            DiEtaCorrespondingRecHits->Fill(DiEtaRecHits);

            // Fill occupancy plot
            if (fabs(DxRecHits) <= 6.0 and abs(DiEtaRecHits) <= 1)
            {
                for (int clsize = 0; clsize < ClusterSizeRecHitTest[ch+1].at(confIndex); clsize++)
                {
                    occupancyIfConfirmedHits->Fill(FirstStripRecHitTest[ch+1].at(confIndex)+clsize,iEtaRecHitTest[ch+1].at(confIndex)-1,ch+1); // Fill for test
                }
            }

            for (unsigned int i_test = 0; i_test < xRecHitTest[ch+1].size(); i_test++)
            {
                DxRecHits = xRecHitTest[ch+1].at(i_test) - xRecHitTest[ch].at(i_ref);
                DiEtaRecHits = iEtaRecHitTest[ch+1].at(i_test) - iEtaRecHitTest[ch].at(i_ref);

                if (fabs(DxRecHits) <= 6.0 and abs(DiEtaRecHits) <= 1)
                {
                    numOK = true;
                }
            }

        }

        if (numOK)
        {
            numerator->Fill(ch+1); // Numberator of corresponding odd chamber +1
            numeratorPerEvt->Fill(nev,ch+1);
        }
    }

    DxRecHits = 999.9;
    DiEtaRecHits = 999;
    confIndex = 999;
    numOK = false;

    // Reference: odd chambers, Under test: even chambers
    for (int ch=1; ch<=29; ch+=2)
    {
        if (xRecHitRef[ch].size()>0)
        {
            denominator->Fill(ch-1); // Denominator of corresponding even chamber -1
            denominatorPerEvt->Fill(nev,ch-1);
        }
        else continue;

        numOK = false;

        for (unsigned int i_ref = 0; i_ref < xRecHitTest[ch].size(); i_ref++)
        {
            if (xRecHitTest[ch-1].size()==0) continue;

            confIndex = 0;
            DxRecHits = xRecHitTest[ch-1].at(0) - xRecHitTest[ch].at(i_ref);
            DiEtaRecHits = iEtaRecHitTest[ch-1].at(0) - iEtaRecHitTest[ch].at(i_ref);

            for (unsigned int i_test = 1; i_test < xRecHitTest[ch-1].size(); i_test++)
            {
                if ( ( fabs(iEtaRecHitTest[ch-1].at(i_test) - iEtaRecHitTest[ch].at(i_ref)) <= fabs(DiEtaRecHits+1) ) and ( abs(xRecHitTest[ch-1].at(i_test) - xRecHitTest[ch].at(i_ref)) < abs(DxRecHits) ) )
                {
                    confIndex = i_test;
                    DxRecHits = xRecHitTest[ch-1].at(i_test) - xRecHitTest[ch].at(i_ref);
                    DiEtaRecHits = iEtaRecHitTest[ch-1].at(i_test) - iEtaRecHitTest[ch].at(i_ref);
                }
            }

            DxCorrespondingRecHits->Fill(-DxRecHits);
            DiEtaCorrespondingRecHits->Fill(-DiEtaRecHits);

            if (fabs(DxRecHits) <= 6.0 and abs(DiEtaRecHits) <= 1)
            {
                for (int clsize = 0; clsize < ClusterSizeRecHitTest[ch-1].at(confIndex); clsize++)
                {
                    occupancyIfConfirmedHits->Fill(FirstStripRecHitTest[ch-1].at(confIndex)+clsize,iEtaRecHitTest[ch-1].at(confIndex)-1,ch-1); // Fill for test
                }
            }

            for (unsigned int i_test = 0; i_test < xRecHitTest[ch-1].size(); i_test++)
            {
                DxRecHits = xRecHitTest[ch-1].at(i_test) - xRecHitTest[ch].at(i_ref);
                DiEtaRecHits = iEtaRecHitTest[ch-1].at(i_test) - iEtaRecHitTest[ch].at(i_ref);

                if (fabs(DxRecHits) <= 6.0 and abs(DiEtaRecHits) <= 1)
                {
                    numOK = true;
                }
            }
        }

        if (numOK)
        {
            numerator->Fill(ch-1); // Numberator of corresponding even chamber -1
            numeratorPerEvt->Fill(nev,ch-1);
        }
    }

    tree->Fill();
}

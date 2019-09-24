#include "Analysis/GEMQC8/interface/HighGranFastEffQC8.h"

using namespace std;
using namespace edm;

HighGranFastEffQC8::HighGranFastEffQC8(const edm::ParameterSet& cfg): GEMBaseValidation(cfg)
{
	time_t rawTime;
	time(&rawTime);
	printf("Begin of HighGranFastEffQC8::HighGranFastEffQC8() at %s\n", asctime(localtime(&rawTime)));
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

	// Numerator and denominator
	num_block   = fs->make<TH3D>("num_block","num_block",48,0,48,10,0,10,30,0,30);
	denom_block = fs->make<TH3D>("denom_block","denom_block",48,0,48,10,0,10,30,0,30);

	num_block_x   = fs->make<TH3D>("num_block_x","num_block_x",500,-50,50,10,0,10,30,0,30);
	denom_block_x = fs->make<TH3D>("denom_block_x","denom_block_x",500,-50,50,10,0,10,30,0,30);

	// Tree branches declaration
	tree = fs->make<TTree>("tree", "Tree for QC8");
	tree->Branch("run",&run,"run/I");
	tree->Branch("lumi",&lumi,"lumi/I");
	tree->Branch("ev",&nev,"ev/I");

	printf("End of HighGranFastEffQC8::HighGranFastEffQC8() at %s\n", asctime(localtime(&rawTime)));
}

void HighGranFastEffQC8::bookHistograms(DQMStore::IBooker & ibooker, edm::Run const & Run, edm::EventSetup const & iSetup )
{
	time_t rawTime;
	time(&rawTime);
	printf("Begin of HighGranFastEffQC8::bookHistograms() at %s\n", asctime(localtime(&rawTime)));
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

	printf("End of HighGranFastEffQC8::bookHistograms() at %s\n", asctime(localtime(&rawTime)));
}

const GEMGeometry* HighGranFastEffQC8::initGeometry(edm::EventSetup const & iSetup)
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

HighGranFastEffQC8::~HighGranFastEffQC8() {
	printf("\n\nFast Efficiency Successfully Done!\n\n");
}

void HighGranFastEffQC8::analyze(const edm::Event& e, const edm::EventSetup& iSetup)
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

	// Arrays: have a chamber fired
	bool fired_ch_test[30];
	bool fired_ch_reference[30];

	bool validEvent[30];

	for (int ch=0; ch<30; ch++)
	{
		fired_ch_test[ch] = false;
		fired_ch_reference[ch] = false;
		validEvent[ch] = true;
	}

	for (int ch=0; ch<30; ch++)
	{
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
	vector<float> xRecHit[30];
	vector<int> iEtaRecHit[30];
	float DxRecHits = 999.9;
	int DiEtaRecHits = 999.9;
	int nblock = 8;
	double xmin = 0., xmax = 0., xmin_upper = 0., xmax_upper = 0., pitch = 0.;
	int ETA = 0, ETA_upper = 0;
	double strip_x_min = 0., strip_x_max = 0., strip_x_min_upper = 0., strip_x_max_upper = 0.;
	int block_min = 0, block_max = 0; //defind as integer in order to round down
	int nx = 0, nx_upper = 0;

	// recHits loop
	for (GEMRecHitCollection::const_iterator rechit = gemRecHits->begin(); rechit != gemRecHits->end(); ++rechit)
	{
		// calculation of chamber id
		GEMDetId hitID((*rechit).rawId());
		int chIdRecHit = hitID.chamberId().chamber() + hitID.chamberId().layer() - 2;

		if (validEvent[chIdRecHit])
		{
			// cluster size plot and selection
			if ((*rechit).clusterSize()<minCLS) continue;
			if ((*rechit).clusterSize()>maxCLS) continue;

			// fired chambers
			fired_ch_test[chIdRecHit] = true;

			// local point of the recHit
			LocalPoint recHitLP = rechit->localPosition();

			// Filling the details of the recHit per chamber
			xRecHit[chIdRecHit].push_back(recHitLP.x());
			iEtaRecHit[chIdRecHit].push_back(hitID.roll());

			// reference hits only if in fiducial area

			// Find which chamber is the one in which we had the hit
			int index=-1;
	  		for (int c=0 ; c<n_ch ; c++)
	  		{
	    		if ((gemChambers[c].id().chamber() + gemChambers[c].id().layer() - 2) == chIdRecHit)
				{
					index = c;
				}
	  		}

			GEMChamber ch = gemChambers[index];

			// Define region 'inside' the ieta of the chamber
			int n_strip = ch.etaPartition(hitID.roll())->nstrips();
			double min_x = ch.etaPartition(hitID.roll())->centreOfStrip(0).x();
			double max_x = ch.etaPartition(hitID.roll())->centreOfStrip(383).x(); 

			// Get information on strips and blocks in the region defined by +/-6cm
			ETA = hitID.roll();
			ETA_upper = hitID.roll() + 1;

			xmin = fabs(ch.etaPartition(ETA)->centreOfStrip(0).x());
			xmax = fabs(ch.etaPartition(ETA)->centreOfStrip(383).x());
			pitch = (xmin + xmax)/n_strip;

			//fired strip
			int fstrip = rechit->firstClusterStrip()+rechit->clusterSize()/2;

			//Min and max strips in that region
			int region_stripmin = fstrip - 6.0/pitch;
			int region_stripmax = fstrip + 6.0/pitch;

			//If strips outside chamber, consider strips at edges			
			if (region_stripmin < 0)   region_stripmin = 0;
			if (region_stripmax > 383) region_stripmax = 383 ;
			if (region_stripmax < region_stripmin) region_stripmax = 0;

			//X position of min and max strips in that region
			strip_x_min = ch.etaPartition(ETA)->centreOfStrip(region_stripmin).x();
			strip_x_max = ch.etaPartition(ETA)->centreOfStrip(region_stripmax).x();

			//X position of min and max strips also for the upper eta 
			if (ETA_upper>1 && ETA_upper<=8)
			{
				xmin_upper = ch.etaPartition(ETA_upper)->centreOfStrip(0).x();
				xmax_upper = ch.etaPartition(ETA_upper)->centreOfStrip(383).x();
				strip_x_min_upper = strip_x_min;
				strip_x_max_upper = strip_x_max;
				if (strip_x_min_upper < xmin_upper) strip_x_min_upper = xmin_upper;
				if (strip_x_max_upper > xmax_upper) strip_x_max_upper = xmax_upper;
			}

			nx = 10 * (strip_x_max - strip_x_min);
			nx_upper = 10 * (strip_x_max_upper - strip_x_min_upper);

			//Min and max blocks in that region
			block_min = region_stripmin/nblock;
			block_max = region_stripmax/nblock;

			//Exclude edges
			if (min_x < max_x)
			{
				min_x = min_x + 4.5;
				max_x = max_x - 4.5;	
			}

			if (max_x < min_x)
			{
				min_x = min_x - 4.5;
				max_x = max_x + 4.5;
			}

			if ((min_x < recHitLP.x()) and (recHitLP.x() < max_x) )
			{
				fired_ch_reference[chIdRecHit] = true;

			}
		}
	}

	// Efficiency check: numerator and denominator

	bool numerator_fired = false; // Flag to check if numerator hits pass the requirements in delta(x) and delta(iEta)

	// Reference: even chambers, Under test: odd chambers
	for (int ch=0; ch<=28; ch+=2)
	{
		numerator_fired = false;

		int counter = 0;

		for (int i=0; i<30; i++)
		{
			if (fired_ch_test[i]==true) counter++;
		}

		if (fired_ch_reference[ch] == true) // We see hits in even chamber 
		{

			//Fill a 3D histogram: blocks / Eta / chamber (denominator)
			for (int block = block_min; block <= block_max; block++)
			{
				denom_block->Fill(block, ETA, ch+1);
				if (ETA < 8) denom_block->Fill(block, ETA+1, ch+1);
				if (ETA > 1) denom_block->Fill(block, ETA-1, ch+1);
			}

			//Fill a 3D histogram: X position / Eta  / chamber (denominator)
			for (int x = 0; x <= nx; x++)
			{
				denom_block_x->Fill(strip_x_min+x*0.1, ETA, ch+1);
				if (ETA > 1) denom_block_x->Fill(strip_x_min+x*0.1, ETA-1, ch+1);
			}

			for (int x = 0; x <= nx_upper; x++)
			{
				if (ETA < 8) denom_block_x->Fill(strip_x_min_upper+x*0.1, ETA+1, ch+1);
			}

			if (fired_ch_test[ch+1] == true) // We see hits in the corresponding odd chamber
			{
				for (unsigned int i_ref_ch = 0; i_ref_ch < xRecHit[ch].size(); i_ref_ch++)
				{
					for (unsigned int i_test_ch = 0; i_test_ch < xRecHit[ch+1].size(); i_test_ch++)
					{
						// Calculate delta x and delta iEta
						DxRecHits = xRecHit[ch+1].at(i_test_ch) - xRecHit[ch].at(i_ref_ch); 
						DiEtaRecHits = iEtaRecHit[ch+1].at(i_test_ch) - iEtaRecHit[ch].at(i_ref_ch);

						// Fill occupancy plot
						if (fabs(DxRecHits) <= 6.0 and abs(DiEtaRecHits) <= 1)
						{

							numerator_fired = true;
						}
					}
				}
			}
			if (numerator_fired == true)
			{

				//Fill a 3D histogram: blocks / Eta / chamber (numerator)
				for (int block = block_min; block <= block_max; block++)
				{
					num_block->Fill(block, ETA, ch+1);
					if (ETA < 8) num_block->Fill(block, ETA+1, ch+1);
					if (ETA > 1) num_block->Fill(block, ETA-1, ch+1);
				}
	
				//Fill a 3D histogram: X position / Eta  / chamber (denominator)
				for (int x = 0; x <= nx; x++)
				{
					num_block_x->Fill(strip_x_min+x*0.1, ETA, ch+1);
					if (ETA > 1) num_block_x->Fill(strip_x_min+x*0.1, ETA-1, ch+1);
				}
		
				for (int x = 0; x <= nx_upper; x++)
				{
					if (ETA < 8) num_block_x->Fill(strip_x_min_upper+x*0.1, ETA+1, ch+1);
				}
			}
		}
	}


	int counter = 0;

	for (int i=0; i<30; i++)
	{
		if (fired_ch_test[i]==true) counter++;
	}

	// Reference: odd chambers, Under test: even chambers
	for (int ch=1; ch<=29; ch+=2)
	{
		counter = 0;

		for (int i=0; i<30; i++)
		{
			if (fired_ch_test[i]==true) counter++;
		}

		numerator_fired = false;

		if (fired_ch_reference[ch] == true) // We see hits in even chamber
		{
			//Fill a 3D histogram: blocks / Eta / chamber (denominator)
			for (int block = block_min; block <= block_max; block++)
			{
				denom_block->Fill(block, ETA, ch-1);
				if (ETA < 8) denom_block->Fill(block, ETA+1, ch-1);
				if (ETA > 1) denom_block->Fill(block, ETA-1, ch-1);
			}

			//Fill a 3D histogram: X position / Eta  / chamber (denominator)
			for (int x = 0; x <= nx; x++)
			{
				denom_block_x->Fill(strip_x_min+x*0.1, ETA, ch-1);
				if (ETA > 1) denom_block_x->Fill(strip_x_min+x*0.1, ETA-1, ch-1);
			}

			for (int x = 0; x <= nx_upper; x++)
			{
				if (ETA < 8) denom_block_x->Fill(strip_x_min_upper+x*0.1, ETA+1, ch-1);
			}

			if (fired_ch_test[ch-1] == true) // We see hits in the corresponding odd chamber
			{
				for (unsigned int i_ref_ch = 0; i_ref_ch < xRecHit[ch].size(); i_ref_ch++)
				{
					for (unsigned int i_test_ch = 0; i_test_ch < xRecHit[ch-1].size(); i_test_ch++)
					{
						// Calculate delta x and delta iEta
						DxRecHits = xRecHit[ch-1].at(i_test_ch) - xRecHit[ch].at(i_ref_ch);
						DiEtaRecHits = iEtaRecHit[ch-1].at(i_test_ch) - iEtaRecHit[ch].at(i_ref_ch);

						// Fill occupancy plot
						if (fabs(DxRecHits) <= 6.0 and abs(DiEtaRecHits) <= 1)
						{

							numerator_fired = true;
						}
					}
				}
			}

			if (numerator_fired == true)
			{
				//Fill a 3D histogram: blocks / Eta / chamber (numerator)
				for (int block = block_min; block <= block_max; block++)
				{
					num_block->Fill(block, ETA, ch-1);
					if (ETA < 8) num_block->Fill(block, ETA+1, ch-1);
					if (ETA > 1) num_block->Fill(block, ETA-1, ch-1);
				}

				//Fill a 3D histogram: X position / Eta  / chamber (denominator)
				for (int x = 0; x <= nx; x++)
				{
					num_block_x->Fill(strip_x_min+x*0.1, ETA, ch-1);
					if (ETA > 1) num_block_x->Fill(strip_x_min+x*0.1, ETA-1, ch-1);
				}
		
				for (int x = 0; x <= nx_upper; x++)
				{
					if (ETA < 8) num_block_x->Fill(strip_x_min_upper+x*0.1, ETA+1, ch-1);
				}
			}
		}
	}

	tree->Fill();

}

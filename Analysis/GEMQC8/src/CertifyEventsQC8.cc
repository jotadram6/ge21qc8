#include "Analysis/GEMQC8/interface/CertifyEventsQC8.h"

using namespace std;
using namespace edm;

CertifyEventsQC8::CertifyEventsQC8(const edm::ParameterSet& cfg): GEMBaseValidation(cfg)
{
	time_t rawTime;
	time(&rawTime);
	printf("Begin of CertifyEventsQC8::CertifyEventsQC8() at %s\n", asctime(localtime(&rawTime)));
	InputTagToken_RH = consumes<GEMRecHitCollection>(cfg.getParameter<edm::InputTag>("recHitsInputLabel"));
	InputTagToken_DG = consumes<GEMDigiCollection>(cfg.getParameter<edm::InputTag>("gemDigiLabel"));
	edm::ParameterSet serviceParameters = cfg.getParameter<edm::ParameterSet>("ServiceParameters");
	theService = new MuonServiceProxy(serviceParameters);
	minCLS = cfg.getParameter<double>("minClusterSize");
	maxCLS = cfg.getParameter<double>("maxClusterSize");
	theUpdator = new KFUpdator();
	time(&rawTime);

	edm::Service<TFileService> fs;

	// Histograms declaration

	// nRecHitsPerEventPerChamber: evolution in time
	nRecHitsPerEvtPerCh = fs->make<TH2D>("nRecHitsPerEvtPerCh","recHits per ch vs event (packages of 1000 evts)",24000,0,24000000,30,0,30);
	nDigisPerEvtPerCh = fs->make<TH2D>("nDigisPerEvtPerCh","digis per ch vs event (packages of 1000 evts)",24000,0,24000000,30,0,30);

	// Tree branches declaration

	tree = fs->make<TTree>("tree", "Tree for QC8");
	tree->Branch("run",&run,"run/I");
	tree->Branch("lumi",&lumi,"lumi/I");
	tree->Branch("ev",&nev,"ev/I");

	printf("End of CertifyEventsQC8::CertifyEventsQC8() at %s\n", asctime(localtime(&rawTime)));
}

void CertifyEventsQC8::bookHistograms(DQMStore::IBooker & ibooker, edm::Run const & Run, edm::EventSetup const & iSetup )
{
	time_t rawTime;
	time(&rawTime);
	printf("Begin of CertifyEventsQC8::bookHistograms() at %s\n", asctime(localtime(&rawTime)));
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

	printf("End of CertifyEventsQC8::bookHistograms() at %s\n", asctime(localtime(&rawTime)));
}

const GEMGeometry* CertifyEventsQC8::initGeometry(edm::EventSetup const & iSetup)
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

CertifyEventsQC8::~CertifyEventsQC8() {
	printf("\n\nCertify Events Successfully Done!\n\n");
}

void CertifyEventsQC8::analyze(const edm::Event& e, const edm::EventSetup& iSetup)
{
	run = e.id().run();
	lumi = e.id().luminosityBlock();
	nev = e.id().event();

	theService->update(iSetup);

	edm::Handle<GEMDigiCollection> digis;
	e.getByToken( this->InputTagToken_DG, digis);
	for (GEMDigiCollection::DigiRangeIterator gemdgIt = digis->begin(); gemdgIt != digis->end(); ++gemdgIt)
	{
		const GEMDetId& gemId = (*gemdgIt).first;
		int chIdRecHit = gemId.chamberId().chamber()+gemId.chamberId().layer()-2;

		const GEMDigiCollection::Range& range = (*gemdgIt).second;
		for ( auto digi = range.first; digi != range.second; ++digi )
		{
			nDigisPerEvtPerCh->Fill(nev,chIdRecHit);
		}
	}

	edm::Handle<GEMRecHitCollection> gemRecHits;
	e.getByToken(this->InputTagToken_RH, gemRecHits);
	if (!gemRecHits.isValid())
	{
		edm::LogError("gemcrValidation") << "Cannot get strips by Token RecHits Token.\n";
		return ;
	}

	// recHits loop
	for (GEMRecHitCollection::const_iterator rechit = gemRecHits->begin(); rechit != gemRecHits->end(); ++rechit)
	{
		// calculation of chamber id
		GEMDetId hitID((*rechit).rawId());
		int chIdRecHit = hitID.chamberId().chamber() + hitID.chamberId().layer() - 2;

		// cluster size plot and selection
		if ((*rechit).clusterSize()<minCLS) continue;
		if ((*rechit).clusterSize()>maxCLS) continue;

		// recHits plots
		nRecHitsPerEvtPerCh->Fill(nev,chIdRecHit);
	}

	tree->Fill();
}

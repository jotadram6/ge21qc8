#include "Analysis/GEMQC8/interface/HotDeadStripsQC8.h"

using namespace std;
using namespace edm;

HotDeadStripsQC8::HotDeadStripsQC8(const edm::ParameterSet& cfg): GEMBaseValidation(cfg)
{
  time_t rawTime;
  time(&rawTime);
  printf("Begin of HotDeadStripsQC8::HotDeadStripsQC8() at %s\n", asctime(localtime(&rawTime)));
  InputTagToken_DG = consumes<GEMDigiCollection>(cfg.getParameter<edm::InputTag>("gemDigiLabel"));
  edm::ParameterSet serviceParameters = cfg.getParameter<edm::ParameterSet>("ServiceParameters");
  theService = new MuonServiceProxy(serviceParameters);
  theUpdator = new KFUpdator();
  time(&rawTime);

  edm::Service<TFileService> fs;

  // Histograms declaration

  digiStrips = fs->make<TH3D>("digiStrips","digi per strip",384,0,384,8,0,8,30,0,30);

  // Tree branches declaration

  tree = fs->make<TTree>("tree", "Tree for QC8");
  tree->Branch("run",&run,"run/I");
  tree->Branch("lumi",&lumi,"lumi/I");
  tree->Branch("ev",&nev,"ev/I");

  printf("End of HotDeadStripsQC8::HotDeadStripsQC8() at %s\n", asctime(localtime(&rawTime)));
}

void HotDeadStripsQC8::bookHistograms(DQMStore::IBooker & ibooker, edm::Run const & Run, edm::EventSetup const & iSetup ) {
  time_t rawTime;
  time(&rawTime);
  printf("Begin of HotDeadStripsQC8::bookHistograms() at %s\n", asctime(localtime(&rawTime)));
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

  printf("End of HotDeadStripsQC8::bookHistograms() at %s\n", asctime(localtime(&rawTime)));
}

const GEMGeometry* HotDeadStripsQC8::initGeometry(edm::EventSetup const & iSetup) {
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

HotDeadStripsQC8::~HotDeadStripsQC8() {
  printf("Done!");
}

void HotDeadStripsQC8::analyze(const edm::Event& e, const edm::EventSetup& iSetup){

  run = e.id().run();
  lumi = e.id().luminosityBlock();
  nev = e.id().event();

  theService->update(iSetup);

  // digis

  edm::Handle<GEMDigiCollection> digis;
  e.getByToken( this->InputTagToken_DG, digis);
  int nNumDigi = 0;
  for (GEMDigiCollection::DigiRangeIterator gemdgIt = digis->begin(); gemdgIt != digis->end(); ++gemdgIt)
  {
    nNumDigi++;
    const GEMDetId& gemId = (*gemdgIt).first;
    const GEMDigiCollection::Range& range = (*gemdgIt).second;
    for ( auto digi = range.first; digi != range.second; ++digi )
    {
      digiStrips->Fill(digi->strip(),gemId.roll()-1,(gemId.chamberId().chamber()+gemId.chamberId().layer()-2)); // Strip#=[0,383] -> OK , Eta#=[1,8] -> -1
    }
  }

  tree->Fill();
}

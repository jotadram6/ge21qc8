#ifndef FastEfficiencyQC8_H
#define FastEfficiencyQC8_H

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/GeometrySurface/interface/Bounds.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "DataFormats/Provenance/interface/Timestamp.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectoryState/interface/PTrajectoryStateOnDet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Geometry/GEMGeometry/interface/GEMChamber.h"
#include "Geometry/GEMGeometry/interface/GEMSuperChamber.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoMuon/CosmicMuonProducer/interface/CosmicMuonSmoother.h"
#include "RecoMuon/StandAloneTrackFinder/interface/StandAloneMuonSmoother.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "Validation/MuonGEMHits/interface/GEMBaseValidation.h"
#include <DataFormats/GEMRecHit/interface/GEMRecHit.h>
#include <DataFormats/GEMRecHit/interface/GEMRecHitCollection.h>
#include <DataFormats/TrackingRecHit/interface/TrackingRecHit.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <Geometry/GEMGeometry/interface/GEMGeometry.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>

#include <algorithm>
#include <iomanip>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>

class FastEfficiencyQC8 : public GEMBaseValidation
{
public:

  explicit FastEfficiencyQC8( const edm::ParameterSet& );
  ~FastEfficiencyQC8();
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  void analyze(const edm::Event& e, const edm::EventSetup&) override;
  const GEMGeometry* initGeometry(edm::EventSetup const & iSetup);
  const GEMGeometry* GEMGeometry_;
  std::vector<GEMChamber> gemChambers;
  int n_ch;
  double maxCLS, minCLS;
  std::vector<std::string> TripEventsPerCh;
  MuonServiceProxy* theService;
  CosmicMuonSmoother* theSmoother;
  KFUpdator* theUpdator;
  edm::EDGetToken InputTagToken_, InputTagToken_RH;

private:

	TH3D *recHits3D;
  TH3D *recHits2DPerLayer;
  TH3D *clusterSize;
	TH3D *occupancyIfConfirmedHits;
	TH1D *DxCorrespondingRecHits;
	TH1D *DiEtaCorrespondingRecHits;

	TH3D *nRecHitsPerEvtPerCh;

  TH1D *numerator;
  TH1D *denominator;

	TH2D *numeratorPerEvt;
	TH2D *denominatorPerEvt;

  TTree *tree;
  int run;
  int lumi;
  int nev;
};

#endif

#ifndef AlignmentQC8_H
#define AlignmentQC8_H

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
#include "DataFormats/GeometryVector/interface/Point3DBase.h"
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
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include <Geometry/GEMGeometry/interface/GEMGeometry.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>

#include <algorithm>
#include <iomanip>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>

class AlignmentQC8 : public GEMBaseValidation
{
public:
	explicit AlignmentQC8( const edm::ParameterSet& );
	~AlignmentQC8();
	void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
	void analyze(const edm::Event& e, const edm::EventSetup&) override;
	int findIndex(GEMDetId id_);
	const GEMGeometry* initGeometry(edm::EventSetup const & iSetup);
	double maxCLS, minCLS, maxRes, trackChi2, trackResY, trackResX, MulSigmaOnWindow;
	std::vector<std::string> SuperChamType;
	std::vector<double> vecChamType;
	bool makeTrack, isMC;

  const GEMGeometry* GEMGeometry_;

	std::vector<GEMChamber> gemChambers;
	int n_ch;
	MuonServiceProxy* theService;
	CosmicMuonSmoother* theSmoother;
	KFUpdator* theUpdator;
	edm::EDGetToken InputTagToken_, InputTagToken_RH, InputTagToken_TR, InputTagToken_TS, InputTagToken_TJ, InputTagToken_TI, InputTagToken_TT, InputTagToken_DG, InputTagToken_US;

private:

	TH1D *goodVStriggeredEvts;
	TH1D *h_resX_eta[15][8]; // 15 SC max, 8 ieta

	TTree *tree;
  int run;
  int lumi;
  int nev;
	float trajPhi;
  float trajTheta;
  float trajX;
  float trajY;
  float trajZ;
  float trajPx;
  float trajPy;
  float trajPz;
  int nRecHitsTraj;
  float chi2Traj;
  int ndofTraj;
	float testTrajHitX[30];
  float testTrajHitY[30];
  float testTrajHitZ[30];
  float testTrajHitXerr[30];
  float testTrajHitYerr[30];
  float confTestHitX[30];
  float confTestHitY[30];
  float confTestHitZ[30];
  float confTestHitXerr[30];
  float confTestHitYerr[30];
  int confTestHitClSize[30];
  int confTestHitiEta[30];
};

#endif

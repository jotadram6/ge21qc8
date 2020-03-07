#ifndef ValidationQC8_H
#define ValidationQC8_H

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

#include "RecoMuon/CosmicMuonProducer/interface/HeaderForQC8.h"

#include <algorithm>
#include <iomanip>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>


class ValidationQC8 : public GEMBaseValidation
{
public:
  explicit ValidationQC8( const edm::ParameterSet& );
  ~ValidationQC8();
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  void analyze(const edm::Event& e, const edm::EventSetup&) override;
  int findIndex(GEMDetId id_);
  int findVFAT(float x, float a, float b);
  const GEMGeometry* initGeometry(edm::EventSetup const & iSetup);
  double maxCLS, minCLS, maxRes, trackChi2, trackResY, trackResX, MulSigmaOnWindow;
  std::vector<std::string> SuperChamType;
  std::vector<double> vecChamType;
  std::vector<std::string> TripEventsPerCh;
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
  TH3D *hitsVFATnum;
  TH3D *hitsVFATdenom;
  TH3D *digiStrips;
  TH2D *digisPerEvtPerCh;
  TH3D *recHits3D;
  TH3D *recHits2DPerLayer;
  TH3D *associatedHits2DPerLayer;
  TH3D *nonAssociatedHits2DPerLayer;
  TH3D *num2DPerLayer;
  TH3D *denom2DPerLayer;
  TH1D *recHitsPerEvt;
  TH2D *nonAssRecHitsPerEvt;
  TH3D *clusterSize;
  TH3D *associatedHitsClusterSize;
  TH3D *nonAssociatedHitsClusterSize;
  TH1D *residualPhi;
  TH1D *residualEta;
  TH1D *recHitsPerTrack;
  TH1D *genMuAngX;
  TH1D *genMuAngY;
  TH1D *trajMuAngX;
  TH1D *trajMuAngY;
  TH1D *deltaMuAngX;
  TH1D *deltaMuAngY;

  TTree *tree;
  int run;
  int lumi;
  int nev;
  double t_begin;
  double t_end;
  int nDigisPerCh[30];
  int nOfNonAssHits[30];
  int nrecHit;
  int nTraj;
  float trajTheta;
  float trajPhi;
  float trajX;
  float trajY;
  float trajZ;
  float trajPx;
  float trajPy;
  float trajPz;
  float testTrajHitX[30];
  float testTrajHitY[30];
  float testTrajHitZ[30];
  float testTrajHitXerr[30];
  float testTrajHitYerr[30];
  float testTrajHitZerr[30];
  float confTestHitX[30];
  float confTestHitY[30];
  float confTestHitZ[30];
  int nTrajHit; // number of trajHits
  int nTrajRecHit; // number of confirmed trajHits
  float trajAngX;
  float trajAngY;

  TTree *genTree;
  float genMuPx;
  float genMuPy;
  float genMuPz;
  float genMuPt;
  float genMuTheta;
  float genMuPhi;
  float genMuX;
  float genMuY;
  float genMuZ;
  float genAngX;
  float genAngY;
};

#endif

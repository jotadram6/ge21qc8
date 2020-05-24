#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/GEMGeometry/interface/GEMGeometry.h>

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "RecoMuon/CosmicMuonProducer/interface/CosmicMuonSmoother.h"
#include "RecoMuon/StandAloneTrackFinder/interface/StandAloneMuonSmoother.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"

using namespace std;

class AlignmentTrackRecoQC8 : public edm::stream::EDProducer<> {
public:
  /// Constructor
  explicit AlignmentTrackRecoQC8(const edm::ParameterSet&);
  /// Destructor
  virtual ~AlignmentTrackRecoQC8() {}
  /// Produce the GEMSegment collection
  void produce(edm::Event&, const edm::EventSetup&) override;
  double maxCLS;
  double minCLS;
  double trackChi2, trackResX, trackResY;
  double MulSigmaOnWindow;
  unsigned int minRecHitsPerTrack;
  std::vector<std::string> g_SuperChamType;
  vector<double> g_vecChamType;
  std::vector<std::string> TripEventsPerCh;
private:
  int iev; // events through
  edm::EDGetTokenT<GEMRecHitCollection> theGEMRecHitToken;
  CosmicMuonSmoother* theSmoother;
  MuonServiceProxy* theService;
  KFUpdator* theUpdator;
  int findSeeds(std::vector<TrajectorySeed> *tmptrajectorySeeds,
                MuonTransientTrackingRecHit::MuonRecHitContainer &seedupRecHits,
                MuonTransientTrackingRecHit::MuonRecHitContainer &seeddnRecHits);
  Trajectory makeTrajectory(TrajectorySeed seed, MuonTransientTrackingRecHit::MuonRecHitContainer &muRecHits, std::vector<GEMChamber> gemChambers);
  const GEMGeometry* gemGeom;
  int nev;
};

AlignmentTrackRecoQC8::AlignmentTrackRecoQC8(const edm::ParameterSet& ps) : iev(0) {
  time_t rawTime;
  time(&rawTime);
  printf("Begin of AlignmentTrackRecoQC8::AlignmentTrackRecoQC8() at %s\n", asctime(localtime(&rawTime)));
  maxCLS = ps.getParameter<double>("maxClusterSize");
  minCLS = ps.getParameter<double>("minClusterSize");
  trackChi2 = ps.getParameter<double>("trackChi2");
  trackResX = ps.getParameter<double>("trackResX");
  trackResY = ps.getParameter<double>("trackResY");
  MulSigmaOnWindow = ps.getParameter<double>("MulSigmaOnWindow");
  minRecHitsPerTrack = ps.getParameter<unsigned int>("minNumberOfRecHitsPerTrack");
  g_SuperChamType = ps.getParameter<vector<string>>("SuperChamberType");
  g_vecChamType = ps.getParameter<vector<double>>("SuperChamberSeedingLayers");
  theGEMRecHitToken = consumes<GEMRecHitCollection>(ps.getParameter<edm::InputTag>("gemRecHitLabel"));
  // register what this produces
  edm::ParameterSet serviceParameters = ps.getParameter<edm::ParameterSet>("ServiceParameters");
  theService = new MuonServiceProxy(serviceParameters);
  edm::ParameterSet smootherPSet = ps.getParameter<edm::ParameterSet>("MuonSmootherParameters");
  theSmoother = new CosmicMuonSmoother(smootherPSet,theService);
  theUpdator = new KFUpdator();
  produces<reco::TrackCollection>();
  produces<TrackingRecHitCollection>();
  produces<reco::TrackExtraCollection>();
  produces<vector<Trajectory> >();
  produces<vector<TrajectorySeed> >();
  printf("End of AlignmentTrackRecoQC8::AlignmentTrackRecoQC8() at %s\n", asctime(localtime(&rawTime)));
}


void AlignmentTrackRecoQC8::produce(edm::Event& ev, const edm::EventSetup& setup)
{
  nev = ev.id().event();

  unique_ptr<reco::TrackCollection >          trackCollection( new reco::TrackCollection() );
  unique_ptr<TrackingRecHitCollection >       trackingRecHitCollection( new TrackingRecHitCollection() );
  unique_ptr<reco::TrackExtraCollection >     trackExtraCollection( new reco::TrackExtraCollection() );
  unique_ptr<vector<Trajectory> >             trajectorys( new vector<Trajectory>() );
  unique_ptr<vector<TrajectorySeed> >         trajectorySeeds( new vector<TrajectorySeed>() );
  TrackingRecHitRef::key_type recHitsIndex = 0;
  TrackingRecHitRefProd recHitCollectionRefProd = ev.getRefBeforePut<TrackingRecHitCollection>();
  reco::TrackExtraRef::key_type trackExtraIndex = 0;
  reco::TrackExtraRefProd trackExtraCollectionRefProd = ev.getRefBeforePut<reco::TrackExtraCollection>();

  theService->update(setup);

  edm::ESHandle<GEMGeometry> gemg;
  setup.get<MuonGeometryRecord>().get(gemg);
  const GEMGeometry* mgeom = &*gemg;
  gemGeom = &*gemg;

  vector<GEMChamber> gemChambers;

  const std::vector<const GEMSuperChamber*>& superChambers_ = mgeom->superChambers();
  for (auto sch : superChambers_)
  {
    int n_lay = sch->nChambers();
    for (int l=0;l<n_lay;l++)
    {
      gemChambers.push_back(*sch->chamber(l+1));
    }
  }

  // get the collection of GEMRecHit
  edm::Handle<GEMRecHitCollection> gemRecHits;
  ev.getByToken(theGEMRecHitToken,gemRecHits);

  if (gemRecHits->size() < minRecHitsPerTrack or gemRecHits->size() > 15)
  {
    ev.put(std::move(trajectorySeeds));
    ev.put(std::move(trackCollection));
    ev.put(std::move(trackingRecHitCollection));
    ev.put(std::move(trackExtraCollection));
    ev.put(std::move(trajectorys));
    return;
  }

  MuonTransientTrackingRecHit::MuonRecHitContainer muRecHits;
  MuonTransientTrackingRecHit::MuonRecHitContainer seedupRecHits;
  MuonTransientTrackingRecHit::MuonRecHitContainer seeddnRecHits;

  unsigned int TCN = 0; // Number of hit chambers
  for (auto ch : gemChambers)
  {
    int index = ch.id().chamber() + ch.id().layer() - 2;
    int nHitOnceFilter = 0;
    for (auto etaPart : ch.etaPartitions())
    {
      GEMDetId etaPartID = etaPart->id();
      GEMRecHitCollection::range range = gemRecHits->get(etaPartID);

      for (GEMRecHitCollection::const_iterator rechit = range.first; rechit!=range.second; ++rechit)
      {
        const GeomDet* geomDet(etaPart);
        if ((*rechit).clusterSize()<minCLS) continue;
        if ((*rechit).clusterSize()>maxCLS) continue;
        muRecHits.push_back(MuonTransientTrackingRecHit::specificBuild(geomDet,&*rechit));

        if ( nHitOnceFilter == 0 ) {
          TCN++;
          nHitOnceFilter = 1;
        }

        if ( g_vecChamType[ index ] == 2 || g_vecChamType[ index ] == 4 ) {
          seedupRecHits.push_back(MuonTransientTrackingRecHit::specificBuild(geomDet,&*rechit));
        } else if ( g_vecChamType[ index ] == 1 || g_vecChamType[ index ] == 3 ) {
          seeddnRecHits.push_back(MuonTransientTrackingRecHit::specificBuild(geomDet,&*rechit));
        }
      }
    }
  }
  if (muRecHits.size() < minRecHitsPerTrack || TCN < minRecHitsPerTrack)
  {
    ev.put(std::move(trajectorySeeds));
    ev.put(std::move(trackCollection));
    ev.put(std::move(trackingRecHitCollection));
    ev.put(std::move(trackExtraCollection));
    ev.put(std::move(trajectorys));
    return;
  }

  vector<TrajectorySeed> trajSeedsBody;
  std::vector<TrajectorySeed> *trajSeeds = &trajSeedsBody;
  findSeeds(trajSeeds, seedupRecHits, seeddnRecHits);
  Trajectory bestTrajectory;
  TrajectorySeed bestSeed;

  float maxChi2 = 10000000.0;

  for (auto seed : *trajSeeds)
  {
    Trajectory smoothed = makeTrajectory(seed, muRecHits, gemChambers);

    if (smoothed.isValid())
    {
      float dProbChiNDF = smoothed.chiSquared()/float(smoothed.ndof());

      if ((dProbChiNDF > 0) && (dProbChiNDF < maxChi2))
      {
        maxChi2 = dProbChiNDF;
        bestTrajectory = smoothed;
        bestSeed = seed;
      }
    }
  }

  if (!bestTrajectory.isValid() || maxChi2 > trackChi2)
  {
    ev.put(std::move(trajectorySeeds));
    ev.put(std::move(trackCollection));
    ev.put(std::move(trackingRecHitCollection));
    ev.put(std::move(trackExtraCollection));
    ev.put(std::move(trajectorys));
    return;
  }

  const FreeTrajectoryState* ftsAtVtx = bestTrajectory.geometricalInnermostState().freeState();

  GlobalPoint pca = ftsAtVtx->position();
  math::XYZPoint persistentPCA(pca.x(),pca.y(),pca.z());
  GlobalVector p = ftsAtVtx->momentum();
  math::XYZVector persistentMomentum(p.x(),p.y(),p.z());

  reco::Track track(bestTrajectory.chiSquared(),
                    bestTrajectory.ndof(true),
                    persistentPCA,
                    persistentMomentum,
                    ftsAtVtx->charge(),
                    ftsAtVtx->curvilinearError());

  //sets the outermost and innermost TSOSs
  TrajectoryStateOnSurface outertsos;
  TrajectoryStateOnSurface innertsos;
  unsigned int innerId, outerId;

  // ---  NOTA BENE: the convention is to sort hits and measurements "along the momentum".
  // This is consistent with innermost and outermost labels only for tracks from LHC collision
  if (bestTrajectory.direction() == alongMomentum) {
    outertsos = bestTrajectory.lastMeasurement().updatedState();
    innertsos = bestTrajectory.firstMeasurement().updatedState();
    outerId = bestTrajectory.lastMeasurement().recHit()->geographicalId().rawId();
    innerId = bestTrajectory.firstMeasurement().recHit()->geographicalId().rawId();
  } else {
    outertsos = bestTrajectory.firstMeasurement().updatedState();
    innertsos = bestTrajectory.lastMeasurement().updatedState();
    outerId = bestTrajectory.firstMeasurement().recHit()->geographicalId().rawId();
    innerId = bestTrajectory.lastMeasurement().recHit()->geographicalId().rawId();
  }

  //build the TrackExtra
  GlobalPoint gv = outertsos.globalParameters().position();
  GlobalVector gp = outertsos.globalParameters().momentum();
  math::XYZVector outmom(gp.x(), gp.y(), gp.z());
  math::XYZPoint outpos(gv.x(), gv.y(), gv.z());
  gv = innertsos.globalParameters().position();
  gp = innertsos.globalParameters().momentum();
  math::XYZVector inmom(gp.x(), gp.y(), gp.z());
  math::XYZPoint inpos(gv.x(), gv.y(), gv.z());

  auto tx = reco::TrackExtra(outpos,
                             outmom,
                             true,
                             inpos,
                             inmom,
                             true,
                             outertsos.curvilinearError(),
                             outerId,
                             innertsos.curvilinearError(),
                             innerId,
                             bestSeed.direction(),
                             bestTrajectory.seedRef());

  //adding rec hits
  Trajectory::RecHitContainer transHits = bestTrajectory.recHits();
  unsigned int nHitsAdded = 0;
  for (Trajectory::RecHitContainer::const_iterator recHit = transHits.begin(); recHit != transHits.end(); ++recHit)
  {
    TrackingRecHit *singleHit = (**recHit).hit()->clone();
    trackingRecHitCollection->push_back( singleHit );
    ++nHitsAdded;
  }

  tx.setHits(recHitCollectionRefProd, recHitsIndex, nHitsAdded);
  recHitsIndex += nHitsAdded;

  trackExtraCollection->push_back(tx);

  reco::TrackExtraRef trackExtraRef(trackExtraCollectionRefProd, trackExtraIndex++ );
  track.setExtra(trackExtraRef);
  trackCollection->push_back(track);

  trajectorys->push_back(bestTrajectory);
  trajectorySeeds->push_back(bestSeed);

  // fill the collection
  // put collection in event
  ev.put(std::move(trajectorySeeds));
  ev.put(std::move(trackCollection));
  ev.put(std::move(trackingRecHitCollection));
  ev.put(std::move(trackExtraCollection));
  ev.put(std::move(trajectorys));
}

int AlignmentTrackRecoQC8::findSeeds(std::vector<TrajectorySeed> *tmptrajectorySeeds,
                                   MuonTransientTrackingRecHit::MuonRecHitContainer &seedupRecHits,
                                   MuonTransientTrackingRecHit::MuonRecHitContainer &seeddnRecHits)
{
  for (auto hit1 : seeddnRecHits){
    for (auto hit2 : seedupRecHits){
      if (hit1->globalPosition().z() < hit2->globalPosition().z())
      {
        // Seeds in the same column
        GEMDetId detId1(hit1->rawId()), detId2(hit2->rawId());
  		  uint32_t chDown = detId1.chamber() + detId1.layer() - 2;
  		  uint32_t chUp = detId2.chamber() + detId2.layer() - 2;
	      if (int(chDown/10) != int(chUp/10)) continue;

        LocalPoint segPos = hit1->localPosition();
        GlobalVector segDirGV(hit2->globalPosition().x() - hit1->globalPosition().x(),
                              hit2->globalPosition().y() - hit1->globalPosition().y(),
                              hit2->globalPosition().z() - hit1->globalPosition().z());

        segDirGV *=10;
        LocalVector segDir = hit1->det()->toLocal(segDirGV);

        int charge= 1;
        LocalTrajectoryParameters param(segPos, segDir, charge);

        AlgebraicSymMatrix mat(5,0);
        mat = hit1->parametersError().similarityT( hit1->projectionMatrix() );
        LocalTrajectoryError error(asSMatrix<5>(mat));

        TrajectoryStateOnSurface tsos(param, error, hit1->det()->surface(), &*theService->magneticField());
        uint32_t id = hit1->rawId();
        PTrajectoryStateOnDet const & seedTSOS = trajectoryStateTransform::persistentState(tsos, id);

        edm::OwnVector<TrackingRecHit> seedHits;
        seedHits.push_back(hit1->hit()->clone());
        seedHits.push_back(hit2->hit()->clone());

        TrajectorySeed seed(seedTSOS,seedHits,alongMomentum);

        tmptrajectorySeeds->push_back(seed);
      }
    }
  }
  return 0;
}


Trajectory AlignmentTrackRecoQC8::makeTrajectory(TrajectorySeed seed, MuonTransientTrackingRecHit::MuonRecHitContainer &muRecHits, std::vector<GEMChamber> gemChambers)
{
  PTrajectoryStateOnDet ptsd1(seed.startingState());
  DetId did(ptsd1.detId());
  const BoundPlane& bp = theService->trackingGeometry()->idToDet(did)->surface();
  TrajectoryStateOnSurface tsos = trajectoryStateTransform::transientState(ptsd1,&bp,&*theService->magneticField());
  TrajectoryStateOnSurface tsosCurrent = tsos;
  TransientTrackingRecHit::ConstRecHitContainer consRecHits;

  TrajectorySeed::range range = seed.recHits();
  int nseed = 0;
  GlobalPoint seedGP[2];
  for (edm::OwnVector<TrackingRecHit>::const_iterator rechit = range.first; rechit!=range.second; ++rechit){
    GEMDetId hitID(rechit->rawId());
    seedGP[nseed] = gemGeom->idToDet((*rechit).rawId())->surface().toGlobal(rechit->localPosition());
    nseed++;
  }
  std::map<double,int> rAndhit;

  for (auto ch : gemChambers)
  {
    std::shared_ptr<MuonTransientTrackingRecHit> tmpRecHit;
    tsosCurrent = theService->propagator("SteppingHelixPropagatorAny")->propagate(tsosCurrent, theService->trackingGeometry()->idToDet(ch.id())->surface());
    if (!tsosCurrent.isValid()) return Trajectory();
    GlobalPoint tsosGP = tsosCurrent.freeTrajectoryState()->position();

    float maxR = 9999;
    int nhit=-1;
    int tmpNhit=-1;
    double tmpR=-1;
    for (auto hit : muRecHits)
    {
      nhit++;
      GEMDetId hitID(hit->rawId());
      if (hitID.chamberId() == ch.id() )
      {
        GlobalPoint hitGP = hit->globalPosition();
        double y_err = hit->localPositionError().yy();
        if (fabs(hitGP.x() - tsosGP.x()) > trackResX * MulSigmaOnWindow) continue;
        if (fabs(hitGP.y() - tsosGP.y()) > trackResY * MulSigmaOnWindow * y_err) continue;
        float deltaR = (hitGP - tsosGP).mag();
        if (maxR > deltaR)
        {
          tmpRecHit = hit;
          maxR = deltaR;
          tmpNhit = nhit;
          tmpR = (hitGP - seedGP[0]).mag();
        }
      }
    }
    if (tmpRecHit)
    {
      rAndhit[tmpR] = tmpNhit;
    }
  }

  if (rAndhit.size() < minRecHitsPerTrack) return Trajectory();
  vector<pair<double,int>> rAndhitV;
  copy(rAndhit.begin(), rAndhit.end(), back_inserter<vector<pair<double,int>>>(rAndhitV));
  for(unsigned int i=0;i<rAndhitV.size();i++)
  {
    consRecHits.push_back(muRecHits[rAndhitV[i].second]);
  }

  if (consRecHits.size() < minRecHitsPerTrack) return Trajectory();
  vector<Trajectory> fitted = theSmoother->trajectories(seed, consRecHits, tsos);
  if ( fitted.size() <= 0 ) return Trajectory();

  Trajectory smoothed = fitted.front();
  return fitted.front();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(AlignmentTrackRecoQC8);

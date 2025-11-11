/*!
 *  \file HoughTrackFinder.cc
 *  \author Christof Roland
 */

//begin

#include "HoughTrackFinder.h"

//#include "AssocInfoContainer.h"                         // for AssocInfoCont...

// sPHENIX Geant4 includes
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrClusterIterationMapv1.h>
#include <trackbase/TrkrDefs.h>  // for getLayer, clu...

#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeedContainer_v1.h>
#include <trackbase_historic/TrackSeed_v2.h>

// sPHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHTimer.h>                              // for PHTimer
#include <phool/getClass.h>

// ROOT includes for debugging
#include <TFile.h>
#include <TH1.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TMath.h>

#include <Eigen/Core>                  // for Matrix
#include <Eigen/Dense>

// BOOST for combi seeding
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

// standard includes
#include <algorithm>
#include <cfloat>
#include <climits>                                     // for UINT_MAX
#include <cmath>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <set>                                          // for set
#include <string>                                       // for string
#include <tuple>
#include <utility>
#include <vector>                                       // for pair, make_pair

// forward declarations
class BbcVertexMap;

class PHCompositeNode;

class PHG4CellContainer;
class PHG4CylinderGeomContainer;
class PHG4Particle;
class PHG4TruthInfoContainer;
class PHTimer;

class sPHENIXSeedFinder;

class SvtxClusterMap;
class SvtxCluster;
class SvtxTrackMap;
class SvtxTrack;
class SvtxTrackState;
class SvtxVertexMap;
class SvtxHitMap;

class TNtuple;
class TFile;

// Boost geometry aliases
namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

// INTERPRET this point as (r, phi, z) in cylindrical coordinates
using point    = bg::model::point<float, 3, bg::cs::cartesian>;
using box      = bg::model::box<point>;
using pointKey = std::pair<point, TrkrDefs::cluskey>;

using myrtree = bgi::rtree<pointKey, bgi::quadratic<16>>;

#define LogDebug(exp) std::cout << "DEBUG: "   << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogError(exp) std::cout << "ERROR: "   << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp

using namespace std;

namespace
{
  inline float dcaBinWidth(float dca)
  {
    const float fine_edge    = 3.0f;
    const float fine_width   = 0.5f;

    const float mid_edge     = 10.0f;
    const float mid_width    = 2.0f;

    const float coarse_width = 20.0f;

    const float ad = std::fabs(dca);
    if (ad < fine_edge)          return fine_width;
    else if (ad < mid_edge)      return mid_width;
    else                         return coarse_width;
  }

  // NEW: clip [low,high] to the region that |dca| belongs to
  inline void dcaSearchRange(float dca,
                             float width,
                             float dcamin_global,
                             float dcamax_global,
                             float& low,
                             float& high)
  {
    const float e1 = 3.0f;   // boundary 1
    const float e2 = 10.0f;  // boundary 2
    const float ad = std::fabs(dca);

    low  = dca - width;
    high = dca + width;

    // Clip to overall global min/max first (safety)
    low  = std::max(low,  dcamin_global);
    high = std::min(high, dcamax_global);

    if (ad < e1)
    {
      // region [-e1, +e1]
      low  = std::max(low,  -e1);
      high = std::min(high, +e1);
    }
    else if (ad < e2)
    {
      // region (e1,e2] or [-e2,-e1)
      if (dca > 0)
      {
        low  = std::max(low,  +e1);
        high = std::min(high, +e2);
      }
      else
      {
        low  = std::max(low,  -e2);
        high = std::min(high, -e1);
      }
    }
    else
    {
      // outer region |dca| > e2
      if (dca > 0)
      {
        low  = std::max(low, +e2);
        // high already clipped to dcamax_global
      }
      else
      {
        high = std::min(high, -e2);
        // low already clipped to dcamin_global
      }
    }
  }
}

vector<TrkrCluster*> clusterpoints;

HoughTrackFinder::HoughTrackFinder(const std::string& name)
  : SubsysReco(name)
{
}

int HoughTrackFinder::GetNodes(PHCompositeNode* topNode)
{
  tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!tGeometry)
  {
    std::cout << PHWHERE << "No acts reco geometry, bailing."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    std::cout << PHWHERE << "No cluster container on node tree, bailing."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int HoughTrackFinder::Init(PHCompositeNode* topNode)
{
  if(topNode==nullptr){
    std::cout << PHWHERE << "No topNode, bailing."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  _nevent = 0;
  if(_write_ntp){
    _tfile = new TFile("./costuple.root", "RECREATE");
    _ntp_cos = new TNtuple("ntp_cos", "cos event info","ev:x:y:z");
    // store (dca, phi, tan(phi), x0)
    _ntp_stub = new TNtuple("ntp_stub", "cos stub info","ev:dca:phi:tanphi:x0");
    // store (dca min/max, phi min/max)
    _ntp_max = new TNtuple("ntp_max", "cos stub info","ev:dcamin:dcamax:phimin:phimax");
    // full track params: use dca, phi, alpha=phi
    _ntp_trk = new TNtuple("ntp_trk", "full track params",
                           "ev:dca:phi:alpha:dist0:nclus");

    const int   nDcaBins  = 240;
    const float dcaMin    = -61;
    const float dcaMax    =  61;

    const int   nPhiBins  = 180;
    const float phiMin    = -TMath::Pi();
    const float phiMax    =  TMath::Pi();

    _hHough = new TH2F("hHough", "Hough space;DCA (cm);#phi (rad)",
                       nDcaBins, dcaMin, dcaMax,
                       nPhiBins, phiMin, phiMax);

                           _track_tree = new TTree("track_clusters", "Tracks and their clusters");

    _track_tree->Branch("ev",    &_trk_ev,    "ev/I");
    _track_tree->Branch("dca",   &_trk_dca,   "dca/F");
    _track_tree->Branch("phi",   &_trk_phi,   "phi/F");
    _track_tree->Branch("alpha", &_trk_alpha, "alpha/F");
    _track_tree->Branch("dist0", &_trk_dist0, "dist0/F");
    _track_tree->Branch("nclus", &_trk_nclus, "nclus/I");

    _track_tree->Branch("clus_r",   &_trk_r);
    _track_tree->Branch("clus_phi", &_trk_phi_clu);
    _track_tree->Branch("clus_z",   &_trk_z);
    _track_tree->Branch("cluskey",  &_trk_cluskey);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int HoughTrackFinder::InitRun(PHCompositeNode* topNode)
{
  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }
  ret = createNodes(topNode);
  return ret;
}

int HoughTrackFinder::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in HoughTrackFinder::createNodes");
  }

  PHNodeIterator dstIter(dstNode);
  PHCompositeNode* svtxNode = dynamic_cast<PHCompositeNode*>(dstIter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  m_seedContainer = findNode::getClass<TrackSeedContainer>(topNode, m_trackMapName);
  if (!m_seedContainer)
  {
    m_seedContainer = new TrackSeedContainer_v1;
    PHIODataNode<PHObject>* trackNode =
        new PHIODataNode<PHObject>(m_seedContainer, m_trackMapName, "PHObject");
    svtxNode->addNode(trackNode);
  }
  if (m_seedContainer){
    cout << "SEED CONTAINER CREATED" << endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

double HoughTrackFinder::phiadd(double phi1, double phi2){
  double s = phi1+phi2;
  if(s>2*M_PI) {return s-2*M_PI;}
  else if(s<0) {return s+2*M_PI;}
  else {return s;}
}

double HoughTrackFinder::phidiff(double phi1, double phi2){
  double d = phi1-phi2;
  if(d>M_PI) {return d-2*M_PI;}
  else if(d<-M_PI) {return d+2*M_PI;}
  else {return d;}
}

// NOTE: search_rtree stores (r, phi, z); we still fit a line in XY
void HoughTrackFinder::get_stub(const myrtree &search_rtree,
                                float pointr, float pointphi, float pointz,
                                int &count, double &slope, double &intercept)
{
  // search window in (r, phi, z)
  const float m1_dr   = 1.5;    // cm in radius
  const float m1_dphi = 0.05;   // rad
  const float m1_dz   = 2.0;    // cm

  std::vector<pointKey> boxclusters;

  // 3D box: (r±dr, phi±dphi, z±dz)
  search_rtree.query(
      bgi::intersects(
        box(
          point(pointr   - m1_dr,
                pointphi - m1_dphi,
                pointz   - m1_dz),
          point(pointr   + m1_dr,
                pointphi + m1_dphi,
                pointz   + m1_dz)
        )
      ),
      std::back_inserter(boxclusters));

  // convert center back to (x,y) for debug
  float cx = pointr * std::cos(pointphi);
  float cy = pointr * std::sin(pointphi);

  std::cout << "!!!!get_stub  center r=" << pointr
            << " phi=" << pointphi
            << " z=" << pointz
            << " (x=" << cx << ", y=" << cy << ")"
            << " r:[" << pointr-m1_dr   << "," << pointr+m1_dr   << "]"
            << " phi:[" << pointphi-m1_dphi << "," << pointphi+m1_dphi << "]"
            << " z:[" << pointz-m1_dz  << "," << pointz+m1_dz   << "]"
            << " nboxclus " << boxclusters.size() << std::endl;

  // still fit a line in XY
  std::vector<std::pair<double,double>> pts;

  int    nhit  = 0;
  double xsum  = 0;
  double x2sum = 0;
  double ysum  = 0;
  double xysum = 0;

  for (auto pbox = boxclusters.begin(); pbox!=boxclusters.end(); ++pbox)
  {
    // R-tree stores (r, phi, z)
    float r   = pbox->first.get<0>();
    float phi = pbox->first.get<1>();
    // float z = pbox->first.get<2>(); // not used for XY fit

    float boxx = r * std::cos(phi);
    float boxy = r * std::sin(phi);

    nhit++;
    pts.emplace_back(boxx, boxy);

    xsum  += boxx;
    ysum  += boxy;
    x2sum += boxx*boxx;
    xysum += boxx*boxy;
  }

  count = nhit;

  const int MIN_STUB_HITS = 2;
  if (nhit < MIN_STUB_HITS)
  {
    slope = 0;
    intercept = 0;
    return;
  }

  const double denominator = (x2sum * nhit) - (xsum*xsum);
  if (denominator == 0)
  {
    slope = 0;
    intercept = 0;
    count = 0;
    return;
  }

  slope     = (xysum * nhit - xsum * ysum) / denominator;
  intercept = (x2sum * ysum - xsum * xysum) / denominator;

  double chi2 = 0.;
  for (const auto &p : pts)
  {
    const double xx = p.first;
    const double yy = p.second;
    const double yfit = slope*xx + intercept;
    const double dy   = yy - yfit;
    chi2 += dy*dy;
    std::cout << " stub pt x: " << xx
              << " y: " << yy
              << " slope: " << slope
              << " intercept: " << intercept
              << " chi2: " << chi2 << std::endl;
  }

  const double ndof = nhit - 2;
  const double chi2ndf = (ndof > 0) ? chi2/ndof : 1e9;
  std::cout << " chi2ndf: " << chi2ndf << " nhit " << nhit << std::endl;

  // still very tight – tune as needed
  const double MAX_STUB_CHI2NDF = 0.02;

  if (!std::isfinite(slope) || chi2ndf > MAX_STUB_CHI2NDF)
  {
    slope = 0;
    intercept = 0;
    count = 0;
    return;
  }
}

//____________________________________________________________________________..
int HoughTrackFinder::process_event(PHCompositeNode* topNode)
{
  if(topNode==nullptr){
    std::cout << PHWHERE << "No topNode, bailing."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
    
  _nevent++;
  std::vector<TrackSeed_v2> clean_chains;

  // Fill rtree of clusters in CYLINDRICAL (r, phi, z)
  bgi::rtree<pointKey, bgi::quadratic<16> > rtree;

  // Hough-space ranges in (dca, phi)
  float phimin =  99999999999.9;
  float phimax = -99999999999.9;
  float dcamin = 99999999999.9;
  float dcamax =-99999999999.9;

  std::cout << "!!!!!!!!!number of clusters in event: " << _cluster_map->size() << endl;
  if(_cluster_map->size()<_min_nclusters){
    std::cout << " not enough clusters in event: " << _cluster_map->size()
              << " < " << _min_nclusters << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // build cylindrical R-tree: (r, phi, z)
  for(const auto& hitsetkey:_cluster_map->getHitSetKeys(TrkrDefs::TrkrId::tpcId)){
    auto range = _cluster_map->getClusters(hitsetkey);
    for( auto clusIter = range.first; clusIter != range.second; ++clusIter ){
      TrkrDefs::cluskey ckey = clusIter->first;
      TrkrCluster *cluster = clusIter->second;

      // get global position (x,y,z)
      const Acts::Vector3 globalpos_d =  tGeometry->getGlobalPosition(ckey, cluster);
      const float gx = globalpos_d.x();
      const float gy = globalpos_d.y();
      const float gz = globalpos_d.z();

      if(_write_ntp)_ntp_cos->Fill(_nevent, gx, gy, gz);

      // convert to cylindrical coordinates
      const float r   = std::sqrt(gx*gx + gy*gy);
      const float phi = std::atan2(gy, gx);

      // remove near-duplicates in (r,phi,z)
      const float eps_r   = 0.01f;   // cm
      const float eps_phi = 0.001f;  // rad
      const float eps_z   = 0.05f;   // cm

      std::vector<pointKey> testduplicate;
      rtree.query(
        bgi::intersects(
          box(point(r   - eps_r,
                    phi - eps_phi,
                    gz  - eps_z),
              point(r   + eps_r,
                    phi + eps_phi,
                    gz  + eps_z))),
        std::back_inserter(testduplicate));

      if (!testduplicate.empty()){
        std::cout<< " duplicate cluster found at r: " << r
                 << " phi: " << phi
                 << " z: " << gz << std::endl;
        continue;
      }

      // store in R-tree as (r, phi, z)
      rtree.insert(std::make_pair(point(r, phi, gz), ckey));
    }
  }

  // cylindrical selection box in (r, phi, z)
  const float RMIN   = 0.0;
  const float RMAX   = 80.0;
  const float PHIMIN = -TMath::Pi();
  const float PHIMAX =  TMath::Pi();
  const float ZMIN   = -120.0;   // or something big enough for your TPC
  const float ZMAX   =  120.0;

  // Get all clusters from rtree
  vector<pointKey> allclusters;
  rtree.query(
    bgi::intersects(box(point(RMIN,   PHIMIN, ZMIN),
                        point(RMAX,   PHIMAX, ZMAX))),
    std::back_inserter(allclusters));
  cout << "!!!!number clus is " << allclusters.size() << endl;

  // R-tree in Hough space: (dca, phi)
  bgi::rtree<pointKey, bgi::quadratic<16> > rtree_stub;

  for (vector<pointKey>::iterator cluster = allclusters.begin();
       cluster != allclusters.end(); ++cluster)
  {
    // tree stores (r, phi, z)
    float pointr   = cluster->first.get<0>();
    float pointphi = cluster->first.get<1>();
    float pointz   = cluster->first.get<2>();

    float pointx = pointr * std::cos(pointphi);
    float pointy = pointr * std::sin(pointphi);

    int    fcount     = 0;
    double fslope     = 0;
    double fintercept = 0;

    std::cout << " Processing cluster at r=" << pointr
              << " phi=" << pointphi
              << " z=" << pointz
              << " (x=" << pointx << ", y=" << pointy << ")"
              << std::endl;

    // search for local stub around (r,phi,z)
    get_stub(rtree, pointr, pointphi, pointz, fcount, fslope, fintercept);

    if (fcount <= 0) continue;
    if (std::isnan(fslope) || !std::isfinite(fslope)) continue;
    if (fslope ==0 && fintercept ==0) continue;

    if (std::isfinite(fslope))
    {
      // convert (slope, intercept) -> (phi, dca) in XY
      const double phi = std::atan(fslope);
      const double dca = fintercept / std::sqrt(fslope*fslope + 1.0);
      if (dca ==0 && phi ==0) continue;

      // insert stub in Hough space (dca, phi)
      rtree_stub.insert(std::make_pair(point(dca, phi, 0.0), 0));

      if(_write_ntp)
      {
        // store: ev, dca, phi, tanphi(=phi for now), x0(=dca)
        _ntp_stub->Fill(_nevent,
                        static_cast<float>(dca),
                        static_cast<float>(phi),
                        static_cast<float>(phi),
                        static_cast<float>(dca));

        _hHough->Fill(dca, phi);
      }

      // update Hough space ranges
      if(phi > phimax) phimax = phi;
      if(phi < phimin) phimin = phi;
      if(dca > dcamax) dcamax = dca;
      if(dca < dcamin) dcamin = dca;

      std::cout << " stub fit count: " << fcount
                << " phi: " << phi << " ( " << phimin << " ; " << phimax << " ) "
                << " dca: " << dca << " ( " << dcamin << " ; " << dcamax << " ) "
                << std::endl;
    }
  }

  // find clusters of (dca, phi) pairs in Hough space
  map<int,pair<float,float>> outtrkmap; // key: multiplicity, value: (dca, phi)
  int nout = 0;

  while(rtree_stub.size()>10)
  {
    std::vector<pointKey> allstubs;
    rtree_stub.query(
      bgi::intersects(
        box(point(dcamin, phimin, -1),
            point(dcamax, phimax,  1))),
      std::back_inserter(allstubs));

    map<int,pair<float,float>> trkmap;

    float phi_width = (phimax - phimin)/30;

    for(vector<pointKey>::iterator stub = allstubs.begin();
        stub!=allstubs.end();++stub)
    {
      float p_dca = stub->first.get<0>();
      float p_phi = stub->first.get<1>();

      float dca_width = dcaBinWidth(p_dca);

      // compute non-overlapping DCA interval for this region
      float dlow = 0.f;
      float dhigh = 0.f;
      dcaSearchRange(p_dca, dca_width, dcamin, dcamax, dlow, dhigh);

      std::cout << " Searching box around stub dca: " << p_dca
                << " phi: " << p_phi
                << " dca_w: " << dca_width
                << " phi_w: " << phi_width << std::endl;

      vector<pointKey> trkcand;
      rtree_stub.query(
        bgi::intersects(
          box(point(dlow,      p_phi - phi_width, -1),
              point(dhigh,     p_phi + phi_width,  1))),
        std::back_inserter(trkcand));

      int   ntrk   = trkcand.size();
      int   count  = 0;
      float dca_sum  = 0;
      float phi_sum  = 0;

      if(ntrk>=5)
      {
        for(vector<pointKey>::iterator ptrk = trkcand.begin();
            ptrk!=trkcand.end();++ptrk)
        {
          float trk_dca = ptrk->first.get<0>();
          float trk_phi = ptrk->first.get<1>();

          cout<< "    stub " << ntrk
              << " dca: " << trk_dca
              << " phi: " << trk_phi
              << endl;

          dca_sum += trk_dca;
          phi_sum += trk_phi;
          count++;
        }
        float mean_dca = (dca_sum/count);
        float mean_phi = (phi_sum/count);
        trkmap[ntrk]   = std::make_pair(mean_dca, mean_phi); 
              
        cout<< " stub in box " << ntrk
            << " dca_center: " << p_dca
            << " phi_center: " << p_phi
            << " dca_w: " << dca_width 
            << " phi_w: " << phi_width
            << " mean_dca " << mean_dca
            << " mean_phi " << mean_phi
            << endl;
      }
    }
    
    if( trkmap.size()>0)
    {
      int size_before = rtree_stub.size();
     
      cout << "mapend: " << trkmap.rbegin()->first
           << " votes | " << (trkmap.rbegin()->second).first
           << " dca | " << trkmap.rbegin()->second.second  << " phi | "<< endl;

      // remove stubs in the most populated Hough cell
      float best_dca = trkmap.rbegin()->second.first;
      float best_phi = trkmap.rbegin()->second.second;

      float dca_width = dcaBinWidth(best_dca);

      float dlow = 0.f;
      float dhigh = 0.f;
      dcaSearchRange(best_dca, dca_width, dcamin, dcamax, dlow, dhigh);

      std::cout<< " best dca: " << best_dca
               << " best phi: " << best_phi
               << " dca_w: " << dca_width
               << " phi_w: " << phi_width << std::endl;

      vector<pointKey> rmcand;
      rtree_stub.query(
        bgi::intersects(
          box(point(dlow,      best_phi - phi_width, -1),
              point(dhigh,     best_phi + phi_width,  1))),
        std::back_inserter(rmcand));

      for(vector<pointKey>::iterator rmstub = rmcand.begin();
          rmstub!=rmcand.end();++rmstub)
      {
        float rmp_dca = rmstub->first.get<0>();
        float rmp_phi = rmstub->first.get<1>();

        cout<< "    rm " <<  " dca: " << rmp_dca
            << " phi: " << rmp_phi 
            << endl;

        rtree_stub.remove(*rmstub);
      }

      int size_after = rtree_stub.size();
      if(size_before == size_after) {break;}

      // store this Hough peak (dca, phi)
      outtrkmap[nout++] = std::make_pair(best_dca, best_phi);
      cout << " tree size after remove: " << rtree_stub.size() << endl;
    }
    else
    {
      break;
    }
    trkmap.clear();
  }

  int numberofseeds = 0;
  bool keep_event = false;

  for (const auto& [key, value] : outtrkmap)
  {
  
    std::cout <<" out ev: " << _nevent << '[' << key << "] = "
              << value.first << " | " << value.second << std::endl;

    // Hough parameters
    float dca = value.first;
    float phi = value.second;

    // Line in standard form a x + b y + c = 0 using (phi, dca)
    // Choose normal n = (sin(phi), -cos(phi)) with |n|=1, c = -dca
    double a = std::sin(phi);
    double b = -std::cos(phi);
    double c = -dca;
      
    // geometry for intersections with circle of radius r
    double r = 100;
    double x0 = -a*c/(a*a+b*b);
    double y0 = -b*c/(a*a+b*b);
    double ax = NAN, ay = NAN, bx = NAN, by = NAN;

    if ((c*c) > (r*r*(a*a+b*b)+0.00000001))
    {
      if(Verbosity()>0) cout << "no points " << endl;
    }
    else if (std::abs(c*c - r*r*(a*a+b*b)) < 0.0000001) 
    {
      if(Verbosity()>0)
      { 
        puts ("1 point");
        cout << x0 << ' ' << y0 << endl;
      }
    }
    else 
    {
      double d = r*r - c*c/(a*a+b*b);
      double mult = std::sqrt(d / (a*a+b*b));
      ax = x0 + b * mult;
      bx = x0 - b * mult;
      ay = y0 - a * mult;
      by = y0 + a * mult;
      if(Verbosity()>0){
        puts ("2 points");
        cout << ax << ' ' << ay << '\n' << bx << ' ' << by << endl;
      }
    }

    float alpha = phi;  // direction angle of the track in XY

    vector<pointKey> trkclusters;
    vector<pointKey> lineclusters;

    // select all clusters in (r,phi,z) region
    rtree.query(
      bgi::intersects(
        box(point(RMIN,   PHIMIN, ZMIN),
            point(RMAX,   PHIMAX, ZMAX))),
      std::back_inserter(lineclusters));

    cout << "number line clus is " << lineclusters.size() << endl;
    
    // Check if track hits MVTX (dist to 0,0 < 2)
    float dist_origin = std::abs(c)/std::sqrt(a*a+b*b); // = |dca|

    if(_max_dist_to_origin>0)
    {
      if(dist_origin>_max_dist_to_origin)
      {
        cout << "dist continue: " << dist_origin
             << " > " << _max_dist_to_origin << endl;
        continue;
      }
    }

    for(vector<pointKey>::iterator clustertrk = lineclusters.begin();
        clustertrk!=lineclusters.end();++clustertrk)
    {
      // tree stores (r, phi, z)
      float rtrk   = clustertrk->first.get<0>();
      float phitrk = clustertrk->first.get<1>();

      float ptx = rtrk * std::cos(phitrk);
      float pty = rtrk * std::sin(phitrk);

      float res = std::abs(a*ptx + b*pty + c)/std::sqrt(a*a + b*b);
      
      cout << " x: " << ptx << " y: " << pty << " res " << res << endl;

      if(res<1.)
      {
        trkclusters.push_back(*clustertrk);
      }
    }

    cout << "number trk clus is " << trkclusters.size() << endl;
    if(trkclusters.size()>=_min_nclusters)
    {
      cout << "setting keep true: " << trkclusters.size()
           << " > " << _min_nclusters << endl;
      keep_event = true;
    }

    // Assemble tracks
    if(_create_tracks)
    {
      if(trkclusters.size()>=20)
      {
        const float nclus = static_cast<float>(trkclusters.size());
        _ntp_trk->Fill(static_cast<float>(_nevent),
                       static_cast<float>(dca),
                       static_cast<float>(phi),
                       alpha,
                       static_cast<float>(dist_origin),
                       nclus);

          _trk_ev    = _nevent;
          _trk_dca   = static_cast<float>(dca);
          _trk_phi   = static_cast<float>(phi);
          _trk_alpha = alpha;
          _trk_dist0 = static_cast<float>(dist_origin);
          _trk_nclus = static_cast<int>(trkclusters.size());

          _trk_r.clear();
          _trk_phi_clu.clear();
          _trk_z.clear();
          _trk_cluskey.clear();

          for(std::vector<pointKey>::const_iterator it = trkclusters.begin();
              it != trkclusters.end(); ++it)
          {
            float r_clu   = it->first.get<0>();
            float phi_clu = it->first.get<1>();
            float z_clu   = it->first.get<2>();

            _trk_r.push_back(r_clu);
            _trk_phi_clu.push_back(phi_clu);
            _trk_z.push_back(z_clu);
            _trk_cluskey.push_back(static_cast<ULong64_t>(it->second));
          }

          _track_tree->Fill();

        auto trackseed = std::make_unique<TrackSeed_v2>();
        for(const auto& cluskeys:trkclusters)
        {
          trackseed->insert_cluster_key(cluskeys.second);
        }

        m_seedContainer->insert(trackseed.get()); 
        cout << "number trk keys is " << trackseed->size_cluster_keys() << endl;
        numberofseeds++;
      }
    }
  }
  
  cout << "number of seeds is " << numberofseeds << endl;
  if(_write_ntp)
  {
    _ntp_max->Fill(_nevent,dcamin,dcamax,phimin,phimax);
  }
  if(!keep_event)
  {
    cout << " ABORT !keep " << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int HoughTrackFinder::Setup(PHCompositeNode* topNode)
{
  cout << "Called Setup" << endl;
  cout << "topNode:" << topNode << endl;
  GetNodes(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

int HoughTrackFinder::End(PHCompositeNode* topNode)
{
  if(topNode==nullptr){
    std::cout << PHWHERE << "No topNode, bailing."
              << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
    
  if(_write_ntp){
    _tfile->cd();
    _ntp_cos->Write();
    _ntp_stub->Write();
    _ntp_max->Write();
    _ntp_trk->Write();
    _hHough->Write();
    _track_tree->Write();
    _tfile->Close();
    cout << "Called End " << endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

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
    _ntp_stub = new TNtuple("ntp_stub", "cos stub info","ev:dca:phi:theta:z0");
    // store (dca min/max, phi min/max)
    _ntp_max = new TNtuple("ntp_max", "cos stub info","ev:dcamin:dcamax:phimin:phimax");
    // full track params: use dca, phi, alpha=phi
    _ntp_trk = new TNtuple("ntp_trk", "full track params", "ev:dca:phi:alpha:dist0:nclus");

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
                                int &count,
                                double &dca,
                                double &phi,
                                double &z0,
                                double &theta)
{
  // search window in (r, phi, z)
  const float m1_dr   = 1.5f;   // cm in radius
  const float m1_dphi = 0.05f;  // rad
  const float m1_dz   = 2.0f;   // cm

  std::vector<pointKey> boxclusters;

  // 3D search box in cylindrical coordinates
  search_rtree.query(
      bgi::intersects(
          box(
              point(pointr   - m1_dr,
                    pointphi - m1_dphi,
                    pointz   - m1_dz),
              point(pointr   + m1_dr,
                    pointphi + m1_dphi,
                    pointz   + m1_dz))),
      std::back_inserter(boxclusters));

  // Debug: center in xy for sanity
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

  // collect 3D points in Cartesian (x,y,z)
  std::vector<Eigen::Vector3d> pts3;
  pts3.reserve(boxclusters.size());

  for (auto pbox = boxclusters.begin(); pbox != boxclusters.end(); ++pbox)
  {
    float r_hit   = pbox->first.get<0>();
    float phi_hit = pbox->first.get<1>();
    float z_hit   = pbox->first.get<2>();

    float x = r_hit * std::cos(phi_hit);
    float y = r_hit * std::sin(phi_hit);

    pts3.emplace_back(x, y, z_hit);
  }

  count = static_cast<int>(pts3.size());

  const int MIN_STUB_HITS = 2;
  if (count < MIN_STUB_HITS)
  {
    dca   = 0.;
    phi   = 0.;
    z0    = 0.;
    theta = 0.;
    return;
  }

  // 1) centroid
  Eigen::Vector3d mean = Eigen::Vector3d::Zero();
  for (const auto &p : pts3) mean += p;
  mean /= static_cast<double>(count);

  // 2) covariance matrix
  Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
  for (const auto &p : pts3)
  {
    Eigen::Vector3d q = p - mean;
    cov += q * q.transpose();
  }

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(cov);
  if (es.info() != Eigen::Success)
  {
    std::cout << " get_stub: Eigen decomposition failed" << std::endl;
    dca   = 0.;
    phi   = 0.;
    z0    = 0.;
    theta = 0.;
    count = 0;
    return;
  }

  // eigenvalues are sorted ascending; largest is index 2
  Eigen::Vector3d dir = es.eigenvectors().col(2);
  dir.normalize();

  // enforce a consistent sign convention (e.g. uz > 0)
  if (dir.z() < 0) dir = -dir;

  double ux = dir.x();
  double uy = dir.y();
  double uz = dir.z();

  double norm_xy = std::sqrt(ux*ux + uy*uy);
  if (norm_xy < 1e-6)
  {
    std::cout << " get_stub: degenerate direction in XY" << std::endl;
    dca   = 0.;
    phi   = 0.;
    z0    = 0.;
    theta = 0.;
    count = 0;
    return;
  }

  // 3) track direction parameters in XY and w.r.t. beam axis
  phi   = std::atan2(uy, ux);                 // azimuth
  theta = std::atan2(norm_xy, uz);            // polar angle to beam (z)

  // 4) DCA and z0 from closest approach in XY
  //    line: r(t) = mean + t * dir
  double x0 = mean.x();
  double y0 = mean.y();
  double zc = mean.z();

  double denom_xy = ux*ux + uy*uy;
  double t_dca    = -(x0*ux + y0*uy) / denom_xy;

  double x_dca = x0 + ux * t_dca;
  double y_dca = y0 + uy * t_dca;
  double z_dca = zc + uz * t_dca;

  double r_dca = std::sqrt(x_dca*x_dca + y_dca*y_dca);

  // sign convention for DCA (left/right of direction)
  double sign = ((x_dca*uy - y_dca*ux) >= 0.0) ? 1.0 : -1.0;
  dca = sign * r_dca;
  z0  = z_dca;

  // 5) chi2 / ndf using 3D distance to the line
  double chi2 = 0.0;
  for (const auto &p : pts3)
  {
    Eigen::Vector3d diff = p - mean;
    double t = diff.dot(dir);
    Eigen::Vector3d proj = mean + t * dir;
    Eigen::Vector3d res  = p - proj;
    chi2 += res.squaredNorm();
  }

  const double ndof        = count - 2;
  const double chi2ndf     = (ndof > 0) ? chi2 / ndof : 1e9;
  const double MAX_CHI2NDF = 0.5;  // tweak as needed

  std::cout << " stub 3D fit: nhit=" << count
            << " dca="    << dca
            << " phi="    << phi
            << " theta="  << theta
            << " z0="     << z0
            << " chi2/ndf=" << chi2ndf
            << std::endl;

  if (!std::isfinite(dca)   || !std::isfinite(phi)   ||
      !std::isfinite(theta) || !std::isfinite(z0)    ||
      chi2ndf > MAX_CHI2NDF)
  {
    std::cout << " get_stub: bad fit, rejecting stub" << std::endl;
    dca   = 0.;
    phi   = 0.;
    z0    = 0.;
    theta = 0.;
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

  float thetamin =  99999999999.9;
  float thetamax = -99999999999.9;

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

    int    fcount   = 0;
    double fdca     = 0.;
    double fphi     = 0.;
    double fz0      = 0.;
    double ftheta   = 0.;

    // ...

    std::cout << " Processing cluster at r=" << pointr
              << " phi=" << pointphi
              << " z="   << pointz
              << " (x="  << pointx << ", y=" << pointy << ")"
              << std::endl;

    // 3D stub fit: returns (dca, phi, z0, theta)
    get_stub(rtree, pointr, pointphi, pointz,
            fcount, fdca, fphi, fz0, ftheta);

if (fcount <= 0) continue;
if (!std::isfinite(fdca) || !std::isfinite(fphi)) continue;
if (fdca == 0.0 && fphi == 0.0) continue;

    if (std::isfinite(fdca))
    {
      // convert (slope, intercept) -> (phi, dca) in XY
 //     const double phi = std::atan(fslope);
   //   const double dca = fintercept / std::sqrt(fslope*fslope + 1.0);
   //   if (dca ==0 && phi ==0) continue;

      // insert stub in Hough space (dca, phi)
      rtree_stub.insert(std::make_pair(point(fdca, fphi, ftheta), 0));

      if(_write_ntp)
      {
        // store: ev, dca, phi, tanphi(=phi for now), x0(=dca)
        _ntp_stub->Fill(_nevent,
                        static_cast<float>(fdca),
                        static_cast<float>(fphi),
                        static_cast<float>(ftheta),
                        static_cast<float>(fz0));

        _hHough->Fill(fdca, fphi);
      }

      // update Hough space ranges
      if (fphi > phimax) phimax = fphi;
      if (fphi < phimin) phimin = fphi;
      if (fdca > dcamax) dcamax = fdca;
      if (fdca < dcamin) dcamin = fdca;
      if (ftheta > thetamax) thetamax = ftheta;
      if (ftheta < thetamin) thetamin = ftheta;

      std::cout << " stub fit count: " << fcount
                << " phi: " << fphi << " ( " << phimin << " ; " << phimax << " ) "
                << " dca: " << fdca << " ( " << dcamin << " ; " << dcamax << " ) "
                << " theta: " << ftheta
                << " z0: " << fz0
                << std::endl;
    }
  }

  // find clusters of (dca, phi) pairs in Hough space
  std::map<int, HoughPeak> outtrkmap;  // key: multiplicity, value: (dca, phi)
  int nout = 0;

  while(rtree_stub.size()>10)
  {
    std::vector<pointKey> allstubs;
    rtree_stub.query(
      bgi::intersects(
        box(point(dcamin, phimin, thetamin),
            point(dcamax, phimax, thetamax))),
      std::back_inserter(allstubs));

    std::map<int, HoughPeak> trkmap;

    float phi_width = (phimax - phimin)/30;
    float theta_width = (thetamax - thetamin)/20;

    for(vector<pointKey>::iterator stub = allstubs.begin();
        stub!=allstubs.end();++stub)
    {
      float p_dca = stub->first.get<0>();
      float p_phi = stub->first.get<1>();
      float p_theta = stub->first.get<2>();

      float dca_width = dcaBinWidth(p_dca);

      // compute non-overlapping DCA interval for this region
      float dlow = 0.f;
      float dhigh = 0.f;
      dcaSearchRange(p_dca, dca_width, dcamin, dcamax, dlow, dhigh);

      std::cout << " Searching box around stub dca: " << p_dca
                << " phi: " << p_phi
                << " theta: " << p_theta
                << " dca_w: " << dca_width
                << " phi_w: " << phi_width
                << " theta_w: "<< theta_width << std::endl;

      vector<pointKey> trkcand;
      rtree_stub.query(
        bgi::intersects(
          box(point(dlow,      p_phi - phi_width, p_theta - theta_width),
              point(dhigh,     p_phi + phi_width, p_theta + theta_width))),
        std::back_inserter(trkcand));

      int   ntrk   = trkcand.size();
      int   count  = 0;
      float dca_sum  = 0;
      float phi_sum  = 0;
      float theta_sum = 0.f; 
      if(ntrk>=3)
      {
        for(vector<pointKey>::iterator ptrk = trkcand.begin();
            ptrk!=trkcand.end();++ptrk)
        {
          float trk_dca = ptrk->first.get<0>();
          float trk_phi = ptrk->first.get<1>();
          float trk_theta = ptrk->first.get<2>();

          cout << "    stub " << ntrk
              << " dca: "   << trk_dca
              << " phi: "   << trk_phi
              << " theta: " << trk_theta
              << endl;


          dca_sum += trk_dca;
          phi_sum += trk_phi;
          theta_sum += trk_theta;
          count++;
        }
        float mean_dca = (dca_sum/count);
        float mean_phi = (phi_sum/count);
        float mean_theta = theta_sum / count; 

        HoughPeak peak;
        peak.dca   = mean_dca;
        peak.phi   = mean_phi;
        peak.theta = mean_theta;
        trkmap[ntrk] = peak;
              
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
          << " votes | " << trkmap.rbegin()->second.dca
          << " dca | "    << trkmap.rbegin()->second.phi
          << " phi | "    << trkmap.rbegin()->second.theta
          << " theta | "  << endl;

      float best_dca   = trkmap.rbegin()->second.dca;
      float best_phi   = trkmap.rbegin()->second.phi;
      float best_theta = trkmap.rbegin()->second.theta;  

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
          box(point(dlow,      best_phi - phi_width,  best_theta - theta_width),
              point(dhigh,     best_phi + phi_width,  best_theta + theta_width))),
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
      outtrkmap[nout++] = trkmap.rbegin()->second;
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

  for (const auto& kv : outtrkmap)
  {
    int       key   = kv.first;
    const HoughPeak& value = kv.second;
  

    std::cout << " out ev: " << _nevent << '[' << key << "] = "
            << value.dca   << " | " 
            << value.phi   << " | "
            << value.theta << std::endl;

    // Hough parameters
    float dca = value.dca;
    float phi = value.phi;
    // float theta = value.theta;

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

      if(res<0.5)
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

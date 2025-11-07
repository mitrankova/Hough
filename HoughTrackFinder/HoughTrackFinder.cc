//============================================================
// HoughTrackFinder.cc — 3D straight-line Hough (r,φ,z)
// Optimized: k-NN local stubs + grid clustering in (k,z0)
// Extended: Extract (rho, theta) in XY plane
//============================================================

#include "HoughTrackFinder.h"

#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/ActsGeometry.h>
#include <g4detectors/PHG4TpcGeomContainer.h>
#include <g4detectors/PHG4TpcGeom.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <TFile.h>
#include <TNtuple.h>

#include <cmath>
#include <algorithm>
#include <map>
#include <set>
#include <iostream>

#ifdef _OPENMP
  #include <omp.h>
#endif

// Boost geometry/index (R-tree)
/*#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

using std::sqrt;
namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;
*/


// ----------------- Private impl for the R-tree -----------------
/*struct RtreeImpl
{
  typedef bg::model::point<float, 3, bg::cs::cartesian> Point3D;
  typedef bg::model::box<Point3D> Box3D;
  typedef std::pair<Point3D, unsigned int> Leaf; // (x,y,z) with index
  typedef bgi::rtree< Leaf, bgi::quadratic<16> > RTree;
  RTree tree;
};*/

//============================================================
// ctor/dtor
//============================================================
HoughTrackFinder::HoughTrackFinder(const std::string &name,
                                   const std::string &outfile)
  : SubsysReco(name)
  , _min_layer(7)
  , _max_layer(55)
  , _knn(15)
  , _stub_window_r(3)
  , _stub_window_phi(0.0053073*8)
  , _stub_window_z(3)
  , _max_z_residual(3.0)
  , _min_stubs(6)
  , _min_hits(20)
  , _slope_bin(0.02f)
  , _z0_bin(2.0f)
  , _max_tracks(16)
  , m_outputFileName(outfile)
  , m_outputFile(0)
  , _ntp_hits(0)
  , _ntp_stubs(0)
  , _nevent(0)
  , _rtree(0)
{
}

HoughTrackFinder::~HoughTrackFinder()
{
  if (_rtree) { delete _rtree; _rtree = 0; }
}

//============================================================
// Init
//============================================================
int HoughTrackFinder::Init(PHCompositeNode * /*topNode*/)
{
  m_outputFile = new TFile(m_outputFileName.c_str(), "RECREATE");
  if (!m_outputFile || m_outputFile->IsZombie())
  {
    std::cerr << "HoughTrackFinder::Init - cannot create " << m_outputFileName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  _ntp_hits   = new TNtuple("hits",   "hits",   "event:x:y:z");
  _ntp_stubs  = new TNtuple("stubs",  "stubs",  "event:d0:z0:phi0:theta:qOverP:npts");


  return Fun4AllReturnCodes::EVENT_OK;
}

//============================================================
// process_event
//============================================================
int HoughTrackFinder::process_event(PHCompositeNode *topNode)
{
  m_hits.clear();
  m_stubs.clear();

  if (_rtree) { delete _rtree; _rtree = 0; }

  int ret = collectHits(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
  if (m_hits.size() < 20) { _nevent++; return Fun4AllReturnCodes::EVENT_OK; }

  buildRtree();
  fitLocalStubs();

  if (m_stubs.empty()) { _nevent++; return Fun4AllReturnCodes::EVENT_OK; }

 // clusterStubs();
//  if (m_tracks.empty()) { _nevent++; return Fun4AllReturnCodes::EVENT_OK; }



if (_ntp_stubs)
{
  for (size_t i = 0; i < m_stubs.size(); ++i)
  {
    const Stub &s = m_stubs[i];
    _ntp_stubs->Fill((float)_nevent,
                     s.d0,
                     s.z0,
                     s.phi0,
                     s.theta,
                     s.qOverP,
                     (float)s.npoints);
  }
}


  
  _nevent++;
  return Fun4AllReturnCodes::EVENT_OK;
}

//============================================================
// End
//============================================================
int HoughTrackFinder::End(PHCompositeNode * /*topNode*/)
{
  if (m_outputFile)
  {
    m_outputFile->cd();
    if (_ntp_hits)   _ntp_hits->Write();
    if (_ntp_stubs)  _ntp_stubs->Write();
    m_outputFile->Close();
    delete m_outputFile; m_outputFile = 0;
  }
  if (_rtree) { delete _rtree; _rtree = 0; }
  return Fun4AllReturnCodes::EVENT_OK;
}

//============================================================
// Helpers
//============================================================

int HoughTrackFinder::collectHits(PHCompositeNode *topNode)
{
  TrkrHitSetContainer *hitmap = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  PHG4TpcGeomContainer *tpcGeom = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  ActsGeometry *geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!hitmap || !tpcGeom || !geometry) return Fun4AllReturnCodes::ABORTEVENT;

  TrkrHitSetContainer::ConstRange allsets = hitmap->getHitSets();
  for (TrkrHitSetContainer::ConstIterator itset = allsets.first; itset != allsets.second; ++itset)
  {
    TrkrDefs::hitsetkey hitsetkey = itset->first;
    TrkrHitSet *hitset = itset->second;
    if (TrkrDefs::getTrkrId(hitsetkey) != TrkrDefs::tpcId) continue;

    unsigned int layer = TrkrDefs::getLayer(hitsetkey);
    if (layer < _min_layer || layer > _max_layer) continue;

    unsigned int side  = TpcDefs::getSide(hitsetkey);
    PHG4TpcGeom *geoLayer = tpcGeom->GetLayerCellGeom(layer);
    if (!geoLayer) continue;

    double radius = geoLayer->get_radius();
    float AdcClockPeriod = geoLayer->get_zstep();

    TrkrHitSet::ConstRange range = hitset->getHits();
    for (TrkrHitSet::ConstIterator hitr = range.first; hitr != range.second; ++hitr)
    {
      TrkrDefs::hitkey hitkey = hitr->first;
      auto *hit = hitr->second;
      int m_adc = hit->getAdc();
      if (m_adc < 40) continue; // adc cut
      unsigned int pad = TpcDefs::getPad(hitkey);
      double phi = geoLayer->get_phicenter(pad, side);
      Acts::Vector3 glob = geometry->getGlobalPositionTpc(hitsetkey, hitkey, phi, radius, AdcClockPeriod);

      HitXYZ h;
      h.x = (float)glob.x();
      h.y = (float)glob.y();
      h.z = (float)glob.z();
      m_hits.push_back(h);

      if (_ntp_hits) _ntp_hits->Fill((float)_nevent, h.x, h.y, h.z);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void HoughTrackFinder::buildRtree()
{
  if (_rtree) { delete _rtree; _rtree = 0; }
  _rtree = new RtreeImpl;
  for (unsigned int i = 0; i < m_hits.size(); ++i)
  {
    const HitXYZ &h = m_hits[i];
    RtreeImpl::Point3D p(h.x, h.y, h.z);
    _rtree->tree.insert(std::make_pair(p, i));
  }
}



void HoughTrackFinder::fitLocalStubs()
{
  if (!_rtree) return;
  m_stubs.clear();
  m_stubs.reserve(m_hits.size());

  const unsigned int N = (unsigned int)m_hits.size();
  std::vector< std::vector<Stub> > thread_stubs;
  unsigned int nthreads=1;
//#ifdef _OPENMP
//  nthreads = (unsigned int)omp_get_max_threads();
//#endif
  thread_stubs.resize(nthreads);

//#ifdef _OPENMP
//#pragma omp parallel for schedule(dynamic)
//#endif
  for (int i=0;i<(int)N;++i)
  {
    const HitXYZ &h = m_hits[(unsigned int)i];
    std::vector<RtreeImpl::Leaf> neigh;
    neigh.reserve(_knn ? _knn : 32);

    if (_knn > 0)
      _rtree->tree.query(bgi::nearest(RtreeImpl::Point3D(h.x,h.y,h.z),(int)_knn),std::back_inserter(neigh));
    else
    {
      float rwin = _stub_window_r;
      float phiwin = _stub_window_phi;
      float zwin = _stub_window_z;

      float r   = sqrt(h.x*h.x + h.y*h.y);
      float phi = atan2(h.y, h.x);
      float dx = fabs(cos(phi) * rwin) + fabs(r * sin(phi) * phiwin);
      float dy = fabs(sin(phi) * rwin) + fabs(r * cos(phi) * phiwin);
      float dz = zwin;
          std::cout<<"Box for hit "<<i<<" : "<<h.x - dx<<" "<<h.y - dy<<" "<<h.z - dz<<" to "
                    <<h.x + dx<<" "<<h.y + dy<<" "<<h.z + dz<<"    point "<<h.x<<" "<<h.y<<" "<<h.z<<"    delta "<<dx<<" "<<dy<<" "<<dz<<std::endl;
      RtreeImpl::Box3D box(
          RtreeImpl::Point3D(h.x - dx, h.y - dy, h.z - dz),
          RtreeImpl::Point3D(h.x + dx, h.y + dy, h.z + dz));

      _rtree->tree.query(bgi::intersects(box),std::back_inserter(neigh));
    }
    std::cout<<"Found "<<neigh.size()<<" neighbors for hit "<<i<<std::endl;
    if (neigh.size()<3) continue;


float d0=0, z0=0, phi0=0, theta=0, qOverP=0;
if (!helix_fit_3D(neigh, m_hits, d0, z0, phi0, theta, qOverP)) continue;


//#ifdef _OPENMP
//    unsigned int tid=(unsigned int)omp_get_thread_num();
//#else
    unsigned int tid=0;
//#endif
    Stub s; 
    s.d0 = d0;
    s.z0 = z0;
    s.phi0 = phi0;
    s.theta = theta;
    s.qOverP = qOverP;
    s.npoints = (int)neigh.size();

    thread_stubs[tid].push_back(s);
  }

  for (unsigned int t=0;t<thread_stubs.size();++t)
    m_stubs.insert(m_stubs.end(),thread_stubs[t].begin(),thread_stubs[t].end());

 /* if (_ntp_stubs)
  {
    for (size_t i=0;i<m_stubs.size();++i)
        _ntp_stubs->Fill((float)_nevent, d0, z0, phi0, theta, qOverP, (float)npoints);
  }*/
}


// Build a best-fit straight line r(t) = r0 + t * vhat from the neighbor hits.
// Use a simple PCA/power-iteration or (good enough) fit vhat as the
// normalized difference of the farthest pair of points in neigh.

void HoughTrackFinder::fit_straight_line_xyz(
  const std::vector<RtreeImpl::Leaf>& neigh,
  const std::vector<HitXYZ>& hits,
  double& r0x, double& r0y, double& r0z,
  double& vx,  double& vy,  double& vz)
{
  // centroid
  double sx=0, sy=0, sz=0; size_t n=neigh.size();
  for (size_t i=0;i<n;++i){ const HitXYZ& h=hits[neigh[i].second]; sx+=h.x; sy+=h.y; sz+=h.z; 
  std::cout<<"  neighbor "<<i<<" : "<<h.x<<" "<<h.y<<" "<<h.z<<std::endl;
  }
  r0x=sx/n; r0y=sy/n; r0z=sz/n;

  // crude direction: farthest-pair heuristic (robust enough for stubs)
  size_t i0=0, i1=0; double maxd2=-1;
  for(size_t i=0;i<n;++i){
    const HitXYZ& a=hits[neigh[i].second];
    for(size_t j=i+1;j<n;++j){
      const HitXYZ& b=hits[neigh[j].second];
      double dx=b.x-a.x, dy=b.y-a.y, dz=b.z-a.z;
      double d2=dx*dx+dy*dy+dz*dz;
      if(d2>maxd2){maxd2=d2;i0=i;i1=j;}
    }
  }
  const HitXYZ& A=hits[neigh[i0].second];
  const HitXYZ& B=hits[neigh[i1].second];
  vx=B.x-A.x; vy=B.y-A.y; vz=B.z-A.z;
  double norm=std::sqrt(vx*vx+vy*vy+vz*vz); if(norm<1e-12){ vx=1; vy=0; vz=0; } else { vx/=norm; vy/=norm; vz/=norm; }
  
}

bool HoughTrackFinder::to_straight_perigee(
  const std::vector<RtreeImpl::Leaf>& neigh,
  const std::vector<HitXYZ>& hits,
  float& d0, float& z0, float& phi0, float& theta, float& qOverP)
{
  double r0x, r0y, r0z, vx, vy, vz;
  std::cout<<"HoughTrackFinder::to_straight_perigee "<<neigh.size()<<" points"<<std::endl;
  fit_straight_line_xyz(neigh, hits, r0x, r0y, r0z, vx, vy, vz);
std::cout<<"  straight line: r0 = ("<<r0x<<", "<<r0y<<", "<<r0z<<") v = ("<<vx<<", "<<vy<<", "<<vz<<")"<<std::endl;
  // direction → angles
  phi0  = (float)std::atan2(vy, vx);
  double vT2 = vx*vx + vy*vy;

  // closest approach to z-axis
  double tstar = 0.0;
  if (vT2 > 1e-12) tstar = -(r0x*vx + r0y*vy)/vT2;
  double xP = r0x + tstar*vx;
  double yP = r0y + tstar*vy;
  double zP = r0z + tstar*vz;

  // perigee params
  d0   = (float)(-xP*std::sin(phi0) + yP*std::cos(phi0));
  z0   = (float)zP;
  double vz_clamped = std::max(-1.0, std::min(1.0, vz));
  theta = (float)std::acos(vz_clamped);   // polar angle to z-axis
  qOverP = 0.0f;                          // straight line ⇒ zero curvature
std::cout<<"  perigee params: d = "<<d0<<" z0 = "<<z0<<" phi0 = "<<phi0<<" theta = "<<theta<<std::endl;
  return std::isfinite(d0) && std::isfinite(z0) && std::isfinite(phi0)
      && std::isfinite(theta);
}

//============================================================
// 3D Helix fit (perigee parameters) with parameter constraints
//============================================================
bool HoughTrackFinder::helix_fit_3D(
  const std::vector<RtreeImpl::Leaf>& neigh,
  const std::vector<HitXYZ>& hits,
  float& d0, float& z0, float& phi0, float& theta, float& qOverP)
{
  std::cout<<"!!!!!!!!HoughTrackFinder::helix_fit_3D "<<neigh.size()<<" points"<<std::endl;
  const size_t n = neigh.size();
  if (n < 4) return false;

  // ------------------------------
  // Constants and constraints
  // ------------------------------
  const double kBField_T = 1.4;           // Tesla
  //const double kMinRadius_cm = 0.5;       // min valid circle radius
  const double kMaxRadius_cm = 140.0;     // max valid radius (~5 m)
  const double kMaxD0_cm = 5.0;          // max transverse impact
  const double kMaxZ0_cm = 20.0;         // max z at perigee
  const double kMinTheta_rad = 0;      // 
  const double kMaxTheta_rad = 10;      // just below π

  // ------------------------------
  // 1. Extract hit coordinates
  // ------------------------------
  std::vector<float> x(n), y(n), z(n);
  for (size_t i = 0; i < n; ++i)
  {
    const HitXYZ& h = hits[neigh[i].second];
    x[i] = h.x;
    y[i] = h.y;
    z[i] = h.z;
    std::cout<<"Hit "<<i<<" : x = "<<x[i]<<"   y = "<<y[i]<<"   z = "<<z[i]<<std::endl;
  }

  // ------------------------------
  // 2. Circle fit (Kåsa method)
  // ------------------------------
  double sum_x = 0, sum_y = 0, sum_x2 = 0, sum_y2 = 0;
  double sum_xy = 0, sum_x3 = 0, sum_y3 = 0, sum_x1y2 = 0, sum_x2y1 = 0;
  for (size_t i = 0; i < n; ++i)
  {
    const double xi = x[i];
    const double yi = y[i];
    const double xi2 = xi*xi;
    const double yi2 = yi*yi;

    sum_x += xi;
    sum_y += yi;
    sum_x2 += xi2;
    sum_y2 += yi2;
    sum_xy += xi*yi;
    sum_x3 += xi2*xi;
    sum_y3 += yi2*yi;
    sum_x1y2 += xi*yi2;
    sum_x2y1 += xi2*yi;
  }

  const double C = n*sum_x2 - sum_x*sum_x;
  const double D = n*sum_xy - sum_x*sum_y;
  const double E = n*sum_y2 - sum_y*sum_y;
  const double G = 0.5*(n*(sum_x3 + sum_x1y2) - (sum_x2 + sum_y2)*sum_x);
  const double H = 0.5*(n*(sum_y3 + sum_x2y1) - (sum_x2 + sum_y2)*sum_y);
  const double denom = (C*E - D*D);
  if (fabs(denom) < 1e-12) return false;

  const double xc = (G*E - D*H)/denom;
  const double yc = (C*H - D*G)/denom;
  const double R = std::sqrt((sum_x2 + sum_y2 - 2*xc*sum_x - 2*yc*sum_y)/n + xc*xc + yc*yc);
std::cout<<"Fitted circle:   xc = "<<xc<<"   yc = "<<yc<<"   R = "<<R<<std::endl;
if (!std::isfinite(R) || R > kMaxRadius_cm) {
  return to_straight_perigee(neigh, hits, d0, z0, phi0, theta, qOverP);
}
  // ------------------------------
  // 3. Compute phi0 and d0
  // ------------------------------
  const double rc = std::sqrt(xc*xc + yc*yc);
  phi0 = std::atan2(yc, xc) + M_PI/2.0;
  d0 = rc - R;
  if (fabs(d0) > kMaxD0_cm) return false;

  // ------------------------------
  // 4. Charge sign by rotation sense
  // ------------------------------
  double cross = 0;
  for (size_t i = 1; i < n; ++i)
    cross += (x[i-1]*y[i] - y[i-1]*x[i]);
  const int qsign = (cross > 0 ? +1 : -1);

  // ------------------------------
  // 5. Compute arc lengths s_i along helix
  // ------------------------------
  std::vector<double> s(n, 0.0);
  for (size_t i = 1; i < n; ++i)
  {
    double phi_prev = std::atan2(y[i-1]-yc, x[i-1]-xc);
    double phi_curr = std::atan2(y[i]-yc, x[i]-xc);
    double dphi = phi_curr - phi_prev;
    while (dphi > M_PI)  dphi -= 2*M_PI;
    while (dphi < -M_PI) dphi += 2*M_PI;
    s[i] = s[i-1] + qsign * R * dphi;
  }

  // ------------------------------
  // 6. Linear fit z(s)
  // ------------------------------
  double sum_s = 0, sum_s2 = 0, sum_z = 0, sum_sz = 0;
  for (size_t i = 0; i < n; ++i)
  {
    sum_s  += s[i];
    sum_s2 += s[i]*s[i];
    sum_z  += z[i];
    sum_sz += s[i]*z[i];
  }
  const double denom2 = (n*sum_s2 - sum_s*sum_s);
  if (fabs(denom2) < 1e-12) return false;

  const double slope = (n*sum_sz - sum_s*sum_z)/denom2; // tan(lambda)
  z0 = (sum_z - slope*sum_s)/n;
  if (fabs(z0) > kMaxZ0_cm) return false;

  const double tan_lambda = slope;
  theta = std::atan(1.0 / tan_lambda);
  if (theta < kMinTheta_rad || theta > kMaxTheta_rad || !std::isfinite(theta)) return false;

  // ------------------------------
  // 7. Compute q/p with cm→m conversion
  // ------------------------------
  const double R_m = R * 0.01; // cm → m
  const double pt = 0.3 * kBField_T * R_m;
  qOverP = (qsign * 1.0) / std::sqrt(pt*pt + (pt*tan_lambda)*(pt*tan_lambda));

  // ------------------------------
  // 8. Validity check
  // ------------------------------
  if (!std::isfinite(qOverP) || fabs(qOverP) > 10.0) return false;

  return true;
}




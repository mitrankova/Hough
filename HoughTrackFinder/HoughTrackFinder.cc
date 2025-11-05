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
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

using std::sqrt;
namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

// ----------------- Private impl for the R-tree -----------------
struct HoughTrackFinder::RtreeImpl
{
  typedef bg::model::point<float, 3, bg::cs::cartesian> Point3D;
  typedef bg::model::box<Point3D> Box3D;
  typedef std::pair<Point3D, unsigned int> Leaf; // (x,y,z) with index
  typedef bgi::rtree< Leaf, bgi::quadratic<16> > RTree;
  RTree tree;
};

//============================================================
// ctor/dtor
//============================================================
HoughTrackFinder::HoughTrackFinder(const std::string &name,
                                   const std::string &outfile)
  : SubsysReco(name)
  , _min_layer(7)
  , _max_layer(55)
  , _knn(15)
  , _stub_window_xy(5.0)
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
  , _ntp_tracks(0)
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
  _ntp_stubs  = new TNtuple("stubs",  "stubs",  "event:slope:z0:rho:theta:npts");
  _ntp_tracks = new TNtuple("tracks", "tracks", "event:slope:z0:rho:theta:nhits");

  return Fun4AllReturnCodes::EVENT_OK;
}

//============================================================
// process_event
//============================================================
int HoughTrackFinder::process_event(PHCompositeNode *topNode)
{
  m_hits.clear();
  m_stubs.clear();
  m_tracks.clear();
  if (_rtree) { delete _rtree; _rtree = 0; }

  int ret = collectHits(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;
  if (m_hits.size() < 20) { _nevent++; return Fun4AllReturnCodes::EVENT_OK; }

  buildRtree();
  fitLocalStubs();
  computeStubXYParameters();
  if (m_stubs.empty()) { _nevent++; return Fun4AllReturnCodes::EVENT_OK; }

  clusterStubs();
  if (m_tracks.empty()) { _nevent++; return Fun4AllReturnCodes::EVENT_OK; }

  assignHitsToTracks();
  computeXYParameters();  // <-- new step

if (_ntp_stubs)
{
  for (size_t i = 0; i < m_stubs.size(); ++i)
  {
    _ntp_stubs->Fill((float)_nevent,
                     m_stubs[i].slope,
                     m_stubs[i].z0,
                     m_stubs[i].rho,
                     m_stubs[i].theta,
                     (float)m_stubs[i].npoints);
  }
}

  // Write output
  if (_ntp_tracks)
  {
    for (size_t i = 0; i < m_tracks.size(); ++i)
    {
      if (m_tracks[i].nhits >= (int)_min_hits)
      {
        _ntp_tracks->Fill((float)_nevent,
                          m_tracks[i].slope, m_tracks[i].z0,
                          (float)m_tracks[i].nhits,
                          m_tracks[i].rho, m_tracks[i].theta);
      }
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
    if (_ntp_tracks) _ntp_tracks->Write();
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
    HoughTrackFinder::RtreeImpl::Point3D p(h.x, h.y, h.z);
    _rtree->tree.insert(std::make_pair(p, i));
  }
}

// linear regression helper
static inline bool linear_fit_rz(const std::vector<std::pair<float,float> > &rz,
                                 float &slope, float &z0)
{
  if (rz.size() < 3) return false;
  double sumr=0,sumz=0,sumr2=0,sumrz=0;
  const unsigned int n=(unsigned int)rz.size();
  for (unsigned int i=0;i<n;++i)
  {
    sumr+=rz[i].first;
    sumz+=rz[i].second;
    sumr2+=rz[i].first*rz[i].first;
    sumrz+=rz[i].first*rz[i].second;
  }
  const double denom = (double)n*sumr2 - sumr*sumr;
  if (fabs(denom)<1e-8) return false;
  slope=(float)(((double)n*sumrz - sumr*sumz)/denom);
  z0=(float)((sumz - slope*sumr)/n);
  return true;
}

void HoughTrackFinder::fitLocalStubs()
{
  if (!_rtree) return;
  m_stubs.clear();
  m_stubs.reserve(m_hits.size());

  const unsigned int N = (unsigned int)m_hits.size();
  std::vector< std::vector<Stub> > thread_stubs;
  unsigned int nthreads=1;
#ifdef _OPENMP
  nthreads = (unsigned int)omp_get_max_threads();
#endif
  thread_stubs.resize(nthreads);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int i=0;i<(int)N;++i)
  {
    const HitXYZ &h = m_hits[(unsigned int)i];
    std::vector<RtreeImpl::Leaf> neigh;
    neigh.reserve(_knn ? _knn : 32);

    if (_knn > 0)
      _rtree->tree.query(bgi::nearest(RtreeImpl::Point3D(h.x,h.y,h.z),(int)_knn),std::back_inserter(neigh));
    else
    {
      float w=_stub_window_xy;
      RtreeImpl::Box3D box(RtreeImpl::Point3D(h.x-w,h.y-w,h.z-w),
                           RtreeImpl::Point3D(h.x+w,h.y+w,h.z+w));
      _rtree->tree.query(bgi::intersects(box),std::back_inserter(neigh));
    }
    if (neigh.size()<4) continue;

    std::vector< std::pair<float,float> > rz;
    rz.reserve(neigh.size());
    for (size_t j=0;j<neigh.size();++j)
    {
      const HitXYZ &hn = m_hits[neigh[j].second];
      rz.push_back(std::make_pair(sqrt(hn.x*hn.x+hn.y*hn.y), hn.z));
    }
    if (rz.size()<4) continue;

    float slope=0,z0=0;
    if (!linear_fit_rz(rz,slope,z0)) continue;
    if (fabs(z0) > Z0_MAX) continue; 

#ifdef _OPENMP
    unsigned int tid=(unsigned int)omp_get_thread_num();
#else
    unsigned int tid=0;
#endif
    Stub s; s.slope=slope; s.z0=z0; s.npoints=(int)rz.size();
    thread_stubs[tid].push_back(s);
  }

  for (unsigned int t=0;t<thread_stubs.size();++t)
    m_stubs.insert(m_stubs.end(),thread_stubs[t].begin(),thread_stubs[t].end());

  if (_ntp_stubs)
  {
    for (size_t i=0;i<m_stubs.size();++i)
      _ntp_stubs->Fill((float)_nevent,m_stubs[i].slope,m_stubs[i].z0,(float)m_stubs[i].npoints);
  }
}

void HoughTrackFinder::clusterStubs()
{
  m_tracks.clear();
  if (m_stubs.empty()) return;

  struct Acc { int count; double sumK,sumZ0; };
  std::map<std::pair<int,int>,Acc> bins;
  const float inv_slope_bin=1.f/_slope_bin, inv_z0_bin=1.f/_z0_bin;

  for (size_t i=0;i<m_stubs.size();++i)
  {
    const float k=m_stubs[i].slope, z0=m_stubs[i].z0;
    int ib=(int)floor(k*inv_slope_bin + (k>=0?0.5f:-0.5f));
    int jb=(int)floor(z0*inv_z0_bin + (z0>=0?0.5f:-0.5f));
    std::pair<int,int> key(ib,jb);
    std::map<std::pair<int,int>,Acc>::iterator it=bins.find(key);
    if (it==bins.end())
    { Acc a={1,k,z0}; bins.insert(std::make_pair(key,a)); }
    else { it->second.count++; it->second.sumK+=k; it->second.sumZ0+=z0; }
  }

  typedef std::pair<std::pair<int,int>,Acc> BinEntry;
  std::vector<BinEntry> sorted; sorted.reserve(bins.size());
  for (std::map<std::pair<int,int>,Acc>::iterator it=bins.begin();it!=bins.end();++it)
    sorted.push_back(*it);

  struct BinCompare {
    bool operator()(const BinEntry &a,const BinEntry &b) const {
      if (a.second.count!=b.second.count) return a.second.count>b.second.count;
      if (a.first.first!=b.first.first) return a.first.first<b.first.first;
      return a.first.second<b.first.second;
    }
  };
  std::sort(sorted.begin(),sorted.end(),BinCompare());

  std::set<std::pair<int,int> > used_bins;
  const int dI[5]={0,1,-1,0,0}, dJ[5]={0,0,0,1,-1};

  for (size_t i=0;i<sorted.size();++i)
  {
    const std::pair<int,int> key=sorted[i].first;
    const Acc &acc=sorted[i].second;
    if (acc.count<(int)_min_stubs) continue;
    bool neighbor=false;
    for (int d=0;d<5;++d){
      std::pair<int,int> nb(key.first+dI[d],key.second+dJ[d]);
      if (used_bins.count(nb)){neighbor=true;break;}
    }
    if (neighbor) continue;

    Track tr;
    tr.slope=(float)(acc.sumK/acc.count);
    tr.z0=(float)(acc.sumZ0/acc.count);
    if (fabs(tr.z0) > Z0_MAX) continue;
    tr.nhits=0; tr.rho=0; tr.theta=0;
    m_tracks.push_back(tr);
    used_bins.insert(key);
    if (m_tracks.size()>=_max_tracks) break;
  }
}

void HoughTrackFinder::assignHitsToTracks()
{
  if (m_tracks.empty()) return;
  const unsigned int N=(unsigned int)m_hits.size();

  for (size_t it=0;it<m_tracks.size();++it)
  {
    int assigned=0;
    const float k=m_tracks[it].slope;
    const float z0=m_tracks[it].z0;
    for (unsigned int i=0;i<N;++i)
    {
      const HitXYZ &h=m_hits[i];
      const float r=sqrt(h.x*h.x+h.y*h.y);
      const float zpred=k*r+z0;
      const float dz=fabs(h.z-zpred);
      if (dz<=_max_z_residual) assigned++;
    }
    m_tracks[it].nhits=assigned;
  }

  std::vector<Track> kept;
  kept.reserve(m_tracks.size());
  for (size_t i=0;i<m_tracks.size();++i)
    if (m_tracks[i].nhits>=(int)_min_hits) kept.push_back(m_tracks[i]);
  m_tracks.swap(kept);
}

//============================================================
// Compute (rho, theta) in XY for each track
//============================================================
void HoughTrackFinder::computeXYParameters()
{
  for (size_t it=0; it<m_tracks.size(); ++it)
  {
    std::vector<float> xs, ys;
    xs.reserve(64); ys.reserve(64);
    const float k = m_tracks[it].slope;
    const float z0 = m_tracks[it].z0;

    for (size_t ih=0; ih<m_hits.size(); ++ih)
    {
      const HitXYZ &h = m_hits[ih];
      const float r = sqrt(h.x*h.x + h.y*h.y);
      const float z_pred = k*r + z0;
      const float dz = fabs(h.z - z_pred);
      if (dz > _max_z_residual) continue;
      xs.push_back(h.x); ys.push_back(h.y);
    }
    if (xs.size()<3) continue;

    // simple y = a*x + b fit
    double sx=0,sy=0,sxx=0,sxy=0;
    int n=(int)xs.size();
    for (int i=0;i<n;++i){sx+=xs[i];sy+=ys[i];sxx+=xs[i]*xs[i];sxy+=xs[i]*ys[i];}
    double denom = n*sxx - sx*sx;
    if (fabs(denom)<1e-8) continue;
    double a = (n*sxy - sx*sy)/denom;
    double b = (sy - a*sx)/n;
    double theta = atan(-1.0/a);
    double rho = b*sin(theta);



    if (fabs(m_tracks[it].rho) > RHO_MAX)  m_tracks[it].rho = 0.f;
if (fabs(m_tracks[it].z0)  > Z0_MAX) m_tracks[it].z0  = 0.f;
if (fabs(m_tracks[it].rho) > RHO_MAX || fabs(m_tracks[it].z0) > Z0_MAX) continue;
    m_tracks[it].theta = (float)theta;
    m_tracks[it].rho   = (float)rho;
    
  }
}


void HoughTrackFinder::computeStubXYParameters()
{
  for (size_t i = 0; i < m_stubs.size(); ++i)
  {
    const Stub &s = m_stubs[i];
    // collect nearby hits again for this stub
    std::vector<float> xs, ys;
    xs.reserve(32); ys.reserve(32);

    // For simplicity: just use global hits around z ≈ s.z0 + s.slope * r
    for (size_t ih = 0; ih < m_hits.size(); ++ih)
    {
      const HitXYZ &h = m_hits[ih];
      const float r = sqrt(h.x*h.x + h.y*h.y);
      const float z_pred = s.slope * r + s.z0;
      const float dz = fabs(h.z - z_pred);
      if (dz > _max_z_residual) continue;
      xs.push_back(h.x);
      ys.push_back(h.y);
    }

    if (xs.size() < 3)
    {
      m_stubs[i].rho = 0;
      m_stubs[i].theta = 0;
      continue;
    }

    // Linear regression y = a*x + b
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    int n = (int)xs.size();
    for (int j = 0; j < n; ++j)
    {
      sx += xs[j];
      sy += ys[j];
      sxx += xs[j] * xs[j];
      sxy += xs[j] * ys[j];
    }

    double denom = n*sxx - sx*sx;
    if (fabs(denom) < 1e-8) continue;

    double a = (n*sxy - sx*sy) / denom;
    double b = (sy - a*sx) / n;

    double theta = atan(-1.0/a);
    double rho   = b * sin(theta);
if (fabs(rho) > RHO_MAX)
{
  m_stubs[i].rho = 9999.f; // mark invalid
  continue;
}
    m_stubs[i].theta = (float)theta;
    m_stubs[i].rho   = (float)rho;
  }
}

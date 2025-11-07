#ifndef HOUGH_TRACK_FINDER_H
#define HOUGH_TRACK_FINDER_H

#include <fun4all/SubsysReco.h>
#include <phool/PHCompositeNode.h>

#include <string>
#include <vector>
#include <utility>
#include <cmath>

class TFile;
class TNtuple;

// ------------------------------------------------------------
// Boost R-tree wrapper
// ------------------------------------------------------------
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

struct RtreeImpl
{
  typedef bg::model::point<float, 3, bg::cs::cartesian> Point3D;
  typedef bg::model::box<Point3D> Box3D;
  typedef std::pair<Point3D, unsigned int> Leaf;
  bgi::rtree<Leaf, bgi::quadratic<16> > tree;
};

// ------------------------------------------------------------
// Main class
// ------------------------------------------------------------
class HoughTrackFinder : public SubsysReco
{
public:
  explicit HoughTrackFinder(const std::string &name = "HoughTrackFinder",
                            const std::string &outfile = "hough_output.root");
  ~HoughTrackFinder() override;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;


void set_knn_neighbors(unsigned int k) { _knn = k; }
  // ----------------------------------------------
  // Data structures (must be public for extern use)
  // ----------------------------------------------
  struct HitXYZ
  {
    float x;
    float y;
    float z;
  };

  struct Stub
  {
    float d0, z0, phi0, theta, qOverP;
    int   npoints;
    float slope;
    float rho;
  };

  struct Track
  {
    float slope, z0, rho, theta;
    int   nhits;
  };

private:
  // ---------- Internal helpers ----------
  int  collectHits(PHCompositeNode *topNode);
  void buildRtree();
  void fitLocalStubs();

  // 3D helix fit
  bool helix_fit_3D(const std::vector<RtreeImpl::Leaf> &neigh,
                    const std::vector<HitXYZ> &hits,
                    float &d0, float &z0, float &phi0,
                    float &theta, float &qOverP);

  // fallback straight-line fit
  bool to_straight_perigee(const std::vector<RtreeImpl::Leaf> &neigh,
                           const std::vector<HitXYZ> &hits,
                           float &d0, float &z0, float &phi0,
                           float &theta, float &qOverP);
  void fit_straight_line_xyz( const std::vector<RtreeImpl::Leaf>& neigh,
                            const std::vector<HitXYZ>& hits,
                            double& r0x, double& r0y, double& r0z,
                            double& vx,  double& vy,  double& vz);

  // Optional stubs
  void computeStubXYParameters() {}
  void assignHitsToTracks() {}
  void computeXYParameters() {}

  // ---------- Configuration ----------
  unsigned int _min_layer;
  unsigned int _max_layer;
  unsigned int _knn;
  float _stub_window_r;
  float _stub_window_phi;
  float _stub_window_z;
  float _max_z_residual;
  unsigned int _min_stubs;
  unsigned int _min_hits;
  float _slope_bin;
  float _z0_bin;
  unsigned int _max_tracks;

  // ---------- Output ----------
  std::string m_outputFileName;
  TFile* m_outputFile;
  TNtuple* _ntp_hits;
  TNtuple* _ntp_stubs;
  TNtuple* _ntp_tracks;

  unsigned int _nevent;

  // ---------- Storage ----------
  std::vector<HitXYZ> m_hits;
  std::vector<Stub> m_stubs;
  std::vector<Track> m_tracks;

  RtreeImpl* _rtree;
};

#endif  // HOUGH_TRACK_FINDER_H

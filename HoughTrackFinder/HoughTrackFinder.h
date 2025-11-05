#ifndef HOUGHTRACKFINDER_H
#define HOUGHTRACKFINDER_H

#include <fun4all/SubsysReco.h>
#include <phool/PHCompositeNode.h>

#include <trackbase/TrkrDefs.h>

#include <string>
#include <vector>
#include <map>
#include <cmath>

class TFile;
class TNtuple;

class TrkrHitSetContainer;
class ActsGeometry;
class PHG4TpcGeomContainer;

//============================================================
//  HoughTrackFinder — 3D R-tree–based straight-line finder
//  (local stub fits + parameter-space grid clustering)
//============================================================
class HoughTrackFinder : public SubsysReco
{
 public:
  explicit HoughTrackFinder(const std::string &name = "HoughTrackFinder",
                            const std::string &outfile = "hough_output.root");
  virtual ~HoughTrackFinder();

  // Fun4All hooks
  virtual int Init(PHCompositeNode *topNode);
  virtual int process_event(PHCompositeNode *topNode);
  virtual int End(PHCompositeNode *topNode);

  // -------- Configuration (public setters if you want) --------
  void set_min_layer(unsigned int l) { _min_layer = l; }
  void set_max_layer(unsigned int l) { _max_layer = l; }
  void set_knn_neighbors(unsigned int k) { _knn = k; }
  void set_stub_window_xy(float w) { _stub_window_xy = w; } // used if KNN disabled
  void set_max_z_residual(float z) { _max_z_residual = z; }
  void set_min_stubs(unsigned int n) { _min_stubs = n; }
  void set_min_hits(unsigned int n) { _min_hits = n; }
  void set_slope_bin(float b) { _slope_bin = b; }
  void set_z0_bin(float b) { _z0_bin = b; }
  void set_max_tracks(unsigned int n) { _max_tracks = n; }

 private:
  // ---- Internal DTOs ----
  struct HitXYZ
  {
    float x;
    float y;
    float z;
  };

  struct Stub
  {
  float slope;   // dz/dr
  float z0;      // intercept
  float rho;     // distance from origin in XY
  float theta;   // azimuth of stub line
  int   npoints;
  };

  struct Track
  {

  float slope;
  float z0;
  int   nhits;
  float rho;
  float theta;
  };

  // ---- Pipeline ----
  int  collectHits(PHCompositeNode *topNode);
  void buildRtree();          // build spatial index of hits
  void fitLocalStubs();       // produce (slope,z0) stubs via k-NN
  void clusterStubs();        // grid clustering in parameter space
  void assignHitsToTracks();  // residual assignment in z vs r
  void computeXYParameters();
  void computeStubXYParameters();

  // ---- Utilities ----
  static inline float calc_r(float x, float y)
  { return (float) std::sqrt(x*x + y*y); }

  // Hough parameter windows
static constexpr float RHO_MIN = -5.0f;   // cm
static constexpr float RHO_MAX =  +5.0f;  // cm
static constexpr float Z0_MIN  = -10.0f;  // cm
static constexpr float Z0_MAX  = +10.0f;  // cm

// Example binning (tune as you like)
static constexpr int   NBINS_RHO = 200;       // 0.05 cm/bin over 10 cm
static constexpr int   NBINS_Z0  = 200;       // 0.10 cm/bin over 20 cm
static constexpr int   NBINS_PHI0 = 360;      // for XY curvature orientation, if used


  // ---- Configuration ----
  unsigned int _min_layer;         // inclusive (TPC layer range filter)
  unsigned int _max_layer;         // inclusive
  unsigned int _knn;               // number of nearest neighbors (if 0 -> box query)
  float        _stub_window_xy;    // used only for box query fallback (cm)
  float        _max_z_residual;    // assignment residual in z (cm)
  unsigned int _min_stubs;         // minimum stubs in a parameter bin to form a track
  unsigned int _min_hits;          // minimum assigned hits to keep a track
  float        _slope_bin;         // bin size for dz/dr
  float        _z0_bin;            // bin size for z0
  unsigned int _max_tracks;        // cap on produced tracks per event

  // ---- State / I/O ----
  std::string m_outputFileName;
  TFile   *m_outputFile;
  TNtuple *_ntp_hits;    // event:x:y:z
  TNtuple *_ntp_stubs;   // event:slope:z0:npts
  TNtuple *_ntp_tracks;  // event:slope:z0:nhits

  unsigned int _nevent;

  // Data
  std::vector<HitXYZ> m_hits;
  std::vector<Stub>   m_stubs;
  std::vector<Track>  m_tracks;

  // R-tree storage (opaque here; defined in .cc)
  struct RtreeImpl;
  RtreeImpl *_rtree; // pointer to avoid including boost headers in the header
};

#endif

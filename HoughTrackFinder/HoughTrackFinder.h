#ifndef HoughTrackFinder_H
#define HoughTrackFinder_H
#include <TH2F.h>
#include <TTree.h>


//begin

#include "trackreco/PHTrackSeeding.h"

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>       // for cluskey
#include <trackbase/ActsGeometry.h>

// TrkrCluster includes
#include <trackbase/TrkrCluster.h>    // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHRandomSeed.h>
#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>              // for PHWHERE

// ROOT includes for debugging
#include <TFile.h>
#include <TMatrixDSymfwd.h>           // for TMatrixDSym
#include <TMatrixTSym.h>              // for TMatrixTSym
#include <TMatrixTUtils.h>
#include <TNtuple.h>
#include <TVector3.h>
#include <TVectorDfwd.h>
#include <TVectorT.h>

// gsl
#include <gsl/gsl_rng.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>

// BOOST for combi seeding
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

// standard includes
#include <algorithm>
#include <cassert>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

class TGeoManager;

#define LogDebug(exp)   std::cout << "DEBUG: "   << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogError(exp)   std::cout << "ERROR: "   << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp

using namespace Eigen;
using namespace std;

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

// forward declarations
class PHCompositeNode;
class PHG4CellContainer;
class PHG4CylinderGeomContainer;
class PHG4HitContainer;
class PHTimer;
class sPHENIXSeedFinder;
class SvtxClusterMap;
class SvtxCluster;
class SvtxTrackMap;
class SvtxTrack;
class SvtxVertexMap;
class SvtxVertex;
class TNtuple;
class TFile;
class TRKR_CLUSTER;
class SvtxHitMap;
class TrackSeedContainer;

// geometry / rtree typedefs (matching the .cc file)
using point    = bg::model::point<float, 3, bg::cs::cartesian>;
using box      = bg::model::box<point>;
using pointKey = std::pair<point, TrkrDefs::cluskey>;
using myrtree  = bgi::rtree<pointKey, bgi::quadratic<16>>;

using cluskey_t = uint64_t;

class HoughTrackFinder : public SubsysReco
{
 public:
  HoughTrackFinder(const std::string &name = "PHRTreeSeeding");

  double chisq(const double *xx);

  virtual ~HoughTrackFinder() {}

  void set_write_debug_ntuple(bool b)   { _write_ntp         = b; }
  void set_create_tracks(bool b)        { _create_tracks     = b; }
  void set_max_distance_to_origin(float val) { _max_dist_to_origin = val; }
  void set_min_nclusters(int n)         { _min_nclusters     = n; }

 protected:
  int Setup(PHCompositeNode *topNode);
  int GetNodes(PHCompositeNode* topNode);
  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:
  // cluster container
  TrkrClusterContainer *_cluster_map = nullptr;

  // geometry
  ActsGeometry *tGeometry {nullptr};

  // helpers
  double phiadd(double phi1, double phi2);
  double phidiff(double phi1, double phi2);
  double pointKeyToTuple(pointKey *pK);
  double costfunction(const double *xx);

  // stub finder: still fits y = m x + b in XY and returns slope/intercept
void get_stub(const myrtree &search_rtree,
              float pointr, float pointphi, float pointz,
              int &count,
              double &dca,
              double &phi,
              double &z0,
              double &theta);


#ifndef __CINT__
 private:
  int createNodes(PHCompositeNode *topNode);

  struct HoughPeak
  {
    float dca;
    float phi;
    float theta;
  };

  std::string        m_trackMapName   = "TpcTrackSeedContainer";
  TrackSeedContainer *m_seedContainer = nullptr;

  unsigned int _nevent            = 0;
  bool         _write_ntp         = true;
  bool         _create_tracks     = true;
  float        _max_dist_to_origin = 0;
  unsigned int _min_nclusters     = 20;

  TNtuple *_ntp_cos  = nullptr;  // cluster positions
  TNtuple *_ntp_stub = nullptr;  // stub Hough params (now dca, phi)
  TNtuple *_ntp_max  = nullptr;  // Hough ranges
  TNtuple *_ntp_trk  = nullptr;  // track params
   TH2F    *_hHough   = nullptr;

     TTree* _track_tree;

  int   _trk_ev;
  float _trk_dca;
  float _trk_phi;
  float _trk_alpha;
  float _trk_dist0;
  int   _trk_nclus;

  std::vector<float>    _trk_r;       
  std::vector<float>    _trk_phi_clu; 
  std::vector<float>    _trk_z;      
  std::vector<ULong64_t> _trk_cluskey;

  TFile   *_tfile    = nullptr;

#endif  // __CINT__
};

#endif

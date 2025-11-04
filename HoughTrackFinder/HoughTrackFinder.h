#ifndef HOUGHTRACKFINDER_H
#define HOUGHTRACKFINDER_H

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TpcDefs.h>
#include <string>
#include <vector>
#include <map>

class PHCompositeNode;
class TrkrHitSetContainer;
class PHG4TpcGeomContainer;
class PHG4CylinderGeomContainer;
class TH2F;
class TFile;
class TTree;

class HoughTrackFinder : public SubsysReco
{
 public:
  explicit HoughTrackFinder(const std::string &name = "HoughTrackFinder", const std::string &filename = "hough_tracks.root");
  ~HoughTrackFinder() override {}

  int Init(PHCompositeNode *) override;
  int process_event(PHCompositeNode *) override;
  int End(PHCompositeNode *) override;

 private:
  struct HitXY {
    double x;
    double y;
  };

  // Internal helper functions
  void fillHough(const std::vector<HitXY> &hits);
  void findPeaks();
  void assignHitsToTracks(const std::vector<HitXY> &hits);
  void clearEvent();


  std::string m_outputFileName;
  TFile *m_outputFile;
  TTree *m_tree;


  // Parameters
  int n_theta_bins = 180;
  int n_rho_bins   = 200;
  double max_rho   = 100.0; // cm
  double peak_threshold_fraction = 0.5;
  double max_hit_distance = 0.5; // cm

  // Hough accumulator
  TH2F* h_hough = nullptr;
  float m_rho = -9999, m_theta = -9999;
  int m_nhits = 0;

  // Found peaks (θ,ρ)
  std::vector<std::pair<double,double>> track_peaks;

  // Hits per track
  std::map<int, std::vector<HitXY>> track_hits;
};

#endif

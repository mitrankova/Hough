#ifndef HOUGH_TRACK_FINDER_H
#define HOUGH_TRACK_FINDER_H

#include <fun4all/SubsysReco.h>
#include <phool/PHCompositeNode.h>

#include <TFile.h>
#include <TH2F.h>
#include <TTree.h>

#include <vector>
#include <map>
#include <string>
#include <utility>

//============================================================
//  HoughTrackFinder
//  3D Hough Transform for straight-line track reconstruction
//  Step 1: 2D Hough (θ,ρ) in XY plane
//  Step 2: 2D Hough (k,z0) in Z–R plane for each XY line
//============================================================

class HoughTrackFinder : public SubsysReco
{
public:
  struct HitXYZ
  {
    double x;
    double y;
    double z;
  };

  explicit HoughTrackFinder(const std::string &name = "HoughTrackFinder",
                            const std::string &filename = "hough3d_output.root");
  ~HoughTrackFinder() override = default;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

private:
  // ---- Core algorithms ----
  void fillHough(const std::vector<HitXYZ> &hits);              // Fill XY accumulator
  void findPeaks();                                             // Find XY peaks
  void assignHitsToTracks(const std::vector<HitXYZ> &hits);     // Assign hits to lines
  void houghKZForTrack(size_t itrk);                            // Per-track (k,z0) Hough
  std::pair<float, float> findKZPeak(TH2F *acc, float frac);    // Find (k,z0) max
  double xyDistanceToLine(double x, double y,
                          double theta, double rho) const;      // Distance from hit to line
  void clearEvent();                                            // Reset per-event memory

  // ---- Output ----
  std::string m_outputFileName;
  TFile *m_outputFile = nullptr;
  TH2F *h_hough = nullptr;  // XY Hough accumulator
  TH2F *h_kz = nullptr;     // per-track (k,z0) accumulator
  TTree *m_tree = nullptr;

  // ---- Per-event variables ----
  int   m_nhits = 0;
  float m_rho = 0;
  float m_theta = 0;
  float m_slope_dzdr = 0;
  float m_phi = 0;
  float m_z0 = 0;

  // ---- Parameters ----
  static constexpr int n_theta_bins = 360;
  static constexpr int n_rho_bins   = 200;
  static constexpr double max_rho   = 5.0;   // ±5 cm
  static constexpr double peak_threshold_fraction = 0.5;
  static constexpr double assign_xy_maxdist = 0.25;  // cm

  // ---- (k,z0) Hough parameters ----
  static constexpr int n_k_bins   = 160;
  static constexpr float k_min    = -3.0f;
  static constexpr float k_max    = 3.0f;
  static constexpr int n_z0_bins  = 200;
  static constexpr float z0_min   = -10.0f;  // ±10 cm
  static constexpr float z0_max   = 10.0f;   // ±10 cm

  // ---- Internal storage ----
  std::vector<std::pair<double, double>> track_peaks;  // (θ,ρ)
  std::map<int, std::vector<HitXYZ>> track_hits;       // hits per XY line
  std::vector<std::pair<float, float>> track_kz;       // (k,z0) per line
};

#endif  // HOUGH_TRACK_FINDER_H

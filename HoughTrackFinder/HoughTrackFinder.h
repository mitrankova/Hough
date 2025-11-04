#ifndef HOUGH_TRACK_FINDER_H
#define HOUGH_TRACK_FINDER_H

#include <fun4all/SubsysReco.h>
#include <phool/PHCompositeNode.h>

#include <vector>
#include <map>
#include <string>
#include <utility>
#include <TFile.h>
#include <TH2F.h>
#include <TTree.h>

class HoughTrackFinder : public SubsysReco
{
public:
  struct HitXY
  {
    double x;
    double y;
  };

  explicit HoughTrackFinder(const std::string &name = "HoughTrackFinder",
                            const std::string &filename = "hough_output.root");

  ~HoughTrackFinder() override = default;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

private:
  void fillHough(const std::vector<HitXY> &hits);
  void findPeaks();
  void assignHitsToTracks(const std::vector<HitXY> &hits);
  void clearEvent();

  // Output
  std::string m_outputFileName;
  TFile *m_outputFile = nullptr;
  TH2F *h_hough = nullptr;
  TTree *m_tree = nullptr;

  // Event variables
  int m_nhits = 0;
  float m_rho = 0;
  float m_theta = 0;

  // Parameters
  static constexpr int n_theta_bins = 360;
  static constexpr int n_rho_bins   = 200;
  static constexpr double max_rho   = 5.0;             // ±5 cm constraint
  static constexpr double peak_threshold_fraction = 0.5;
  static constexpr double max_hit_distance = 0.3;      // cm tolerance when grouping

  // Internal storage
  std::vector<std::pair<double, double>> track_peaks;  // (θ, ρ)
  std::map<int, std::vector<HitXY>> track_hits;
};

#endif

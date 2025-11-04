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

#include <cmath>
#include <iostream>

//============================================================
// Constructor
//============================================================
HoughTrackFinder::HoughTrackFinder(const std::string &name, const std::string &filename)
  : SubsysReco(name)
  , m_outputFileName(filename)
{}

//============================================================
// Init
//============================================================
int HoughTrackFinder::Init(PHCompositeNode* /*topNode*/)
{
  std::cout << "HoughTrackFinder::Init - Creating output file: " << m_outputFileName << std::endl;

  m_outputFile = new TFile(m_outputFileName.c_str(), "RECREATE");
  if (!m_outputFile || m_outputFile->IsZombie())
  {
    std::cerr << "HoughTrackFinder::Init - ERROR: Cannot create output file " << m_outputFileName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Only ±5 cm around origin
  h_hough = new TH2F("h_hough", "Hough accumulator;#theta [rad];#rho [cm]",
                     n_theta_bins, -M_PI/2, M_PI/2,
                     n_rho_bins, -max_rho, max_rho);

  m_tree = new TTree("phase_tree", "Hough Line Reconstruction");
  m_tree->Branch("nhits", &m_nhits, "m_nhits/I");
  m_tree->Branch("rho", &m_rho, "m_rho/F");
  m_tree->Branch("theta", &m_theta, "m_theta/F");

  return Fun4AllReturnCodes::EVENT_OK;
}

//============================================================
// process_event
//============================================================
int HoughTrackFinder::process_event(PHCompositeNode* topNode)
{
  auto *hitmap = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!hitmap)
  {
    std::cerr << "HoughTrackFinder::process_event - No TRKR_HITSET node found!" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  auto *tpcGeom = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  auto *geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!tpcGeom || !geometry)
  {
    std::cerr << "HoughTrackFinder::process_event - Missing geometry nodes!" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  clearEvent();

  std::vector<HitXY> hits;

  //------------------------------------------------------------
  // 1. Collect (x,y) positions from all TPC hits
  //------------------------------------------------------------
  TrkrHitSetContainer::ConstRange all_hitsets = hitmap->getHitSets();
  for (auto hitsetiter = all_hitsets.first; hitsetiter != all_hitsets.second; ++hitsetiter)
  {
    auto hitsetkey = hitsetiter->first;
    TrkrHitSet* hitset = hitsetiter->second;
    if (TrkrDefs::getTrkrId(hitsetkey) != TrkrDefs::tpcId)
      continue;

    unsigned int layer = TrkrDefs::getLayer(hitsetkey);
    unsigned int side = TpcDefs::getSide(hitsetkey);
    auto *geoLayer = tpcGeom->GetLayerCellGeom(layer);

    for (auto hitr = hitset->getHits().first; hitr != hitset->getHits().second; ++hitr)
    {
      auto hitkey = hitr->first;
      auto pad = TpcDefs::getPad(hitkey);
      double phi = geoLayer->get_phicenter(pad, side);
      double radius = geoLayer->get_radius();
      float AdcClockPeriod = geoLayer->get_zstep();
      auto glob = geometry->getGlobalPositionTpc(hitsetkey, hitkey, phi, radius, AdcClockPeriod);

      hits.push_back({glob.x(), glob.y()});
    }
  }

  if (hits.empty())
  {
    std::cout << "HoughTrackFinder::process_event - No hits found in event." << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  //------------------------------------------------------------
  // 2. Fill Hough accumulator
  //------------------------------------------------------------
  fillHough(hits);

  //------------------------------------------------------------
  // 3. Find peaks
  //------------------------------------------------------------
  findPeaks();

  //------------------------------------------------------------
  // 4. Assign hits to tracks
  //------------------------------------------------------------
  assignHitsToTracks(hits);

  //------------------------------------------------------------
  // 5. Output summary
  //------------------------------------------------------------
  std::cout << "HoughTrackFinder: found " << track_peaks.size() << " track candidates" << std::endl;
  for (size_t i = 0; i < track_peaks.size(); ++i)
  {
    m_nhits = track_hits[i].size();
    m_theta = track_peaks[i].first;
    m_rho   = track_peaks[i].second;
    m_tree->Fill();

    std::cout << "  Track " << i
              << " θ=" << m_theta << " ρ=" << m_rho
              << " with " << m_nhits << " hits" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//============================================================
// Fill Hough accumulator (θ, ρ) within ±5 cm constraint
//============================================================
void HoughTrackFinder::fillHough(const std::vector<HitXY> &hits)
{
  h_hough->Reset();

  for (const auto &hit : hits)
  {
    for (int itheta = 0; itheta < n_theta_bins; ++itheta)
    {
      double theta = -M_PI/2 + itheta * M_PI / n_theta_bins;
      double rho = hit.x * std::cos(theta) + hit.y * std::sin(theta);

      // Only fill if line passes within ±5 cm of origin
      if (std::fabs(rho) < max_rho)
        h_hough->Fill(theta, rho);
    }
  }
}

//============================================================
// Find peaks (local maxima) above threshold
//============================================================
void HoughTrackFinder::findPeaks()
{
  track_peaks.clear();

  int ntheta = h_hough->GetNbinsX();
  int nrho   = h_hough->GetNbinsY();
  double maxvote = h_hough->GetMaximum();
  double threshold = peak_threshold_fraction * maxvote;

  for (int ix = 2; ix < ntheta - 1; ++ix)
  {
    for (int iy = 2; iy < nrho - 1; ++iy)
    {
      double val = h_hough->GetBinContent(ix, iy);
      if (val < threshold) continue;

      bool isLocalMax = true;
      for (int dx = -1; dx <= 1 && isLocalMax; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
        {
          if (dx == 0 && dy == 0) continue;
          if (h_hough->GetBinContent(ix + dx, iy + dy) > val)
          {
            isLocalMax = false;
            break;
          }
        }

      if (isLocalMax)
      {
        double theta = h_hough->GetXaxis()->GetBinCenter(ix);
        double rho   = h_hough->GetYaxis()->GetBinCenter(iy);
        track_peaks.push_back({theta, rho});
      }
    }
  }
}

//============================================================
// Assign hits to nearest line candidate
//============================================================
void HoughTrackFinder::assignHitsToTracks(const std::vector<HitXY> &hits)
{
  track_hits.clear();

  for (const auto &hit : hits)
  {
    double bestDist = 1e9;
    int bestTrack = -1;

    for (size_t i = 0; i < track_peaks.size(); ++i)
    {
      auto [theta, rho] = track_peaks[i];
      double d = std::fabs(hit.x * std::cos(theta) + hit.y * std::sin(theta) - rho);
      if (d < bestDist)
      {
        bestDist = d;
        bestTrack = i;
      }
    }

    // Keep only hits within a narrow band around each line
    if (bestTrack >= 0 && bestDist < max_hit_distance)
      track_hits[bestTrack].push_back(hit);
  }
}

//============================================================
// Clear event
//============================================================
void HoughTrackFinder::clearEvent()
{
  track_peaks.clear();
  track_hits.clear();
}

//============================================================
// End of job
//============================================================
int HoughTrackFinder::End(PHCompositeNode* /*topNode*/)
{
  m_outputFile->cd();
  h_hough->Write();
  m_tree->Write();
  m_outputFile->Close();
  return Fun4AllReturnCodes::EVENT_OK;
}

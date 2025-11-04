#include "HoughTrackFinder.h"

#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TpcDefs.h>

#include <trackbase/ActsGeometry.h>  
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <TH2F.h>
#include <TFile.h>


#include <TTree.h>


#include <cmath>
#include <iostream>
#include <algorithm>

HoughTrackFinder::HoughTrackFinder(const std::string &name, const std::string &filename)
  : SubsysReco(name)
  , m_outputFileName(filename)
  , m_outputFile(nullptr)
  , m_tree(nullptr)
{}

int HoughTrackFinder::Init(PHCompositeNode* /*topNode*/)
{
std::cout << "ClusterPhaseAnalysis::Init - Creating output file: " << m_outputFileName << std::endl;

  m_outputFile = new TFile(m_outputFileName.c_str(), "RECREATE");
  if (!m_outputFile || m_outputFile->IsZombie())
  {
    std::cout << "ClusterPhaseAnalysis::Init - Error: Cannot create output file" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  h_hough = new TH2F("h_hough","Hough accumulator;theta [rad];rho [cm]",
                     n_theta_bins, -M_PI/2, M_PI/2,
                     n_rho_bins, -max_rho, max_rho);

  m_tree = new TTree("phase_tree", "Cluster Phase Analysis Tree");

  // Event-level branches
  m_tree->Branch("nhits", &m_nhits, "m_nhits/I");
  m_tree->Branch("rho", &m_rho, "m_rho/F");
  m_tree->Branch("theta", &m_theta, "m_theta/F");

  return Fun4AllReturnCodes::EVENT_OK;
}

int HoughTrackFinder::process_event(PHCompositeNode* topNode)
{
  TrkrHitSetContainer* hitmap =
      findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!hitmap)
  {
    std::cout << "No TRKR_HITSET node found!" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  auto *tpcGeom = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  auto *geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  std::cout << "HoughTrackFinder::process_event - Processing event" << std::endl;
  clearEvent();

  std::vector<HitXY> hits;

  // ------------------------------------------------------------------
  // 1. Collect hits (replace getX/getY if you have to decode layer/phi)
  // ------------------------------------------------------------------
  TrkrHitSetContainer::ConstRange all_hitsets = hitmap->getHitSets();
  for (TrkrHitSetContainer::ConstIterator hitsetiter = all_hitsets.first;
       hitsetiter != all_hitsets.second;
       ++hitsetiter)
  {
    auto m_hitsetkey = hitsetiter->first;
    TrkrHitSet* hitset = hitsetiter->second;

    auto m_hitlayer = TrkrDefs::getLayer(m_hitsetkey);
    auto m_side = TpcDefs::getSide(m_hitsetkey);
    std::cout<<"HoughTrackFinder::process_event --- Processing hitset key "<<m_hitsetkey
             <<" layer "<<m_hitlayer<<" side "<<m_side<<std::endl;

    auto det = TrkrDefs::getTrkrId(m_hitsetkey);
    if (det != TrkrDefs::tpcId) continue;

    auto hitrangei = hitset->getHits();
    for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
         hitr != hitrangei.second;
         ++hitr)
    {
      auto hitkey = hitr->first;
     // auto *hit = hitr->second;
     // m_adc = hit->getAdc();
         auto    m_hitpad = TpcDefs::getPad(hitkey);

        auto *geoLayer = tpcGeom->GetLayerCellGeom(m_hitlayer);
        auto phi = geoLayer->get_phicenter(m_hitpad, m_side);
        auto radius = geoLayer->get_radius();
        float AdcClockPeriod = geoLayer->get_zstep();
        auto glob = geometry->getGlobalPositionTpc(m_hitsetkey, hitkey, phi, radius, AdcClockPeriod);


      double x = glob.x();
      double y = glob.y();

      hits.push_back({x, y});
    }
  }
  

  if (hits.empty())
  {
    std::cout << "No hits found in event." << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // ------------------------------------------------------------------
  // 2. Fill Hough accumulator
  // ------------------------------------------------------------------
  fillHough(hits);

  // ------------------------------------------------------------------
  // 3. Find peaks in accumulator (possible tracks)
  // ------------------------------------------------------------------
  findPeaks();

  // ------------------------------------------------------------------
  // 4. Assign hits to each track candidate
  // ------------------------------------------------------------------
  assignHitsToTracks(hits);

  // ------------------------------------------------------------------
  // 5. Print summary
  // ------------------------------------------------------------------
  std::cout << "Found " << track_peaks.size() << " track candidates." << std::endl;
  for (size_t i = 0; i < track_peaks.size(); ++i)
  {
    std::cout << "  Track " << i
              << " theta=" << track_peaks[i].first
              << " rho=" << track_peaks[i].second
              << " with " << track_hits[i].size() << " hits." << std::endl;
    m_nhits = track_hits[i].size();
    m_theta = track_peaks[i].first;
    m_rho   = track_peaks[i].second;
    m_tree->Fill();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void HoughTrackFinder::fillHough(const std::vector<HitXY> &hits)
{
  h_hough->Reset();
std::cout << "Filling Hough accumulator with " << hits.size() << " hits." << std::endl;
  for (const auto &hit : hits)
  {
    for (int itheta = 0; itheta < n_theta_bins; ++itheta)
    {
      double theta = -M_PI/2 + itheta * M_PI / n_theta_bins;
      double rho = hit.x * std::cos(theta) + hit.y * std::sin(theta);
      h_hough->Fill(theta, rho);
    }
  }
}

void HoughTrackFinder::findPeaks()
{
    std::cout << "Finding peaks in Hough accumulator." << std::endl;
  track_peaks.clear();

  int ntheta = h_hough->GetNbinsX();
  int nrho   = h_hough->GetNbinsY();
  double maxvote = h_hough->GetMaximum();
  double threshold = peak_threshold_fraction * maxvote;

  for (int ix = 2; ix < ntheta; ++ix)
  {
    for (int iy = 2; iy < nrho; ++iy)
    {
      double val = h_hough->GetBinContent(ix, iy);
      if (val < threshold) continue;

      bool isLocalMax = true;
      for (int dx = -1; dx <= 1; ++dx)
      {
        for (int dy = -1; dy <= 1; ++dy)
        {
          if (dx == 0 && dy == 0) continue;
          if (h_hough->GetBinContent(ix+dx, iy+dy) > val)
          {
            isLocalMax = false;
            break;
          }
        }
        if (!isLocalMax) break;
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

void HoughTrackFinder::assignHitsToTracks(const std::vector<HitXY> &hits)
{
    std::cout << "Assigning hits to tracks." << std::endl;
  track_hits.clear();

  for (const auto &hit : hits)
  {
    double bestDist = 1e9;
    int bestTrack = -1;
    for (size_t i=0; i<track_peaks.size(); ++i)
    {
      auto [theta, rho] = track_peaks[i];
      double d = std::fabs(hit.x * std::cos(theta) + hit.y * std::sin(theta) - rho);
      if (d < bestDist) { bestDist = d; bestTrack = i; }
    }
    if (bestTrack >= 0 && bestDist < max_hit_distance)
      track_hits[bestTrack].push_back(hit);
  }
}

void HoughTrackFinder::clearEvent()
{
    std::cout << "Clearing event data." << std::endl;
  track_peaks.clear();
  track_hits.clear();
}

int HoughTrackFinder::End(PHCompositeNode* /*topNode*/)
{
  m_outputFile->cd();
  h_hough->Write();
  m_tree->Write();
  m_outputFile->Close();
  return Fun4AllReturnCodes::EVENT_OK;
}

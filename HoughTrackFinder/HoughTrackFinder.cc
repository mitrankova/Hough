//============================================================
// HoughTrackFinder.cc  — 3D straight-line Hough (r,φ,z)
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

#include <cmath>
#include <iostream>
#include <algorithm>

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
  std::cout << "HoughTrackFinder::Init - Creating output file: "
            << m_outputFileName << std::endl;

  m_outputFile = new TFile(m_outputFileName.c_str(), "RECREATE");
  if (!m_outputFile || m_outputFile->IsZombie())
  {
    std::cerr << "ERROR: cannot create " << m_outputFileName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  h_hough = new TH2F("h_hough", "Hough XY;#theta [rad];#rho [cm]",
                     n_theta_bins, -M_PI/2, M_PI/2,
                     n_rho_bins, -max_rho, max_rho);

  m_tree = new TTree("phase_tree", "3D Hough Line Reconstruction");
  m_tree->Branch("nhits", &m_nhits, "nhits/I");
  m_tree->Branch("rho", &m_rho, "rho/F");
  m_tree->Branch("theta", &m_theta, "theta/F");
  m_tree->Branch("dzdr", &m_slope_dzdr, "dzdr/F");
  m_tree->Branch("phi", &m_phi, "phi/F");
  m_tree->Branch("z0", &m_z0, "z0/F");

  return Fun4AllReturnCodes::EVENT_OK;
}

//============================================================
// process_event
//============================================================
int HoughTrackFinder::process_event(PHCompositeNode* topNode)
{
  auto *hitmap = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!hitmap) return Fun4AllReturnCodes::ABORTEVENT;

  auto *tpcGeom = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  auto *geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!tpcGeom || !geometry) return Fun4AllReturnCodes::ABORTEVENT;

  clearEvent();
  std::vector<HitXYZ> hits;

  //------------------------------------------------------------
  // 1. Collect hits (x,y,z)
  //------------------------------------------------------------
  for (auto hitsetiter = hitmap->getHitSets().first;
       hitsetiter != hitmap->getHitSets().second; ++hitsetiter)
  {
    auto hitsetkey = hitsetiter->first;
    TrkrHitSet* hitset = hitsetiter->second;
    if (TrkrDefs::getTrkrId(hitsetkey) != TrkrDefs::tpcId) continue;

    unsigned int layer = TrkrDefs::getLayer(hitsetkey);
    unsigned int side  = TpcDefs::getSide(hitsetkey);
    auto *geoLayer = tpcGeom->GetLayerCellGeom(layer);

    for (auto hitr = hitset->getHits().first; hitr != hitset->getHits().second; ++hitr)
    {
      auto hitkey = hitr->first;
      auto pad = TpcDefs::getPad(hitkey);
      double phi = geoLayer->get_phicenter(pad, side);
      double radius = geoLayer->get_radius();
      float AdcClockPeriod = geoLayer->get_zstep();

      auto glob = geometry->getGlobalPositionTpc(hitsetkey, hitkey, phi, radius, AdcClockPeriod);
      hits.push_back({glob.x(), glob.y(), glob.z()});
    }
  }

  if (hits.empty()) return Fun4AllReturnCodes::EVENT_OK;

  //------------------------------------------------------------
  // 2D Hough in XY
  //------------------------------------------------------------
  fillHough(hits);
  findPeaks();
  assignHitsToTracks(hits);

  //------------------------------------------------------------
  // 3D: (k,z0) per XY line
  //------------------------------------------------------------
  track_kz.clear();
  for (size_t i = 0; i < track_peaks.size(); ++i) houghKZForTrack(i);

  //------------------------------------------------------------
  // Save
  //------------------------------------------------------------
  for (size_t i = 0; i < track_peaks.size(); ++i)
  {
    m_nhits = track_hits[i].size();
    if (m_nhits < 3) continue;

    m_theta = track_peaks[i].first;
    m_rho   = track_peaks[i].second;

    const auto kz = (i < track_kz.size()) ? track_kz[i] : std::pair<float,float>{0.f,0.f};
    m_slope_dzdr = kz.first;
    m_z0         = kz.second;
    m_phi        = std::atan(m_slope_dzdr);

    m_tree->Fill();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//============================================================
// fillHough
//============================================================
void HoughTrackFinder::fillHough(const std::vector<HitXYZ> &hits)
{
  h_hough->Reset();
  for (const auto &h : hits)
  {
    for (int it = 0; it < n_theta_bins; ++it)
    {
      double theta = -M_PI/2 + it * M_PI / n_theta_bins;
      double rho = h.x * std::cos(theta) + h.y * std::sin(theta);
      if (std::fabs(rho) <= max_rho) h_hough->Fill(theta, rho);
    }
  }
}

//============================================================
// findPeaks in XY
//============================================================
void HoughTrackFinder::findPeaks()
{
  track_peaks.clear();
  const int nx = h_hough->GetNbinsX();
  const int ny = h_hough->GetNbinsY();
  const double maxvote = h_hough->GetMaximum();
  const double thr = peak_threshold_fraction * maxvote;

  for (int ix=2; ix<nx-1; ++ix)
  for (int iy=2; iy<ny-1; ++iy)
  {
    const double v = h_hough->GetBinContent(ix,iy);
    if (v<thr) continue;
    bool localmax=true;
    for (int dx=-1; dx<=1 && localmax; ++dx)
    for (int dy=-1; dy<=1; ++dy)
    {
      if(dx==0 && dy==0) continue;
      if(h_hough->GetBinContent(ix+dx,iy+dy)>v) localmax=false;
    }
    if(!localmax) continue;

    double th = h_hough->GetXaxis()->GetBinCenter(ix);
    double rh = h_hough->GetYaxis()->GetBinCenter(iy);
    if (std::fabs(rh) <= max_rho) track_peaks.push_back({th, rh});
  }
}

//============================================================
// assignHitsToTracks
//============================================================
void HoughTrackFinder::assignHitsToTracks(const std::vector<HitXYZ> &hits)
{
  track_hits.clear();
  for(size_t i=0;i<track_peaks.size();++i)
    track_hits[(int)i]={};

  for(const auto &h:hits)
  {
    int best=-1; double bestd=1e9;
    for(size_t i=0;i<track_peaks.size();++i)
    {
      const auto &p=track_peaks[i];
      double d=xyDistanceToLine(h.x,h.y,p.first,p.second);
      if(d<bestd){bestd=d;best=(int)i;}
    }
    if(best>=0 && bestd<=assign_xy_maxdist)
      track_hits[best].push_back(h);
  }
}

double HoughTrackFinder::xyDistanceToLine(double x,double y,double theta,double rho) const
{
  return std::fabs(x*std::cos(theta)+y*std::sin(theta)-rho);
}

//============================================================
// 2D Hough in (k,z0) for each XY line
//============================================================
void HoughTrackFinder::houghKZForTrack(size_t itrk)
{
  if(h_kz){delete h_kz; h_kz=nullptr;}
  h_kz=new TH2F("h_kz","k-z0 accumulator;k=dz/dr;z0[cm]",
                n_k_bins,k_min,k_max,n_z0_bins,z0_min,z0_max);

  const auto &hits=track_hits[(int)itrk];
  if(hits.size()<3) return;

  std::vector<std::pair<double,double>> rz;
  rz.reserve(hits.size());
  for(const auto &h:hits)
    rz.emplace_back(std::sqrt(h.x*h.x+h.y*h.y),h.z);

  for(const auto &p:rz)
  {
    const double r=p.first;
    const double z=p.second;
    for(int ik=1;ik<=n_k_bins;++ik)
    {
      const double k=h_kz->GetXaxis()->GetBinCenter(ik);
      const double z0=z-k*r;
      if(std::fabs(z0)<=10.0) h_kz->Fill(k,z0);
    }
  }

  auto best=findKZPeak(h_kz,0.5);
  if(itrk>=track_kz.size()) track_kz.resize(itrk+1,{0.f,0.f});
  track_kz[itrk]=best;
}

std::pair<float,float> HoughTrackFinder::findKZPeak(TH2F* acc,float frac)
{
  if(!acc) return {0.f,0.f};
  const double maxv=acc->GetMaximum();
  const double thr=std::max(1.0,frac*maxv);
  int nx=acc->GetNbinsX(),ny=acc->GetNbinsY();
  double bestv=-1; int bix=1,biy=1;

  for(int ix=2;ix<nx;++ix)
  for(int iy=2;iy<ny;++iy)
  {
    double v=acc->GetBinContent(ix,iy);
    if(v<thr) continue;
    bool loc=true;
    for(int dx=-1;dx<=1&&loc;++dx)
    for(int dy=-1;dy<=1;++dy)
    {
      if(dx==0&&dy==0) continue;
      if(acc->GetBinContent(ix+dx,iy+dy)>v) loc=false;
    }
    if(loc&&v>bestv){bestv=v;bix=ix;biy=iy;}
  }

  float k=acc->GetXaxis()->GetBinCenter(bix);
  float z0=acc->GetYaxis()->GetBinCenter(biy);
  return {k,z0};
}

//============================================================
// clearEvent / End
//============================================================
void HoughTrackFinder::clearEvent()
{
  if(h_hough)h_hough->Reset("ICES");
  if(h_kz)h_kz->Reset("ICES");
  track_peaks.clear(); track_hits.clear(); track_kz.clear();
  m_nhits=0; m_rho=0; m_theta=0; m_slope_dzdr=0; m_phi=0; m_z0=0;
}

int HoughTrackFinder::End(PHCompositeNode*)
{
  if(m_outputFile)
  {
    m_outputFile->cd();
    if(h_hough)h_hough->Write();
    if(m_tree)m_tree->Write();
    m_outputFile->Close();
    delete m_outputFile; m_outputFile=nullptr;
  }
  if(h_kz){delete h_kz;h_kz=nullptr;}
  return Fun4AllReturnCodes::EVENT_OK;
}

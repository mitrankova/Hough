#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLine.h>
#include <TStyle.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <cmath>

void PlotTracks() {
  gStyle->SetOptStat(0);

  // --- Open input files ---
  TFile *f_resid = TFile::Open("/sphenix/user/mitrankova/Hough/output_resid/TPC_Au_Au_ZeroField_1mrad_1evt_skip8_aligned_9evts_75570-0_resid.root");
  TFile *f_hough = TFile::Open("/sphenix/user/mitrankova/Hough/macro/hough_tracks.root");

  if (!f_resid || f_resid->IsZombie() || !f_hough || f_hough->IsZombie()) {
    std::cerr << "Error: failed to open one of the files!" << std::endl;
    return;
  }

  // --- Get residual tree ---
  TTree *t_resid = (TTree*)f_resid->Get("residualtree");
  if (!t_resid) {
    std::cerr << "Error: residual_tree not found in output_resid.root" << std::endl;
    return;
  }

  std::vector<float> *clusgx = nullptr;
  std::vector<float> *clusgy = nullptr;
  t_resid->SetBranchAddress("clusgx", &clusgx);
  t_resid->SetBranchAddress("clusgy", &clusgy);

  // --- Collect all cluster points ---
  std::vector<double> vx, vy;

  Long64_t nentries = t_resid->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {
    t_resid->GetEntry(i);
    if (!clusgx || !clusgy) continue;
    size_t nclus = std::min(clusgx->size(), clusgy->size());
    for (size_t j = 0; j < nclus; ++j) {
      vx.push_back(clusgx->at(j));
      vy.push_back(clusgy->at(j));
    }
  }

  // --- Get Hough parameters ---
  TTree *t_phase = (TTree*)f_hough->Get("phase_tree");
  if (!t_phase) {
    std::cerr << "Error: phase_tree not found in hough_tracks.root" << std::endl;
    return;
  }

  float rho, theta;
  t_phase->SetBranchAddress("rho", &rho);
  t_phase->SetBranchAddress("theta", &theta);

  std::vector<std::pair<double,double>> hough_lines;
  Long64_t nlines = t_phase->GetEntries();
  for (Long64_t i = 0; i < nlines; ++i) {
    t_phase->GetEntry(i);
    // convert to radians if stored in degrees
    double th = theta;
    if (std::fabs(theta) > 2 * TMath::Pi()) th = theta * TMath::DegToRad();
    hough_lines.emplace_back(rho, th);
  }

  // --- Draw everything ---
  TCanvas *c1 = new TCanvas("c1", "Cluster Positions and Hough Lines", 900, 800);
  TGraph *g = new TGraph(vx.size(), vx.data(), vy.data());
  g->SetTitle("Cluster positions and Hough lines;clusgx;clusgy");
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.5);
  g->SetMarkerColor(kBlue);
  g->Draw("AP");

  double xmin = g->GetXaxis()->GetXmin();
  double xmax = g->GetXaxis()->GetXmax();
  if (xmin == xmax) { xmin = -200; xmax = 200; }
  double ymin = g->GetYaxis()->GetXmin();
  double ymax = g->GetYaxis()->GetXmax();
  if (ymin == ymax) { ymin = -200; ymax = 200; }

  std::vector<TLine*> line_objs;
  int color = 2;
  for (auto [r, th] : hough_lines) {
    double cosT = std::cos(th);
    double sinT = std::sin(th);
    if (std::fabs(sinT) < 1e-6) continue;

    double y1 = (r - xmin * cosT) / sinT;
    double y2 = (r - xmax * cosT) / sinT;

    TLine *line = new TLine(xmin, y1, xmax, y2);
    line->SetLineColor(color);
    line->SetLineWidth(2);
    line->Draw("same");
     line_objs.push_back(line);

    color++;
    if (color > 9) color = 2; // cycle colors
  }

  c1->Update();
  c1->SaveAs("Hough_Tracks_Plot.png");
    TFile *fout = new TFile("combined_output.root", "RECREATE");
  fout->cd();
  g->Write();  // save TGraph
  for (size_t i = 0; i < line_objs.size(); ++i) {
    TString name = Form("hough_line_%zu", i);
    line_objs[i]->Write(name);
  }
  c1->Write("canvas_with_lines");
  fout->Close();
}

/*
 * This macro shows a minimum working example of running the tracking
 * hit unpackers with some basic seeding algorithms to try to put together
 * tracks. There are some analysis modules run at the end which package
 * hits, clusters, and clusters on tracks into trees for analysis.
 */

#include <fun4all/Fun4AllUtils.h>
#include <G4_ActsGeom.C>
#include <G4_Global.C>
#include <G4_Magnet.C>
#include <G4_Mbd.C>
#include <GlobalVariables.C>
#include <QA.C>
#include <Trkr_Clustering.C>
#include <Trkr_LaserClustering.C>
#include <Trkr_Reco.C>
#include <Trkr_RecoInit.C>
#include <Trkr_TpcReadoutInit.C>

#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/recoConsts.h>

#include <cdbobjects/CDBTTree.h>

#include <tpccalib/PHTpcResiduals.h>

#include <mvtxrawhitqa/MvtxRawHitQA.h>
#include <inttrawhitqa/InttRawHitQA.h>
#include <trackingqa/InttClusterQA.h>
#include <trackingqa/MicromegasClusterQA.h>
#include <trackingqa/MvtxClusterQA.h>
#include <trackingqa/TpcClusterQA.h>
#include <tpcqa/TpcRawHitQA.h>
#include <trackingqa/SiliconSeedsQA.h>
#include <trackingqa/TpcSeedsQA.h>
#include <trackingqa/TpcSiliconQA.h>

#include <trackingdiagnostics/TrackResiduals.h>
#include <trackingdiagnostics/TrkrNtuplizer.h>
#include <trackingdiagnostics/KshortReconstruction.h>
#include <kfparticle_sphenix/KFParticle_sPHENIX.h>
#include </sphenix/user/mitrankova/Hough/HoughTrackFinder/HoughTrackFinder.h>

#include <stdio.h>

R__LOAD_LIBRARY(libkfparticle_sphenix.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libcdbobjects.so)
R__LOAD_LIBRARY(libmvtx.so)
R__LOAD_LIBRARY(libintt.so)
R__LOAD_LIBRARY(libtpc.so)
R__LOAD_LIBRARY(libmicromegas.so)
R__LOAD_LIBRARY(libTrackingDiagnostics.so)
R__LOAD_LIBRARY(libtrackingqa.so)
R__LOAD_LIBRARY(libtpcqa.so)

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>
R__LOAD_LIBRARY(/sphenix/user/mitrankova/Hough/HoughTrackFinder/install/lib/libHoughTrackFinder.so)

class SkipFirstN : public SubsysReco {
 public:
  explicit SkipFirstN(int n) : SubsysReco("SkipFirstN"), target_(n) {}
  int process_event(PHCompositeNode*) override {
    if (count_ < target_) { ++count_; return Fun4AllReturnCodes::ABORTEVENT; }
    return Fun4AllReturnCodes::EVENT_OK;
  }
 private:
  int target_ = 0;
  int count_  = 0;
};




void Fun4All_raw_hit_TPC(
    const int nEvents = 9,
    const std::string filelist = "/sphenix/user/mitrankova/Hough/filelists/filelist_00075570_00000.list",
    const std::string outfilename = "TPC_Au_Au_ZeroField_1mrad_1evt_skip8_aligned",
    const bool convertSeeds = true,
    const int nSkip = 8)
{
  auto se = Fun4AllServer::instance();
  se->Verbosity(2);

  se->registerSubsystem(new SkipFirstN(nSkip));

  auto rc = recoConsts::instance();

   //input manager for QM production raw hit DST file
  std::ifstream ifs(filelist);
  std::string filepath;

  int i = 0;
  int runnumber = std::numeric_limits<int>::quiet_NaN();
  int segment = std::numeric_limits<int>::quiet_NaN();
  bool process_endpoints = false;
  
  while(std::getline(ifs,filepath))
    {
     std::cout << "Adding DST with filepath: " << filepath << std::endl; 
     if(i==0)
	   {
	    std::pair<int, int>
	      runseg = Fun4AllUtils::GetRunSegment(filepath);
	    runnumber = runseg.first;
	    segment = runseg.second;
	    rc->set_IntFlag("RUNNUMBER", runnumber);
	    rc->set_uint64Flag("TIMESTAMP", runnumber);
        
	   }
     if(filepath.find("ebdc") != std::string::npos)
     {
	    if(filepath.find("_0_") != std::string::npos or
	     filepath.find("_1_") != std::string::npos)
	    {
	      process_endpoints = true;
	    }
      }
      std::string inputname = "InputManager" + std::to_string(i);
      auto hitsin = new Fun4AllDstInputManager(inputname);
      hitsin->fileopen(filepath);
      se->registerInputManager(hitsin);
      i++;
    }
  
  rc->set_IntFlag("RUNNUMBER", runnumber);
  rc->set_IntFlag("RUNSEGMENT", segment);

  Enable::QA = false;
  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");//2024p022 newcdbtag 2024p018
  rc->set_uint64Flag("TIMESTAMP", runnumber);



  G4TRACKING::convert_seeds_to_svtxtracks = convertSeeds;
  std::cout << "Converting to seeds : " << G4TRACKING::convert_seeds_to_svtxtracks << std::endl;
  //TRACKING::reco_t0=-8.5;

  std::cout<< " run: " << runnumber
	   << " samples: " << TRACKING::reco_tpc_maxtime_sample
	   << " pre: " << TRACKING::reco_tpc_time_presample
	   << " vdrift: " << G4TPC::tpc_drift_velocity_reco
	   << std::endl;

 TRACKING::pp_mode = false;

  // distortion calibration mode
  /*
   * set to true to enable residuals in the TPC with
   * TPC clusters not participating to the ACTS track fit
   */
 // std::string outdir  = "/sphenix/tg/tg01/hf/mitrankova/AuAu_2025_ZeroField/75570_TPC/";
  std::string outdir  = "/sphenix/user/mitrankova/Hough/output_resid/";
  TString outfile = outdir + outfilename + "_" + to_string(nEvents) + "evts_" + runnumber + "-" + segment;
  std::string theOutfile = outfile.Data();

  FlagHandler *flag = new FlagHandler();
  se->registerSubsystem(flag);

  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");

  Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  TpcReadoutInit( runnumber );
  G4TPC::REJECT_LASER_EVENTS=true;
  G4TPC::ENABLE_MODULE_EDGE_CORRECTIONS = false;
  //Flag for running the tpc hit unpacker with zero suppression on
  TRACKING::tpc_zero_supp = true;

  //MVTX
  Enable::MVTX_APPLYMISALIGNMENT = true;
  ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT;

  //to turn on the default static corrections, enable the two lines below
  G4TPC::ENABLE_STATIC_CORRECTIONS = false;
  G4TPC::USE_PHI_AS_RAD_STATIC_CORRECTIONS = false;

  //to turn on the average corrections derived from simulation, enable the three lines below
  //note: these are designed to be used only if static corrections are also applied
  //G4TPC::ENABLE_AVERAGE_CORRECTIONS = true;
  //G4TPC::USE_PHI_AS_RAD_AVERAGE_CORRECTIONS = false;
  //G4TPC::average_correction_filename = std::string(getenv("CALIBRATIONROOT")) + "/distortion_maps/average_minus_static_distortion_inverted_10-new.root";


  G4MAGNET::magfield = "0.01";
  G4MAGNET::magfield_tracking = G4MAGNET::magfield;
  G4MAGNET::magfield_rescale = 1;
  
  
  TrackingInit();

  for(int felix=0; felix < 6; felix++)
  {
    Mvtx_HitUnpacking(std::to_string(felix));
  }
  for(int server = 0; server < 8; server++)
  {
    Intt_HitUnpacking(std::to_string(server));
  }
  ostringstream ebdcname;
  for(int ebdc = 0; ebdc < 24; ebdc++)
    {
      if(!process_endpoints)
	{
	  ebdcname.str("");
	  if(ebdc < 10)
	    {
	      ebdcname<<"0";
	    }
	  ebdcname<<ebdc;
	  Tpc_HitUnpacking(ebdcname.str());
	}
      
      else if(process_endpoints)
	{
	  for(int endpoint = 0; endpoint <2; endpoint++)
	    {
	      ebdcname.str("");
	      if(ebdc < 10)
		{
		  ebdcname<<"0";
		}
	      ebdcname<<ebdc <<"_"<<endpoint;
	      Tpc_HitUnpacking(ebdcname.str());
	    }
	}
    }

  Micromegas_HitUnpacking();

  Mvtx_Clustering();

  Intt_Clustering();

  Tpc_LaserEventIdentifying();

  auto tpcclusterizer = new TpcClusterizer;
  tpcclusterizer->Verbosity(0);
  tpcclusterizer->set_do_hit_association(G4TPC::DO_HIT_ASSOCIATION);
  tpcclusterizer->set_rawdata_reco();
  tpcclusterizer->set_reject_event(G4TPC::REJECT_LASER_EVENTS);
  se->registerSubsystem(tpcclusterizer);


  Micromegas_Clustering();

  Reject_Laser_Events();
  /*
   * Begin Track Seeding
   */


  auto silicon_Seeding = new PHActsSiliconSeeding;
  silicon_Seeding->Verbosity(0);
  silicon_Seeding->setStrobeRange(-5,5);
  // these get us to about 83% INTT > 1
  silicon_Seeding->setinttRPhiSearchWindow(0.4);
  silicon_Seeding->setinttZSearchWindow(2.0);
  silicon_Seeding->seedAnalysis(false);
  se->registerSubsystem(silicon_Seeding);

  auto merger = new PHSiliconSeedMerger;
  merger->Verbosity(0);
  se->registerSubsystem(merger);

  /*
   * Tpc Seeding
   */
  auto seeder = new PHCASeeding("PHCASeeding");
  double fieldstrength = std::numeric_limits<double>::quiet_NaN();  // set by isConstantField if constant
  bool ConstField = isConstantField(G4MAGNET::magfield_tracking, fieldstrength);
  if (ConstField)
  {
    seeder->useConstBField(true);
    seeder->constBField(fieldstrength);
  }
  else
  {
    seeder->set_field_dir(-1 * G4MAGNET::magfield_rescale);
    seeder->useConstBField(false);
    seeder->magFieldFile(G4MAGNET::magfield_tracking);  // to get charge sign right
  }
  seeder->Verbosity(0);
  seeder->SetLayerRange(7, 55);
  seeder->SetSearchWindow(2.,0.05); // z-width and phi-width, default in macro at 1.5 and 0.05
  seeder->SetClusAdd_delta_window(3.0,0.06); //  (0.5, 0.005) are default; sdzdr_cutoff, d2/dr2(phi)_cutoff
  //seeder->SetNClustersPerSeedRange(4,60); // default is 6, 6
  seeder->SetMinHitsPerCluster(0);
  seeder->SetMinClustersPerTrack(3);
  seeder->useFixedClusterError(true);
  seeder->set_pp_mode(true);
  se->registerSubsystem(seeder);

  // expand stubs in the TPC using simple kalman filter
  auto cprop = new PHSimpleKFProp("PHSimpleKFProp");
  cprop->set_field_dir(G4MAGNET::magfield_rescale);
  if (ConstField)
  {
    cprop->useConstBField(true);
    cprop->setConstBField(fieldstrength);
  }
  else
  {
    cprop->magFieldFile(G4MAGNET::magfield_tracking);
    cprop->set_field_dir(-1 * G4MAGNET::magfield_rescale);
  }
  cprop->useFixedClusterError(true);
  cprop->set_max_window(5.);
  cprop->Verbosity(0);
  cprop->set_pp_mode(true);
  se->registerSubsystem(cprop);

  // Always apply preliminary distortion corrections to TPC clusters before silicon matching
  // and refit the trackseeds. Replace KFProp fits with the new fit parameters in the TPC seeds.
  auto prelim_distcorr = new PrelimDistortionCorrection;
  prelim_distcorr->set_pp_mode(true);
  prelim_distcorr->Verbosity(0);
  se->registerSubsystem(prelim_distcorr);

  /*
   * Track Matching between silicon and TPC
   */
  // The normal silicon association methods
  // Match the TPC track stubs from the CA seeder to silicon track stubs from PHSiliconTruthTrackSeeding
  auto silicon_match = new PHSiliconTpcTrackMatching;
  silicon_match->Verbosity(0);
  silicon_match->set_pp_mode(TRACKING::pp_mode);
  if(G4TPC::ENABLE_AVERAGE_CORRECTIONS)
  {
    // for general tracking
    // Eta/Phi window is determined by 3 sigma window
    // X/Y/Z window is determined by 4 sigma window
    silicon_match->window_deta.set_posQoverpT_maxabs({-0.014,0.0331,0.48});
    silicon_match->window_deta.set_negQoverpT_maxabs({-0.006,0.0235,0.52});
    silicon_match->set_deltaeta_min(0.03);
    silicon_match->window_dphi.set_QoverpT_range({-0.15,0,0}, {0.15,0,0});
    silicon_match->window_dx.set_QoverpT_maxabs({3.0,0,0});
    silicon_match->window_dy.set_QoverpT_maxabs({3.0,0,0});
    silicon_match->window_dz.set_posQoverpT_maxabs({1.138,0.3919,0.84});
    silicon_match->window_dz.set_negQoverpT_maxabs({0.719,0.6485,0.65});
    silicon_match->set_crossing_deltaz_max(30);
    silicon_match->set_crossing_deltaz_min(2);

    // for distortion correction using SI-TPOT fit and track pT>0.5
    if (G4TRACKING::SC_CALIBMODE)
    {
      silicon_match->window_deta.set_posQoverpT_maxabs({0.016,0.0060,1.13});
      silicon_match->window_deta.set_negQoverpT_maxabs({0.022,0.0022,1.44});
      silicon_match->set_deltaeta_min(0.03);
      silicon_match->window_dphi.set_QoverpT_range({-0.15,0,0}, {0.09,0,0});
      silicon_match->window_dx.set_QoverpT_maxabs({2.0,0,0});
      silicon_match->window_dy.set_QoverpT_maxabs({1.5,0,0});
      silicon_match->window_dz.set_posQoverpT_maxabs({1.213,0.0211,2.09});
      silicon_match->window_dz.set_negQoverpT_maxabs({1.307,0.0001,4.52});
      silicon_match->set_crossing_deltaz_min(1.2);
    }
  }
  silicon_match->zeroField(true);
 // se->registerSubsystem(silicon_match);

  // Match TPC track stubs from CA seeder to clusters in the micromegas layers
  auto mm_match = new PHMicromegasTpcTrackMatching;
  mm_match->Verbosity(0);
  mm_match->set_rphi_search_window_lyr1(3.);
  mm_match->set_rphi_search_window_lyr2(15.0);
  mm_match->set_z_search_window_lyr1(30.0);
  mm_match->set_z_search_window_lyr2(3.);

  mm_match->set_min_tpc_layer(38);             // layer in TPC to start projection fit
  mm_match->set_test_windows_printout(false);  // used for tuning search windows only
  mm_match->zeroField(true);
  //se->registerSubsystem(mm_match);

  /*
   * End Track Seeding
   */

  /*
   * Either converts seeds to tracks with a straight line/helix fit
   * or run the full Acts track kalman filter fit
   */
  if (G4TRACKING::convert_seeds_to_svtxtracks)
  {
    auto converter = new TrackSeedTrackMapConverter;
    // Default set to full SvtxTrackSeeds. Can be set to
    // SiliconTrackSeedContainer or TpcTrackSeedContainer
    converter->setTrackSeedName("TpcTrackSeedContainer");
    converter->setFieldMap(G4MAGNET::magfield_tracking);
    converter->Verbosity(0);
    se->registerSubsystem(converter);
  }
  else
  {
    auto deltazcorr = new PHTpcDeltaZCorrection;
    deltazcorr->Verbosity(0);
    se->registerSubsystem(deltazcorr);

    // perform final track fit with ACTS
    auto actsFit = new PHActsTrkFitter;
    actsFit->Verbosity(0);
    actsFit->commissioning(G4TRACKING::use_alignment);
    // in calibration mode, fit only Silicons and Micromegas hits
    actsFit->fitSiliconMMs(G4TRACKING::SC_CALIBMODE);
    actsFit->setUseMicromegas(G4TRACKING::SC_USE_MICROMEGAS);
    actsFit->set_pp_mode(TRACKING::pp_mode);
    actsFit->set_use_clustermover(true);  // default is true for now
    actsFit->useActsEvaluator(false);
    actsFit->useOutlierFinder(false);
    actsFit->setFieldMap(G4MAGNET::magfield_tracking);
    actsFit->setDirectNavigation(true);
    se->registerSubsystem(actsFit);

    auto cleaner = new PHTrackCleaner();
    cleaner->Verbosity(0);
    cleaner->set_pp_mode(TRACKING::pp_mode);
    se->registerSubsystem(cleaner);

    if (G4TRACKING::SC_CALIBMODE)
    {
      /*
      * in calibration mode, calculate residuals between TPC and fitted tracks,
      * store in dedicated structure for distortion correction
      */
      auto residuals = new PHTpcResiduals;
      const TString tpc_residoutfile = theOutfile + "_PhTpcResiduals.root";
      residuals->setOutputfile(tpc_residoutfile.Data());
      residuals->setUseMicromegas(G4TRACKING::SC_USE_MICROMEGAS);

      // matches Tony's analysis
      residuals->setMinPt( 0.2 );

      // reconstructed distortion grid size (phi, r, z)
      residuals->setGridDimensions(36, 48, 80);
      se->registerSubsystem(residuals);
    }

  }

  auto finder = new PHSimpleVertexFinder;
  finder->Verbosity(0);
  finder->zeroField(true);
  finder->setDcaCut(0.5);
  finder->setTrackPtCut(0.3);
  finder->setBeamLineCut(1);
  finder->setTrackQualityCut(1000);
  finder->setNmvtxRequired(3);
  finder->setOutlierPairCut(0.1);
  se->registerSubsystem(finder);

  if (!G4TRACKING::convert_seeds_to_svtxtracks)
  {
    // Propagate track positions to the vertex position
    auto vtxProp = new PHActsVertexPropagator;
    vtxProp->Verbosity(0);
    vtxProp->fieldMap(G4MAGNET::magfield_tracking);
    se->registerSubsystem(vtxProp);

  }


      HoughTrackFinder *ana = new HoughTrackFinder("HoughTrackFinder");
    se->registerSubsystem(ana);

   TString residoutfile = theOutfile + "_resid.root";
   std::string residstring(residoutfile.Data());
   std::cout << "!!!!!!!!!!Residual output file: " << residstring << std::endl;

 auto resid = new TrackResiduals("TrackResiduals");
   resid->outfileName(residstring);
   resid->alignment(false);
   resid->vertexTree();
//   // adjust track map name
//   if(G4TRACKING::SC_CALIBMODE && !G4TRACKING::convert_seeds_to_svtxtracks)
//   {
//     resid->trackmapName("SvtxSiliconMMTrackMap");
//     if( G4TRACKING::SC_USE_MICROMEGAS )
//     { resid->set_doMicromegasOnly(true); }
//   }
  resid->zeroField();
   resid->clusterTree();
   resid->hitTree();
   resid->convertSeeds(G4TRACKING::convert_seeds_to_svtxtracks);

   resid->Verbosity(0);
   se->registerSubsystem(resid);

   
   
  //auto ntuplizer = new TrkrNtuplizer("TrkrNtuplizer");
  //se->registerSubsystem(ntuplizer);
   
  /*
    // To write an output DST
    TString dstfile = theOutfile;
   std::string theDSTFile = dstfile.Data();
   Fun4AllOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", theDSTFile.c_str());
   out->AddNode("Sync");
   out->AddNode("EventHeader");
   out->AddNode("TRKR_CLUSTER");
   out->AddNode("SiliconTrackSeedContainer");
   out->AddNode("TpcTrackSeedContainer");
   out->AddNode("SvtxTrackSeedContainer");
   out->AddNode("SvtxTrackMap");
   out->AddNode("SvtxVertexMap");
   out->AddNode("MbdVertexMap");
   out->AddNode("GL1RAWHIT");
   se->registerOutputManager(out);

  */
  

  



 // se->skip(nSkip);
  se->run(nEvents);
  se->End();
  se->PrintTimer();


  if (Enable::QA)
  {
    TString qaname = theOutfile + "_qa.root";
    std::string qaOutputFileName(qaname.Data());
    QAHistManagerDef::saveQARootFile(qaOutputFileName);
  }
  CDBInterface::instance()->Print();
  delete se;
  std::cout << "Finished" << std::endl;
  gSystem->Exit(0);
}
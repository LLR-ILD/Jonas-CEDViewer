/**
 *    @author: Szymon Daraszewicz, DESY summerstudent (original).
 *    @author:  iLCSoft/CEDViewer author list (DSTViewer).
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
 **/
// -- C++ STL headers.

// -- ROOT headers.
#include "TError.h"

// -- LCIO headers.
////#include "EVENT/MCParticle.h"
////#include "EVENT/ReconstructedParticle.h"

// -- Marlin headers.
#include "marlin/Exceptions.h"
////#include "ClusterShapes.h"
////#include "HelixClass.h"

// -- Includes for detector drawing.
////#include "DD4hep/DD4hepUnits.h"
////#include "DD4hep/Detector.h"
#include "DDMarlinCED.h"

// -- Header for this processor and other project-specific headers.
#include "ColorMap.h"
#include "DDDSTViewer.h"
#include "viewer_util.h"

// -- Using-declarations and global constants.
// Only in .cc files, never in .h header files!
const char kScale              = viewer_util::kLog;
const unsigned int kColorSteps = 256;
const double kEnMin            = 0.;
const double kEnMax            = 100.;
// Thickness of drawn lines.
const int kHelixSize           = 1;
const int kLineSize            = 1;
const int kHitSize             = 2;
const int kHitMarker           = 1;
// Scaling factor for the momenta at IP line length.
const double kMomScale         = 100;

enum ColorMaps {
  kHotColorMap = 1,
  kColdColorMap,
  kJetColorMap,
  kCyclicColorMap,
  kGreyColorMap,
  kBlueColorMap,
};

enum ConeColor {
  kConeWhite     = 0x99ffff,
  kConeBlack     = 0x000000,
  kConeRed       = 0xdd0000,
  kConeGreen     = 0x00dd00,
  kConeYellow    = 0xdddd00,
  kConeFuchsia   = 0xdd00dd,
  kConeBlue      = 0x0000dd,
  kConeDarkBlue  = 0x000066,
  kConeOrange    = 0xff9900,
  kConePink      = 0xee82ee,
  kConeViolet    = 0x660099,
  kConePurple    = 0x800080,
  kConeSilver    = 0xc0c0c0,
  kConeGray      = 0x808080,
  kConeSkyBlue   = 0x87ceeb,
  kConeKhaki     = 0xf0e68c,
};

enum Layers {
  // PFO momentum at interaction point.
  kMomentum         =  1,  // Key: 1.
  kMomBelowECut     = 11,  // Key: !.
  // Jet IP layer.
  kIPJet2           =  2,  // Key: 2.
  kIPJet3           =  3,  // Key: 3.
  kIPJet4           =  4,  // Key: 4.
  kIPJet5           =  5,  // Key: 5.
  kIPJet6           =  6,  // Key: 6.
  // Jet layers.
  kJetDefault       =  7,  // Key: 7.
  kJet2             =  7,  // Key: 7.
  kJet3             =  7,  // Key: 7.
  kJet4             =  7,  // Key: 7.
  kJet5             =  7,  // Key: 7.
  kJet6             =  7,  // Key: 7.
  // PFO layers.
  kTPC              =  8,  // Key: 8.
  kCluster          =  9,  // Key: 9.
  kHit              =  0,  // Key: 0.
  kTPCBelowECut     = 19,  // Key: (.
  kClusterBelowECut = 10,  // Key: ).
  kHitBelowECut     = 13,  // Key: #.
  // MC Particle Layers.
  kMCHiggs          = 20,  // Key: t.
  kMCHiggsDecay     = 22,  // Key: u.
  kMCTau            = 23,  // Key: i.
  kMCOther          = 24,  // Key: o.
};
// ----------------------------------------------------------------------------

DDDSTViewer aDDDSTViewer;

void DDDSTViewer::writeLayerDescription(void){
  DDMarlinCED::add_layer_description("Mom at IP",        kMomentum);
  DDMarlinCED::add_layer_description("Mom below E",      kMomBelowECut);
  DDMarlinCED::add_layer_description("2 IP Jets",        kIPJet2);
  DDMarlinCED::add_layer_description("3 IP Jets",        kIPJet3);
  DDMarlinCED::add_layer_description("4 IP Jets",        kIPJet4);
  DDMarlinCED::add_layer_description("5 IP Jets",        kIPJet5);
  DDMarlinCED::add_layer_description("6 IP Jets",        kIPJet6);
  DDMarlinCED::add_layer_description("Default Jet",      kJetDefault);
  DDMarlinCED::add_layer_description("2 Jets",           kJet2);
  DDMarlinCED::add_layer_description("3 Jets",           kJet3);
  DDMarlinCED::add_layer_description("4 Jets",           kJet4);
  DDMarlinCED::add_layer_description("5 Jets",           kJet5);
  DDMarlinCED::add_layer_description("6 Jets",           kJet6);
  DDMarlinCED::add_layer_description("TPC",              kTPC);
  DDMarlinCED::add_layer_description("Clusters",         kCluster);
  DDMarlinCED::add_layer_description("Hits",             kHit);
  DDMarlinCED::add_layer_description("TPC below E",      kTPCBelowECut);
  DDMarlinCED::add_layer_description("Clusters below E", kClusterBelowECut);
  DDMarlinCED::add_layer_description("Hits below E",     kHitBelowECut);
  DDMarlinCED::add_layer_description("Higgs",            kMCHiggs);
  DDMarlinCED::add_layer_description("XX, ZH->ZXX",      kMCHiggsDecay);
  DDMarlinCED::add_layer_description("Tau, ZH->tautauH", kMCTau);
  DDMarlinCED::add_layer_description("Other primary MC", kMCOther);
}

DDDSTViewer::DDDSTViewer() : Processor("DDDSTViewer") {
  _description = "CED based event display for DST files with dd4hep geometry.";

  registerInputCollection(
    LCIO::RECONSTRUCTEDPARTICLE,
    "ParticleCollection",
    "Name of the particle collection.",
    rp_col_name_,
    std::string("RecoParticles"));

  registerInputCollection(
    LCIO::MCPARTICLE,
   "MCCollection" ,
   "Name of the Monte Carlo particle collection (if it should be drawn)" ,
   mc_col_name_ ,
   std::string("MCParticlesSkimmed"));

  StringVec  default_jet_col_names{
    // There are no jet collections in the DST files.
    //"Durham_2Jets",
    //"Durham_3Jets"
    //"Durham_4Jets"
    //"Durham_5Jets"
    //"Durham_6Jets"
  };
  registerInputCollections(
    LCIO::RECONSTRUCTEDPARTICLE,
    "JetCollections" ,
    "Names of the jet collections." ,
    jet_col_names_,
    default_jet_col_names);

  registerProcessorParameter(
    "WaitForKeyboard",
    "Wait for Keyboard before proceed.",
    wait_for_keyboard_,
    true);

  registerProcessorParameter(
    "PtDrawCut",
    "Transverse mom. threshold to divide between low and high energy drawing.",
    pt_draw_cut_,
    (double) 0.0);

  IntVec default_h_decays{};  // The default choice is empty -> draw all events.
  registerProcessorParameter(
    "HiggsDecays",
    "If non-empty: List of MC Higgs decays for which the event is drawn."
    "Knows: 3 4 5 13 15 21 22 23 24 (c s b mu tau gluon photon Z W).",
    use_h_decays_,
    default_h_decays);

  // Add drawing related options.
  StringVec detailled_surfaces_example{};
  // detailled_surfaces_example.push_back("NameToBeDrawnDetailled");
  registerProcessorParameter(
    "DetailledDrawing" ,
    "List of detector names to be printed in more detail.",
    detailled_drawn_detector_surfaces_,
    detailled_surfaces_example,
    1);

  registerProcessorParameter(
    "DrawSurfaces",
    "Draw the geometry as a set of individual surfaces (if available) instead "
    "of simplified structures.",
    is_drawn_surfaces_,
    false);

  registerProcessorParameter(
    "HalfConeOpening",
    "Angle [rad] up to which the cones are drawn.",
    half_cone_opening_,
    0.15);
}

void DDDSTViewer::init() {
  DDMarlinCED::init(this);
  printParameters();
}

void DDDSTViewer::processRunHeader(LCRunHeader* run) {
  n_run_ = run->getRunNumber();
}

/**
 *  Mute annoying ROOT TColor warning.
 *
 * There is a ROOT TColor warning printed for each detector element in each
 * event, polluting the output. This is from the way the color of a detector
 * element is initialized in the DDMarlinCED class of the MarlinUtil package.
 * In the `getVisAttributes` method, the TColor initializer is (knowingly)
 * used on already existing color ids, to return their rgb values.
 * I assume this warning was added to TColor later than the time of writing of
 * DDMarlinCED. To avoid touching/personalizing the MarlinUtil package on top,
 * the following wrapper is utilized.
 * Quick and dirty, probably temporary fix.
 **/
void wrapperDrawDD4hepDetector(dd4hep::Detector& the_detector,
    bool is_drawn_surfaces_, StringVec detailled_drawn_detector_surfaces_) {
  Int_t old_root_error_ignore_level = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  DDMarlinCED::drawDD4hepDetector(the_detector, is_drawn_surfaces_,
      detailled_drawn_detector_surfaces_);
  gErrorIgnoreLevel = old_root_error_ignore_level;
}

void DDDSTViewer::processEvent(LCEvent* event) {
  n_event_ = event->getEventNumber();
  streamlog_out(DEBUG) << "  Processing event no " << n_event_
    << " - run " << n_run_ << std::endl << std::endl;
  if (skipUnwantedHiggsDecay(use_h_decays_, mc_col_name_, event)) {
    throw marlin::SkipEventException(this);
  }
  DDMarlinCED::newEvent(this, event);
  writeLayerDescription();
  // Get the detector instance (now dd4hep, replaced gear).
  dd4hep::Detector& the_detector = dd4hep::Detector::getInstance();
  // For now, we opt to draw the geometry as simplified structures instead of
  // individual surfaces. The list of detector names to draw in more detail is
  // empty.
  wrapperDrawDD4hepDetector(the_detector, is_drawn_surfaces_,
      detailled_drawn_detector_surfaces_);
  DDCEDPickingHandler &p_handler = DDCEDPickingHandler::getInstance();
  p_handler.update(event);
  // Add a line to disentangle DD4hep output from output actually added by this
  // processor itself.
  streamlog_out(MESSAGE) << std::endl
    << "-------------------------------------------------" << std::endl
    << std::endl;
// Drawing the reconstructed particles.
  streamlog_out(DEBUG) << "  Drawing RPs from collection: "
      << rp_col_name_ << "." << std::endl ;
  EVENT::LCCollection* rp_col = nullptr;
  try {
    rp_col = event->getCollection(rp_col_name_);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "RP collection " << rp_col_name_
      << " is not available!" << std::endl;
  }
  if (rp_col) {
    // Draw the helix of charged tracks (other than muons) only in the TPC
    // within the tracker volume.
    double helix_min_r = 0.0;
    double helix_max_r = viewer_util::getTrackerExtent(the_detector)[0];
    double helix_max_z = viewer_util::getTrackerExtent(the_detector)[1];
    viewer_util::AnyParticle ev_sum;
    // Obtain the detector's B-field for the per-event helix calculation.
    double b_vector[3];
    the_detector.field().combinedMagnetic(dd4hep::Position(0,0,0), b_vector);
    double b_field_z = b_vector[2] / dd4hep::tesla;
    for (int i=0; i < rp_col->getNumberOfElements(); ++i) {
      EVENT::ReconstructedParticle* rp =
        static_cast<EVENT::ReconstructedParticle*>(rp_col->getElementAt(i));
      // Conveniently store the necessary particle information in a struct.
      viewer_util::AnyParticle ap(rp);
      ev_sum += ap;
      // Reference point of the particle. So far, only origin is implemented.
      std::vector<double> ref_pt = {0, 0, 0};
      // Define the layers the rp components are printed to (E-dependent).
      int mom_layer = kMomentum;
      int tpc_layer = kTPC;
      int hit_layer = kHit;
      int clu_layer = kCluster;
      if (ap.getPt() < pt_draw_cut_) {
        mom_layer = kMomBelowECut;
        tpc_layer = kTPCBelowECut;
        hit_layer = kHitBelowECut;
        clu_layer = kClusterBelowECut;
      }
    // Two kinds of objects should be drawn for each reconstructed particle:
    // The track (in the case of charged RPs) and the cluster.
    // Lets start with the track.
      viewer_util::DSTColor track_color = returnTrackColor(rp->getType());
      DDMarlinCED::drawHelix(b_field_z, ap.charge,
          ref_pt[0], ref_pt[1], ref_pt[2],
          ap.x, ap.y, ap.z, tpc_layer,
          kHelixSize, track_color.hexa(),
          helix_min_r, helix_max_r, helix_max_z, rp->id());
      // For the momentum lines, both length and color are momentum dependent.
      viewer_util::DSTColor mom_line_color = returnRGBClusterColor(
          ap.getP(), kEnMin, kEnMax, kScale, kColdColorMap);
      ced_line_ID(ref_pt[0], ref_pt[1], ref_pt[2],
          kMomScale*ap.x, kMomScale*ap.y, kMomScale*ap.z, mom_layer,
          kLineSize, mom_line_color.hexa(), rp->id());
    // Now the clusters.
      if (rp->getClusters().size() > 0) {
        // Only draw a single cluster per particle.
        EVENT::Cluster* cluster = rp->getClusters()[0];
        // As input to different functions, both float and double array are
        // needed.
        //const float* const_ctr = cluster->getPosition();
        float const_ctr[3] = {cluster->getPosition()[0],
            cluster->getPosition()[1], cluster->getPosition()[2]};
        float center_f[3] = {const_ctr[0], const_ctr[1], const_ctr[2]};
        double center_d[3] = {const_ctr[0], const_ctr[1], const_ctr[2]};
        float phi = cluster->getIPhi();
        float theta = cluster->getITheta();
        if (phi == 0 && theta == 0) {
          phi = atan2(center_d[1], center_d[0]);
          theta = atan(sqrt(center_d[0]*center_d[0] + center_d[1]*center_d[1])
                       / center_d[2]);
        }
        // for particles centered at the origin
        if (ref_pt != std::vector<double>{0, 0, 0}){
          streamlog_out(WARNING) << "This is not yet implemented. The rotation "
            "operation so far is based on the reference point of the track "
            "being the center of the detector." << std::endl;
          return;
        }
        double rotate[] = {0.0, theta*180/M_PI, phi*180/M_PI};
        float en_cluster = cluster->getEnergy();
        double cl_scale = returnClusterSize(en_cluster, kEnMin, kEnMax);
        double cl_size[] = {cl_scale, cl_scale, 4*cl_scale};
        viewer_util::DSTColor clu_color = returnRGBClusterColor(
            en_cluster, kEnMin, kEnMax, kScale, kJetColorMap);
        if (true) {   // Enables picking (and adds the small hit dot).
          ced_hit_ID(center_d[0], center_d[1], center_d[2], kHitMarker,
              hit_layer, (int)(sqrt(2)*cl_size[0]/4),
              clu_color.hexa(), cluster->id());
        }
        if (true) {  // 3 2D ellipses of the cluster extend.
          clu_color.setAlpha(0x88);
          ced_cluellipse_r_ID((float) cl_size[0], (float) cl_size[2],
              center_f, rotate, clu_layer, clu_color.hexa(),
              cluster->id());  // Originally: BACKUP_LAYER.
        }
        if (true) {  // A filled ellipse. Combined with the previous draw, this
        // creates a nice ellipsoid (3D).
          clu_color.setAlpha(0xCC);
          ced_ellipsoid_r_ID(cl_size, center_d, rotate, clu_layer,
              clu_color.hexa(), cluster->id());  // Originally: BACKUP_LAYER2.
        }
        ////if (false) {  // Sides of a barrel around the ellipsoid (visual aid for
        ////// the cluster direction). I don't think it's nice.
        ////  int cylinder_sides = 30;
        ////  ced_geocylinder_r(cl_size[0]/2, cl_size[2], center_d, rotate,
        ////      cylinder_sides, clu_color.hexa(), clu_layer);
        ////}
        ////if (false) {  // A line with arrow. For what, I don't know.
        //// // The line.
        ////  float clu_dir[3] = {
        ////    std::sin(theta)*std::cos(phi),
        ////    std::sin(theta)*std::sin(phi),
        ////    std::cos(theta),
        ////  };
        ////  int    l_width = 2 * kLineSize;
        ////  double l_length = 0.5 * cl_size[2];
        ////  double arrow_length = 0.2 * l_length;
        ////  float x_end = center_f[0] + (l_length - arrow_length) * clu_dir[0];
        ////  float y_end = center_f[1] + (l_length - arrow_length) * clu_dir[1];
        ////  float z_end = center_f[2] + (l_length - arrow_length) * clu_dir[2];
        ////  ced_line_ID(center_f[0], center_f[1], center_f[2],
        ////      x_end, y_end, z_end, clu_layer,  // Originally: CLU_LAYER.
        ////      l_width, clu_color.hexa(), cluster->id());
        //// // The direction arrow.
        ////  float x_arrow = center_f[0] + (l_length) * clu_dir[0];
        ////  float y_arrow = center_f[1] + (l_length) * clu_dir[1];
        ////  float z_arrow = center_f[2] + (l_length) * clu_dir[2];
        ////  ced_line_ID(x_end, y_end, z_end,
        ////      x_arrow, y_arrow, z_arrow, clu_layer,  // Originally: CLU_LAYER.
        ////      l_width, clu_color.hexa(), cluster->id());
        ////}
      }
    }
    unsigned int n_ticks = 6;
    showLegendSpectrum(kScale, kJetColorMap, kEnMin, kEnMax, n_ticks);
    streamlog_out(MESSAGE) << std::endl
      << "Total Energy and Momentum Balance of Event" << std::endl
      << "  Energy = " << ev_sum.E
      << "  Charge = " << ev_sum.charge << std::endl
      << "  PX     = " << ev_sum.x
      << "  PY     = " << ev_sum.y
      << "  PZ     = " << ev_sum.z << std::endl << std::endl;
    streamlog_out(DEBUG) << "  Setup properties: " << std::endl
      << "    B-Field = " << b_field_z << " T." << std::endl << std::endl;
  }

// Now we draw the jets (if any ).
  for (size_t i_col = 0; i_col < jet_col_names_.size(); ++i_col) {
    streamlog_out(DEBUG) << "  Drawing jets from collection: "
      << jet_col_names_[i_col] << "." << std::endl ;
    EVENT::LCCollection* jet_col = nullptr;
    try {
      jet_col = event->getCollection(jet_col_names_[i_col]);
    } catch (DataNotAvailableException &e) {
      streamlog_out(ERROR) << "Jet collection " << jet_col_names_[i_col]
        << " is not available!" << std::endl;
    }
    if (jet_col) {
      int jet_layer = returnJetLayer(jet_col_names_[i_col]);
      for (int i_elem = 0; i_elem < jet_col->getNumberOfElements(); ++i_elem) {
        EVENT::ReconstructedParticle* jet =
            static_cast<EVENT::ReconstructedParticle*>(
                  jet_col->getElementAt(i_elem));
        viewer_util::AnyParticle ap(jet);
        double ref_pt[3] = {0., 0., 0. };
        double rotate[3] = {0., ap.getTheta()*180./M_PI, ap.getPhi()*180./M_PI};
        const ReconstructedParticleVec & rp_vec = jet->getParticles();
        double jet_pt = 0.0;
        for (unsigned int i_part = 0; i_part < rp_vec.size(); ++i_part) {
          viewer_util::AnyParticle pp(rp_vec[i_part]);
          double mom_transvers_to_jet = pp.getP()
              - (pp.x*ap.x + pp.y*ap.y + pp.z*ap.z) / sqrt(ap.getP());
          jet_pt += mom_transvers_to_jet;
          viewer_util::DSTColor part_color = returnJetColor(i_elem);
          int layer_ip = returnIpLayer(jet_col_names_[i_col]);;
          ced_line_ID(ref_pt[0], ref_pt[1], ref_pt[2],
              kMomScale*pp.x, kMomScale*pp.y, kMomScale*pp.z,
              layer_ip, kLineSize, part_color.hexa(), rp_vec[i_part]->id());
        }
        viewer_util::DSTColor jet_color = returnJetColor(i_elem);
        int tpc_hit = 104;
        for (int k = 1; k < tpc_hit; ++k) {
          ced_line_ID((k-1)/4*ap.x, (k-1)/4*ap.y, (k-1)/4*ap.z,
            k/4*ap.x, k/4*ap.y, k/4*ap.z,
            jet_layer, kLineSize, jet_color.hexa(), jet->id());
        }
        ced_hit_ID((tpc_hit-1)/4*ap.x, (tpc_hit-1)/4*ap.y, (tpc_hit-1)/4*ap.z,
            kHitMarker, jet_layer, kHitSize, jet_color.hexa(), jet->id());
        // Draw a cone. Have the cone length depend on the particle momentum.
        float cone_base = 50 + 20 * jet_pt;
        float cone_height = 25 * ap.getP();
        float jet_color_rgba[4];
        jet_color.setAlpha(0x33);
        jet_color.update_rgba(jet_color_rgba);
        ced_cone_r_ID(cone_base, cone_height, ref_pt, rotate,
            jet_layer, jet_color_rgba, jet->id());
      }
    }
  }
// My MC Particle addition:
// Adapted from the Jet loop above.
  streamlog_out(DEBUG) << "  Drawing MC from collection: " << mc_col_name_
    << "." << std::endl;
  EVENT::LCCollection* mc_col = nullptr;
  try {
    mc_col = event->getCollection(mc_col_name_);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "MC collection " << mc_col_name_
      << " is not available!" << std::endl;
  }
  if (mc_col) {
    for (int j = 0; j < mc_col->getNumberOfElements(); ++j) {
      EVENT::MCParticle* mcp = static_cast<EVENT::MCParticle*>(
          mc_col->getElementAt(j));
      // Only display the some of the MC particles.
      // Also, define colors depending on the PDG.
      viewer_util::DSTColor mc_color = kConeBlue;
      int mc_layer = kMCOther;
      // Draw all primary MC particles (with special colors for Higgs and tau).
      if (mcp->getParents().size() == 0) {
        if(mcp->getPDG() == 25) {
          mc_color.setHexa(kConeGray);
          mc_layer = kMCHiggs;
        } else if(abs(mcp->getPDG()) == 15) {
          mc_color.setHexa(kConeGreen);
          mc_layer = kMCTau;
        // Tiny energy photons etc are not of interest.
        } else if(mcp->getEnergy() < 2) {
          continue;
        }
      // Draw the direct children of the Higgs.
      } else if (mcp->getParents()[0]->getPDG() == 25 && mcp->getPDG() != 25) {
        mc_color.setHexa(kConeViolet);
        mc_layer = kMCHiggsDecay;
        streamlog_out(MESSAGE) << "H -> " << mcp->getPDG() << "." << std::endl;
      // Ignore all remaining MC particles.
      } else {
        continue;
      }
      viewer_util::AnyParticle ap(mcp);
      double ref_pt[3] = {0., 0., 0. };
      double rotate[3] = {0., ap.getTheta()*180./M_PI, ap.getPhi()*180./M_PI};
      // Trajectory of the MC tracks is not necessary.
      // But we need a line for picking.
      int tpc_hit = 104;
      for (int k = 1; k < tpc_hit; ++k) {
        ced_line_ID((k-1)/4*ap.x, (k-1)/4*ap.y, (k-1)/4*ap.z,
          k/4*ap.x, k/4*ap.y, k/4*ap.z,
          mc_layer, kLineSize, mc_color.hexa(), mcp->id());
      }
      ced_hit_ID((tpc_hit-1)/4*ap.x, (tpc_hit-1)/4*ap.y, (tpc_hit-1)/4*ap.z,
          kHitMarker, mc_layer, kHitSize, mc_color.hexa(), mcp->id());
      // Draw a cone. Have the cone length depend on the particle momentum.
      ////float cone_base = 20 + 2 * ap.getP();
      // Longer cones were considered (to aid seeing which clusters are inside
      // the cone). At least with the current layer choices, a cone of any
      // transparency that reaches into the cluster region hides the clusters.
      // As this counteracts the idea of simplyfication, we stay with short
      // cones.
      float cone_height =  25 * ap.getP();
      float cone_base = 2 * tan(half_cone_opening_) * cone_height;
      float mc_color_rgba[4];
      mc_color.setAlpha(0x30);
      mc_color.update_rgba(mc_color_rgba);
      ced_cone_r_ID(cone_base, cone_height, ref_pt, rotate,
          mc_layer, mc_color_rgba, mcp->id());
    }
  }
  DDMarlinCED::draw(this, wait_for_keyboard_);
}

void DDDSTViewer::showLegendSpectrum(char scale, int color_map,
    double en_min, double en_max, unsigned int ticks) {
  // Legend colour matrix.
  unsigned int** rgb_matrix = new unsigned int*[kColorSteps];
  for (unsigned int i=0; i<kColorSteps; ++i){
    rgb_matrix[i] = new unsigned int[3];  // To avoid seg fault.
    ColorMap::selectColorMap(color_map)(rgb_matrix[i], i, 0.0, kColorSteps);
  }
  ced_legend(en_min, en_max, kColorSteps, rgb_matrix, ticks, scale);
}

viewer_util::DSTColor DDDSTViewer::returnRGBClusterColor(double energy,
    double cutoff_min, double cutoff_max, char scale, int color_map) {
  if (color_map < 0 || color_map > 6) {
    streamlog_out(ERROR) << "Wrong color_map parameter!" << std::endl;
  }
  // Implicit conversion of the double returned from the function to an int in
  // the range of kColorSteps.
  int color_delta = viewer_util::convertScales(energy, cutoff_max, cutoff_min,
    kColorSteps, 0, viewer_util::ScaleMapping(scale));
  unsigned int rgb[] = {0, 0, 0};
  ColorMap::selectColorMap(color_map)(rgb, color_delta, 0, kColorSteps);
  viewer_util::DSTColor color;
  color.setRGB(rgb[0], rgb[1], rgb[2]);
  return color;
}

double DDDSTViewer::returnClusterSize(double en_cluster,
     double cutoff_min, double cutoff_max) {
  // None-zero to ensure visibility of small clusters.
  int clu_size_min =  30; ////10;
  // Not related to the energy scale of the event, but the size of the ECAL.
  int clu_size_max = 300; ////120;
  return viewer_util::convertScales(en_cluster, cutoff_max, cutoff_min,
    clu_size_max, clu_size_min, viewer_util::ScaleMapping(kScale));
}

int DDDSTViewer::returnIpLayer(std::string jet_col_name) {
  int layer = kIPJet2; // The default value.
  if        (jet_col_name == "Durham_2Jets") {
    layer = kIPJet2;
  } else if (jet_col_name == "Durham_3Jets") {
    layer = kIPJet3;
  } else if (jet_col_name == "Durham_4Jets") {
    layer = kIPJet4;
  } else if (jet_col_name == "Durham_5Jets") {
    layer = kIPJet5;
  } else if (jet_col_name == "Durham_6Jets") {
    layer = kIPJet6;
  }
  return layer;
}

int DDDSTViewer::returnJetLayer(std::string jet_col_name) {
  int layer = kJetDefault; // The default value.
  if        (jet_col_name == "Durham_2Jets") {
    layer = kJet2;
  } else if (jet_col_name == "Durham_3Jets") {
    layer = kJet3;
  } else if (jet_col_name == "Durham_4Jets") {
    layer = kJet4;
  } else if (jet_col_name == "Durham_5Jets") {
    layer = kJet5;
  } else if (jet_col_name == "Durham_6Jets") {
    layer = kJet6;
  }
  return layer;
}

viewer_util::DSTColor DDDSTViewer::returnTrackColor(int particle_type) {
  switch (particle_type) {
    case   211:  // Pion+.
      return kConeRed;
    case  -211:  // Pion-.
      return kConeOrange;
    case    22:  // Photon.
      return kConeYellow;
    case    11:  // Electron.
      return kConeDarkBlue;
    case   -11:  // Positron.
      return kConeViolet;
    case    13:  // Muon.
      return kConeGreen;
    case   -13:  // Anti-muon.
      return kConeGreen;
    case  2112:  // Neutron.
    case -2112:  // Anti-neutron.
      return kConeWhite;
    default:
      streamlog_out(DEBUG) << "  Unconsidered particle type " << particle_type
      << ". A default color is chosen." << std::endl;
      return kConeBlack;
  }
}

viewer_util::DSTColor DDDSTViewer::returnJetColor(int col_number) {
  ////int transparency = 0x66;
  viewer_util::DSTColor color;
  switch (col_number) {
  case 0:
    color.setHexa(kConeYellow);
    break;
  case 1:
    color.setHexa(kConeRed);
    break;
  case 2:
    color.setHexa(kConeWhite);
    break;
  case 3:
    color.setHexa(kConeFuchsia);
    break;
  case 4:
    color.setHexa(kConeBlack);
    break;
  case 5:
    color.setHexa(kConeGreen);
    break;
  default:
    color.setHexa(kConeBlue);
    break;
  }
  ////color.addAlpha(-transparency);
  return color;
}

bool DDDSTViewer::skipUnwantedHiggsDecay(
    IntVec h_decays_to_draw, std::string mc_col_name, EVENT::LCEvent* event) {
  if (h_decays_to_draw.size() == 0) return false;
  std::cout << "ENTER";
  EVENT::LCCollection* mc_col = nullptr;
  try {
    mc_col = event->getCollection(mc_col_name);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "MC collection " << mc_col_name
      << " is not available!" << std::endl;
  }
  if (mc_col) {
    for (int j = 0; j < mc_col->getNumberOfElements(); ++j) {
      EVENT::MCParticle* mcp = static_cast<EVENT::MCParticle*>(
          mc_col->getElementAt(j));
      if (mcp->getParents().size() != 1) continue;
      if (mcp->getParents()[0]->getPDG() == 25 && mcp->getPDG() != 25) {
        for (auto pdg_to_use : h_decays_to_draw) {
          std::cout << pdg_to_use;
          if (mcp->getPDG() == pdg_to_use) return false;
        }
      }
    }
    // The only place for which true should be returned: The MC collection
    // exists, Higgs decay modes are specified, but the event does not have one
    // of the required decay modes.
    return true;
  }
  return false;
}
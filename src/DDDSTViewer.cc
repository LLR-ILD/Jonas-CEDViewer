/**
 *    @author: Szymon Daraszewicz, DESY summerstudent (original).
 *    @author:  iLCSoft/CEDViewer author list (DSTViewer).
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
 **/
// -- C++ STL headers.

// -- ROOT headers.

// -- LCIO headers.
////#include "EVENT/MCParticle.h"
////#include "EVENT/ReconstructedParticle.h"

// -- Marlin headers.
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
  kConeWhite     = 0x555555,
  kConeBlack     = 0x000000,
  kConeRed       = 0x550000,
  kConeGreen     = 0x008000,
  kConeYellow    = 0x555500,
  kConeFuchsia   = 0x550055,
  kConeBlue      = 0x000055,
  kConeOrange    = 0x55a500,
  kConeViolet    = 0xee82ee,
  kConePurple    = 0x800080,
  kConeSilver    = 0xc0c0c0,
  kConeGold      = 0x55d700,
  kConeGray      = 0x808080,
  kConeAqua      = 0x005555,
  kConeSkyBlue   = 0x87ceeb,
  kConeLightBlue = 0xa558e6,
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
    1);

  registerProcessorParameter(
    "EDrawCut",
    "Energy threshold to divide between low and high energy drawing.",
    e_draw_cut_,
    (double) 0.0);
}

void DDDSTViewer::init() {
  DDMarlinCED::init(this);
  printParameters();
}

void DDDSTViewer::processRunHeader(LCRunHeader* run) {
  n_run_ = run->getRunNumber();
}

void DDDSTViewer::processEvent(LCEvent* event) {
  n_event_ = event->getEventNumber();
  streamlog_out(DEBUG) << "Processing event no " << n_event_
    << " - run " << n_run_ << std::endl;
  DDMarlinCED::newEvent(this, event);
  writeLayerDescription();
  // Get the detector instance (now dd4hep, replaced gear).
  dd4hep::Detector& the_detector = dd4hep::Detector::getInstance();
  // For now, we opt to draw the geometry as simplified structures instead of
  // individual surfaces. The list of detector names to draw in more detail is
  // empty.
  DDMarlinCED::drawDD4hepDetector(the_detector, false, StringVec{});
  DDCEDPickingHandler &p_handler = DDCEDPickingHandler::getInstance();
  p_handler.update(event);
// Drawing the reconstructed particles.
  streamlog_out(DEBUG) << "Drawing RPs from collection: "
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
      if (ap.E < e_draw_cut_) {
        mom_layer = kMomBelowECut;
        tpc_layer = kTPCBelowECut;
        hit_layer = kHitBelowECut;
        clu_layer = kClusterBelowECut;
      }
    // Two kinds of objects should be drawn for each reconstructed particle:
    // The track (in the case of charged RPs) and the cluster.
    // Lets start with the track.
      int track_color = returnTrackColor(rp->getType());
      DDMarlinCED::drawHelix(b_field_z, ap.charge,
          ref_pt[0], ref_pt[1], ref_pt[2],
          ap.x, ap.y, ap.z, tpc_layer,
          kHelixSize, track_color,
          helix_min_r, helix_max_r, helix_max_z, rp->id());
      // For the momentum lines, both length and color are momentum dependent.
      int mom_line_color = returnRGBClusterColor(ap.getP(), kEnMin, kEnMax,
          kScale, kColdColorMap);
      ced_line_ID(ref_pt[0], ref_pt[1], ref_pt[2],
          kMomScale*ap.x, kMomScale*ap.y, kMomScale*ap.z, mom_layer,
          kLineSize, mom_line_color, rp->id());
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
          streamlog_out(DEBUG3)	<< "This is not yet implemented. The rotation "
            "operation so far is based on the reference point of the track "
            "being the center of the detector." << std::endl;
          return;
        }
        double rotate[] = {0.0, theta*180/M_PI, phi*180/M_PI};
        float en_cluster = cluster->getEnergy();
        double cl_scale = returnClusterSize(en_cluster, kEnMin, kEnMax);
        double cl_size[] = {cl_scale, cl_scale, 4*cl_scale};
        int color = returnRGBClusterColor(en_cluster, kEnMin, kEnMax, kScale,
            kJetColorMap);
        if (true) {   // Enables picking (and adds the small hit dot).
          ced_hit_ID(center_d[0], center_d[1], center_d[2], kHitMarker,
            hit_layer, (int)(sqrt(2)*cl_size[0]/4), color, cluster->id());
        }
        if (true) {  // 3 2D ellipses of the cluster extend.
          int rgba = addTransparencyToColor(color, 0x66);
          ced_cluellipse_r_ID((float) cl_size[0], (float) cl_size[2],
              center_f, rotate,
              clu_layer, rgba, cluster->id());  // Originally: BACKUP_LAYER.
        }
        if (true) {  // A filled ellipse. Combined with the previous draw, this
        // creates a nice ellipsoid (3D).
          int rgba = addTransparencyToColor(color, 0xCC);
          ced_ellipsoid_r_ID(cl_size, center_d, rotate,
              clu_layer, rgba, cluster->id());  // Originally: BACKUP_LAYER2.
        }
        ////if (false) {  // Sides of a barrel around the ellipsoid (visual aid for
        ////// the cluster direction). I don't think it's nice.
        ////  int cylinder_sides = 30;
        ////  ced_geocylinder_r(cl_size[0]/2, cl_size[2], center_d, rotate,
        ////      cylinder_sides, color, clu_layer);
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
        ////      l_width, color, cluster->id());
        //// // The direction arrow.
        ////  // Have the arrow colored a bit differently (greener). DOn't see why..
        ////  int color_arrow = color + 0x009900;
        ////  float x_arrow = center_f[0] + (l_length) * clu_dir[0];
        ////  float y_arrow = center_f[1] + (l_length) * clu_dir[1];
        ////  float z_arrow = center_f[2] + (l_length) * clu_dir[2];
        ////  ced_line_ID(x_end, y_end, z_end,
        ////      x_arrow, y_arrow, z_arrow, clu_layer,  // Originally: CLU_LAYER.
        ////      l_width, color_arrow, cluster->id());
        ////}
      }
    }
    unsigned int n_ticks = 6;
    showLegendSpectrum(kScale, kJetColorMap, kEnMin, kEnMax, n_ticks);
    streamlog_out(MESSAGE) << std::endl
      << "Total Energy and Momentum Balance of Event" << std::endl
      << "Energy  = " << ev_sum.E
      << " PX     = " << ev_sum.x
      << " PY     = " << ev_sum.y
      << " PZ     = " << ev_sum.z
      << " Charge = " << ev_sum.charge << std::endl << std::endl ;
    streamlog_out(DEBUG2) << "Setup properties: " << std::endl
      << "B-Field = " << b_field_z << " T." << std::endl ;
  }

// Now we draw the jets (if any ).
  for (size_t i_col = 0; i_col < jet_col_names_.size(); ++i_col) {
    streamlog_out(DEBUG) << "Drawing jets from collection: "
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
          int part_color = returnJetColor(i_elem);
          int layer_ip = returnIpLayer(jet_col_names_[i_col]);;
          ced_line_ID(ref_pt[0], ref_pt[1], ref_pt[2],
              kMomScale*pp.x, kMomScale*pp.y, kMomScale*pp.z,
              layer_ip, kLineSize, part_color, rp_vec[i_part]->id());
        }
        int jet_color = returnJetColor(i_elem);
        int tpc_hit = 104;
        for (int k = 1; k < tpc_hit; ++k) {
          ced_line_ID((k-1)/4*ap.x, (k-1)/4*ap.y, (k-1)/4*ap.z,
            k/4*ap.x, k/4*ap.y, k/4*ap.z,
            jet_layer, kLineSize, jet_color, jet->id());
        }
        ced_hit_ID((tpc_hit-1)/4*ap.x, (tpc_hit-1)/4*ap.y, (tpc_hit-1)/4*ap.z,
            kHitMarker, jet_layer, kHitSize, jet_color, jet->id());
        // Draw a cone. Have the cone length depend on the particle momentum.
        float cone_base = 50 + 20 * jet_pt;
        float cone_height = 25 * ap.getP();
        float rgba_jet_color[4];
        viewer_util::fromIntToRGBA(jet_color, rgba_jet_color);
        ced_cone_r_ID(cone_base, cone_height, ref_pt, rotate,
            jet_layer, rgba_jet_color, jet->id());
      }
    }
  }
// My MC Particle addition:
// Adapted from the Jet loop above.
  streamlog_out(DEBUG) << "Drawing MC from collection: " << mc_col_name_
    << "." << std::endl ;
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
      int mc_color = kConeGray;
      int mc_layer = kMCOther;
      // Draw all primary MC particles (with special colors for Higgs and tau).
      if (mcp->getParents().size() == 0) {
        if(mcp->getPDG() == 25) {
          mc_color = kConeGold;
          mc_layer = kMCHiggs;
        } else if(abs(mcp->getPDG()) == 15) {
          mc_color = kConeGreen;
          mc_layer = kMCTau;
        // Tiny energy photons etc are not of interest.
        } else if(mcp->getEnergy() < 2) {
          continue;
        }
      // Draw the direct children of the Higgs.
      } else if (mcp->getParents()[0]->getPDG() == 25 && mcp->getPDG() != 25) {
        mc_color = kConeViolet;
        mc_layer = kMCHiggsDecay;
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
          mc_layer, kLineSize, mc_color, mcp->id());
      }
      ced_hit_ID((tpc_hit-1)/4*ap.x, (tpc_hit-1)/4*ap.y, (tpc_hit-1)/4*ap.z,
          kHitMarker, mc_layer, kHitSize, mc_color, mcp->id());
      // Draw a cone. Have the cone length depend on the particle momentum.
      float cone_base = 20 + 2 * ap.getP();
      float cone_height =   25 * ap.getP();
      float rgba_jet_color[4];
      viewer_util::fromIntToRGBA(mc_color, rgba_jet_color);
      ced_cone_r_ID(cone_base, cone_height, ref_pt, rotate,
          mc_layer, rgba_jet_color, mcp->id());
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

int DDDSTViewer::returnRGBClusterColor(double energy,
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
  int hex_color = ColorMap::RGB2HEX(rgb[0], rgb[1], rgb[2]);
  return hex_color;
}

double DDDSTViewer::returnClusterSize(double en_cluster,
     double cutoff_min, double cutoff_max) {
  // None-zero to ensure visibility of small clusters.
  int clu_size_min = 10;
  // Not related to the energy scale of the event, but the size of the ECAL.
  int clu_size_max = 120;
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

int DDDSTViewer::addTransparencyToColor(int color, int transparency) {
  return (transparency<<24) + color;
}

enum TrackColor {
  kTrackBlack    = 0x999999,
  kTrackRed      = 0x990000,
  kTrackOrange   = 0x996600,
  kTrackYellow   = 0xffff00,
  kTrackGreen    = 0x00ff00,
  kTrackDarkBlue = 0x660066,
  kTrackViolet   = 0x660099,
  kTrackWhite    = 0x99FFFF,
};

int DDDSTViewer::returnTrackColor(int particle_type) {
  switch (particle_type) {
    case   211:  // Pion+.
      return kTrackRed;
    case  -211:  // Pion-.
      return kTrackOrange;
    case    22:  // Photon.
      return kTrackYellow;
    case    11:  // Electron.
      return kTrackDarkBlue;
    case   -11:  // Positron.
      return kTrackViolet;
    case    13:  // Muon.
      return kTrackGreen;
    case   -13:  // Anti-muon.
      return kTrackGreen;
    case  2112:  // Neutron.
      return kTrackWhite;
    case -2112:  // Anti-neutron.
      return kTrackWhite;
    default:
      streamlog_out(DEBUG) << "Unconsidered particle type " << particle_type
      << ". A default color is chosen." << std::endl;
      return kTrackBlack;
  }
}

int DDDSTViewer::returnJetColor(int col_number) {
  int transparency = 0x66;
  switch (col_number) {
  case 0:
    return addTransparencyToColor(kConeYellow, transparency);
  case 1:
    return addTransparencyToColor(kConeRed, transparency);
  case 2:
    return addTransparencyToColor(kConeWhite, transparency);
  case 3:
    return addTransparencyToColor(kConeFuchsia, transparency);
  case 4:
    return addTransparencyToColor(kConeBlack, transparency);
  case 5:
    return addTransparencyToColor(kConeGreen, transparency);
  default:
    return addTransparencyToColor(kConeBlue, transparency);
  }
}
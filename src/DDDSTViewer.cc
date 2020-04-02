/**
 *    @author: Szymon Daraszewicz, DESY summerstudent (original).
 *    @author:  iLCSoft/CEDViewer author list (DSTViewer).
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
 **/
// TODO: Clean up includes.
#include "DD4hep/Detector.h"
#include "DDMarlinCED.h"
#include "viewer_util.h"

#include "DDDSTViewer.h"
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/MCParticle.h>
#include <iostream>
#include "ClusterShapes.h"

#include "HelixClass.h"
#include <math.h>
#include "ColorMap.h"

enum layers {
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
    "Durham_2Jets",
    "Durham_3Jets"
    "Durham_4Jets"
    "Durham_5Jets"
    "Durham_6Jets"
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
    0.0f);
}

void DDDSTViewer::init() {
  DDMarlinCED::init(this);
  printParameters();
}

void DDDSTViewer::processRunHeader( LCRunHeader* run) {
  n_run_ = run->getRunNumber();
}

void DDDSTViewer::processEvent(LCEvent* event) {
  n_event_ = event->getEventNumber();
  streamlog_out(DEBUG) << "Processing event no " << n_event_
    << " - run " << n_run_ << std::endl;
  DDMarlinCED::newEvent(this);
  DDDSTViewer::writeLayerDescription();
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
      << rp_col_name_ << std::endl ;
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
    double helix_max_r = viewer_util::getTrackerExtent(the_detector)[0];
    double helix_max_z = viewer_util::getTrackerExtent(the_detector)[1];
    viewer_util::AnyParticle ev_sum;
    // Obtain the detector's B-field for the per-event helix calculation.
    double* b_vector = new double[3];
    the_detector.field().combinedMagnetic(dd4hep::Position(0,0,0), b_vector);
    double b_field_z = b_vector[2] / dd4hep::tesla;
    delete[] b_vector;
    // Define the scale for object color and size.
    double en_min = 0.0, en_max = 100.0;  // Clusters,
    double ptot_min = 0.0, ptot_max = 25.0;  // Tracks.
    for (int i=0; i < rp_col->getNumberOfElements(); ++i) {
      EVENT::ReconstructedParticle* rp =
        static_cast<EVENT::ReconstructedParticle*>(rp_col->getElementAt(i));
      // Conveniently store the necessary particle information in a struct.
      viewer_util::AnyParticle ap(rp);
      ev_sum += ap;
      // Reference point of the particle. So far, only origin is implemented.
      std::vector<double> ref_pt = {0, 0, 0};
      // Define the layers the rp's components are printed to (E-dependent).
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
      int type = rp->getType();
      int color = returnTrackColor(type);
      int track_size = 1;

      DDMarlinCED::drawHelix(b_field_z, ap.charge, ref_pt[0], ref_pt[1], ref_pt[2], ap.x, ap.y, ap.z, tpc_layer,
          track_size, color, 0.0, helix_max_r, helix_max_z, rp->id() );

      //				/** Draw momentum lines from the ip */
      char Mscale = 'b'; // 'b': linear, 'a': log
      int McolorMap = 2; //hot: 3
      int McolorSteps = 256;
      int Mcolor = returnRGBClusterColor(ap.getP(), ptot_min, ptot_max, McolorSteps, Mscale, McolorMap);
      int LineSize = 1;
      float momScale = 100;
      ced_line_ID(ref_pt[0], ref_pt[1], ref_pt[2], momScale*ap.x, momScale*ap.y, momScale*ap.z, mom_layer, LineSize, Mcolor, rp->id()); //the right id?




      ClusterVec clusterVec = rp->getClusters();
      int nClusters = (int)clusterVec.size();
      if (nClusters > 0 ) {
	// std::cout 	<< "nCluster > 0" << std::endl;
	// std::cout 	<<  clusterVec.size() <<std::endl;

	Cluster * cluster = clusterVec[0];

	double * center  = new double[3];
	float * center_r  = new float[3];
	center[0] = cluster->getPosition()[0];
	center[1] = cluster->getPosition()[1];
	center[2] = cluster->getPosition()[2];

	center_r[0] = (float)cluster->getPosition()[0];
	center_r[1] = (float)cluster->getPosition()[1];
	center_r[2] = (float)cluster->getPosition()[2];

	float phi = cluster->getIPhi();
	float theta = cluster->getITheta();

	if( phi ==0. && theta==0.){
	  // use the cluster postion:
	  theta = atan( sqrt( center[0]*center[0] + center[1]*center[1] ) / center[2]  ) ;
	  phi = atan2( center[1] , center[0] ) ;
	}

	float eneCluster = cluster->getEnergy();

  double rotate[] = {0.0, 0.0, 0.0};
	// for particles centered at the origin
	if (ref_pt == std::vector<double>{0, 0, 0}){
	  rotate[0] = 0.0;
	  rotate[1] = theta*180/M_PI;
	  rotate[2] = phi*180/M_PI;
	}
	else {
	  streamlog_out( DEBUG3 )	<< "This is not yet implemented" << std::endl;
	  return;
	}


	double sizes[] = {100.0, 100.0, 100.0};
	int scale_z = 4;
	sizes[0] = returnClusterSize(eneCluster, en_min, en_max);
	sizes[1] = sizes[0];
	sizes[2] = scale_z*returnClusterSize(eneCluster, en_min, en_max);

	char scale = 'a'; // 'b': linear, 'a': log
	int colorMap = 3; //jet: 3
	int colorSteps = 256;
	int color = returnRGBClusterColor(eneCluster, en_min, en_max, colorSteps, scale, colorMap);

	//	int hit_type = 1 | HIT_LAYER;

	int cylinder_sides = 30;
	//problem with fisheye view
	ced_geocylinder_r(sizes[0]/2, sizes[2], center, rotate, cylinder_sides, color, clu_layer);



	//ced_hit(center[0],center[1],center[2], hit_type, (int)(sqrt(2)*sizes[0]/4), color);
	ced_hit_ID(center[0],center[1],center[2], 1, hit_layer, (int)(sqrt(2)*sizes[0]/4), color, cluster->id()); //hauke


	int transparency = 0x66;
	int rgba = addAlphaChannelToColor(color, transparency);

	// Originally: BACKUP_LAYER.
  ced_cluellipse_r_ID((float)sizes[0], (float)sizes[2], center_r, rotate, clu_layer, rgba, cluster->id()); //hauke


	transparency = 0xCC;
	rgba = addAlphaChannelToColor(color, transparency);

	//ced_ellipsoid_r(sizes, center, rotate, BACKUP_LAYER2, rgba);
	// Originally: BACKUP_LAYER2.
	ced_ellipsoid_r_ID(sizes, center, rotate, clu_layer, rgba, cluster->id()); //hauke

	/*
	 * End points for the arrow
	 */

	int sizeLine = 2;
	//int ml_line = 9;
	float radius = (float)0.5*sizes[2];
	float arrow = 0.2*radius;
	int color_arrow = color + 0x009900; //colour of the needle
  ///std::vector<double> end
	float x_end = center[0] + (radius - arrow) * std::sin(theta)*std::cos(phi);
	float y_end = center[1] + (radius - arrow) * std::sin(theta)*std::sin(phi);
	float z_end = center[2] + (radius - arrow) * std::cos(theta);

	float xArrowEnd = center[0] + (radius) * std::sin(theta)*std::cos(phi);
	float yArrowEnd = center[1] + (radius) * std::sin(theta)*std::sin(phi);
	float zArrowEnd = center[2] + (radius) * std::cos(theta);

	// this is the direction arrow
	// Originally: CLU_LAYER.
	ced_line_ID(center[0], center[1], center[2], x_end, y_end, z_end, clu_layer, sizeLine, color,cluster->id()); //hauke
	// Originally: CLU_LAYER.
  ced_line_ID(x_end, y_end, z_end, xArrowEnd, yArrowEnd, zArrowEnd, clu_layer, sizeLine, color_arrow, cluster->id());//hauke

      }
    }


    char scale = 'a'; // 'b': linear, 'a': log
    int colorMap = 3; //jet: 3
    unsigned int colorSteps = 256;
    unsigned int ticks = 6; //middle
    DDDSTViewer::showLegendSpectrum(colorSteps, scale, colorMap, en_min, en_max, ticks);

    streamlog_out( MESSAGE ) << std::endl
			     << "Total Energy and Momentum Balance of Event" << std::endl
			     << "Energy = " << ev_sum.E
			     << " PX = " << ev_sum.x
			     << " PY = " << ev_sum.y
			     << " PZ = " << ev_sum.z
			     << " Charge = " << ev_sum.charge << std::endl
			     << std::endl ;

    streamlog_out( DEBUG2 ) << "Setup properties" << std::endl
			    << "B_Field = " << b_field_z << " T" << std::endl ;
  }

// Now we draw the jets (if any ).
  for (size_t i = 0; i < jet_col_names_.size(); ++i) {
    streamlog_out(DEBUG) << " drawing jets from collection "
      << jet_col_names_[i] << std::endl;
    LCCollection* col = event->getCollection( jet_col_names_[i] );
    int nelem = col->getNumberOfElements();
    for (int j=0; j < nelem; ++j) {
      ReconstructedParticle * jet = dynamic_cast<ReconstructedParticle*>( col->getElementAt(j) );
      streamlog_out( DEBUG )  <<   "     - jet energy " << jet->getEnergy() << std::endl ;
      const double * mom = jet->getMomentum() ;
      //const double jet_ene = jet->getEnergy();
      gear::Vector3D v( mom[0], mom[1] , mom[2] ) ;
      int layer = 0;
      const ReconstructedParticleVec & pv = jet->getParticles();
      float pt_norm = 0.0;
      for (unsigned int k = 0; k<pv.size(); ++k) {
        const double * pm = pv[k]->getMomentum();
        gear::Vector3D pp( pm[0], pm[1] , pm[2] ) ;
        gear::Vector3D ju = v.unit();
        gear::Vector3D pt = pp - (ju.dot(pp))*ju;
        pt_norm += pt.r();
        layer = returnJetLayer(jet_col_names_[i]);
        int color = returnJetColor(jet_col_names_[i], j);
        int LineSize = 1;
        // Reference point of the particle. So far, only origin is implemented.
        std::vector<double> ref_pt = {0, 0, 0};
        float momScale = 100;
        int layerIp = returnIpLayer(jet_col_names_[i]);;
        ced_line_ID(ref_pt[0], ref_pt[1], ref_pt[2], momScale*pm[0], momScale*pm[1], momScale*pm[2], layerIp, LineSize, color, pv[k]->id()); //hauke
      }
      double center_c[3] = {0., 0., 0. };
      double rotation_c[3] = { 0.,  v.theta()*180./M_PI , v.phi()*180./M_PI };
      float rgba_color[4] = { float(0.2+0.2*i), float(0.2+0.2*i), float(1.0-0.25*i), float(0.3)};
      if (i==0){
        rgba_color[0] = 1.0;
        rgba_color[1] = 0.3;
        rgba_color[2] = 0.1;
        rgba_color[3] = 0.3;
      }
      double scale_pt = 20;
      double scale_mom = 25;
      double min_pt = 50;
      const double *pm=jet->getMomentum();
      int i;
      for(i=1;i<104;i++){
        ced_line_ID((i-1)/4*pm[0], (i-1)/4*pm[1], (i-1)/4*pm[2], i/4*pm[0], i/4*pm[1], i/4*pm[2], layer, 1, viewer_util::fromRGBAToInt(rgba_color), jet->id());
      }
      ced_hit_ID((i-1)/4*pm[0], (i-1)/4*pm[1], (i-1)/4*pm[2],
        0, layer, 2, viewer_util::fromRGBAToInt(rgba_color), jet->id()); //hauke
      ced_cone_r_ID( min_pt + scale_pt*pt_norm , scale_mom*v.r() , center_c, rotation_c, layer, rgba_color,jet->id()); //hauke
    }
  }

// My MC Particle addition:
// Adapted from the Jet loop above.
  streamlog_out(DEBUG) << " drawing MC from collection " << mc_col_name_ << std::endl ;
  EVENT::LCCollection* mc_col = nullptr;
  try {
    mc_col = event->getCollection(mc_col_name_);
  } catch (DataNotAvailableException &e) {
    streamlog_out(ERROR) << "MC collection " << mc_col_name_
      << " is not available!" << std::endl;
  }
  if (mc_col) {
    for (int j=0; j < mc_col->getNumberOfElements(); ++j) {
      EVENT::MCParticle * mcp = dynamic_cast<EVENT::MCParticle*>(mc_col->getElementAt(j) );
      // Only display the some of the MC particles. Also, define colors depending on the PDG.
      int i_col = 0;
      int mc_layer = kMCOther;
      // Draw all primary MC particles (with special colors for Higgs and tau).
      if (mcp->getParents().size() == 0) {
        if(mcp->getPDG() == 25) {
          i_col = 1;
          mc_layer = kMCHiggs;
        } else if(abs(mcp->getPDG()) == 15) {
          i_col = 2;
          mc_layer = kMCTau;
        // Tiny energy photons etc are not of interest.
        } else if(mcp->getEnergy() < 2) {
          continue;
        }
      // Draw the direct children of the Higgs.
      } else if (mcp->getParents()[0]->getPDG() == 25 && mcp->getPDG() != 25) {
        i_col = 3;
        mc_layer = kMCHiggsDecay;
      // Ignore all remaining MC particles.
      } else {
        continue;
      }
      streamlog_out( DEBUG )  <<   "     - mcp energy " << mcp->getEnergy() << std::endl ;
      const double * mom = mcp->getMomentum() ;
      gear::Vector3D v( mom[0], mom[1] , mom[2] ) ;
      //const EVENT::MCParticleVec & pv = mcp->getParticles();
      float pt_norm = 0.0;
      double center_c[3] = {0., 0., 0. };
      double rotation_c[3] = { 0.,  v.theta()*180./M_PI , v.phi()*180./M_PI };
      // This i was from the loop over mcp collections.
      float rgba_color[4] = { float(0.2+0.2*i_col), float(0.2+0.2*i_col), float(1.0-0.25*i_col), float(0.3)};
      if (i_col == 0) {
        rgba_color[0] = 1.0;
        rgba_color[1] = 0.3;
        rgba_color[2] = 0.1;
        rgba_color[3] = 0.3;
      }
      double scale_pt = 20;
      double scale_mom = 25;
      double min_pt = 50;
      const double *pm=mcp->getMomentum();
      int color = int(rgba_color[2]*(15*16+15)) + int(rgba_color[1]*(15*16+15))*16*16+ int(rgba_color[0]*(15*16+15))*16*16*16*16;
      // Trajectory of the MC tracks is not necessary.
      // But we need a line for picking.
      int tpc_hit = 104;
      for (int i = 1; i < tpc_hit; ++i) {
        ced_line_ID((i-1)/4*pm[0], (i-1)/4*pm[1], (i-1)/4*pm[2],
                  i/4*pm[0], i/4*pm[1], i/4*pm[2], mc_layer, 1, color, mcp->id());
      }
      ced_hit_ID((tpc_hit-1)/4*pm[0], (tpc_hit-1)/4*pm[1], (tpc_hit-1)/4*pm[2],
          0, mc_layer, 2, color, mcp->id());
      // Draw a cone.
      ced_cone_r_ID( min_pt + scale_pt*pt_norm , scale_mom*v.r() , center_c, rotation_c, mc_layer, rgba_color,mcp->id()); //hauke
    }
  }
  // This refreshes the view?
  DDMarlinCED::draw(this, wait_for_keyboard_);
}


/**
 * PRELIM */
int DDDSTViewer::returnTrackColor(int type) {

  int kcol = 0x999999; //default black - unknown

  if (type==211) 			kcol = 0x990000; //red - pi^+
  else if (type==-211) 	kcol = 0x996600; //orange - pi^-
  else if (type==22) 		kcol = 0x00ff00; //yellow - photon
  else if (type==11) 		kcol = 0x660066; //dark blue - e^-
  else if (type==-11) 	kcol = 0x660099; //violet - e^+
  else if (type==2112)	kcol = 0x99FFFF; //white - n0
  else if (type==-2112)	kcol = 0x99FFFF; //white n bar
  else {
    streamlog_out( DEBUG ) << "Unassigned type of colour: default" << std::endl;
  }
  return kcol;
}

void DDDSTViewer::showLegendSpectrum(const unsigned int color_steps, char scale, int colorMap, float en_min, float en_max, unsigned int ticks){

  const unsigned int numberOfColours = 3;
  unsigned int** rgb_matrix = new unsigned int*[color_steps]; //legend colour matrix;

  for (unsigned int i=0; i<color_steps; ++i){
    rgb_matrix[i] = new unsigned int[numberOfColours];
    ColorMap::selectColorMap(colorMap)(rgb_matrix[i], (float)i, 0.0, (float)color_steps);
    //			std::cout << "----------------------" << std::endl;
    //			std::cout << "i = " << i << std::endl;
    //			std::cout << "red = " << rgb_matrix[i][0] << std::endl;
    //			std::cout << "green = " << rgb_matrix[i][1] << std::endl;
    //			std::cout << "blue = " << rgb_matrix[i][2] << std::endl;
  }
  ced_legend(en_min, en_max, color_steps, rgb_matrix, ticks, scale);
}

int DDDSTViewer::returnRGBClusterColor(float energy, float cutoff_min, float cutoff_max, int color_steps, char scale, int color_map){
  if (color_map < 0 || color_map > 6) {
    streamlog_out(ERROR) << "Wrong color_map parameter!" << std::endl;
  }
  // Implicit conversion of the double returned from the function to an int in
  // the range of color_steps.
  int color_delta = viewer_util::convertScales(energy, cutoff_max, cutoff_min,
    color_steps, 0, viewer_util::ScaleMapping(scale));
  unsigned int rgb[] = {0, 0, 0};
  ColorMap::selectColorMap(color_map)(rgb, color_delta, 0, color_steps);
  int hex_color = ColorMap::RGB2HEX(rgb[0], rgb[1], rgb[2]);
  return hex_color;
}



int DDDSTViewer::returnClusterSize(float eneCluster, float cutoff_min, float cutoff_max){
  if (cutoff_min > cutoff_max) {
    streamlog_out(ERROR) << "cutoff_min < cutoff_max" << std::endl;
  }
  if (eneCluster < 0.0) {
    streamlog_out(ERROR) << "eneCluster is negative!" << std::endl;
  }
  if (cutoff_min < 0.0) {
    streamlog_out(ERROR) << "eneCluster is negative!" << std::endl;
  }

  int size = 0; //default size: zero

		// sizes
  int size_min = 10;
  int size_max = 120;

  // Input values in log-scale
  float log_ene = std::log(eneCluster+1);
  float log_min = std::log(cutoff_min+1);
  float log_max = std::log(cutoff_max+1);

  int size_steps = size_max - size_min; // default: 90 step sizes
  float log_delta = log_max - log_min;
  float log_step = log_delta/size_steps;

  int size_delta = (int) ((log_ene-log_min)/log_step); // which size bin does the value go to?

  if (size_delta >= size_steps){
    size = size_max;
  }
  else if (size_delta < size_min){
    size = size_min;
  }
  else {
    size = size_min + size_delta;
  }

  /**
   * Check the output */
  if (size <=0){
    std::cout << "Error in 'DDDSTViewer::returnClusterSize': return size is negative!" << std::endl;
  }
  return size;
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

int DDDSTViewer::addAlphaChannelToColor(int color, int alpha_channel){
  int rgba = 0xEE000000;
  rgba = (alpha_channel<<24) + color;
  return rgba;
}

int DDDSTViewer::returnJetColor(std::string jet_col_name, int col_number) {
  int white = 0x55555555;
  int black = 0x55000000;
  int red = 0x55550000;
  int green = 0x55008000;
  int yellow = 0x55555500;
  int fuchsia = 0x55550055;
  //int blue = 0x55000055;
  //int orange = 0x5555a500;
  //int violet = 0x55ee82ee;
  //int purple = 0x55800080;
  //int silver = 0x55c0c0c0;
  //int gold = 0x5555d700;
  //int gray = 0x55808080;
  //int aqua = 0x55005555;
  //int skyblue = 0x5587ceeb;
  //int lightblue = 0x55a558e6;
  //int khaki = 0x55f0e68c;
  // Assign the color based on the index in the parameter JetCollections in
  // the steering file.
  int color = 0x00000000;
  switch (col_number) {
  case 0:
    color = yellow;
    break;
  case 1:
    color = red;
    break;
  case 2:
    color = white;
    break;
  case 3:
    color = fuchsia;
    break;
  case 4:
    color = black;
    break;
  case 5:
    color = green;
    break;
  }
  return color;
}
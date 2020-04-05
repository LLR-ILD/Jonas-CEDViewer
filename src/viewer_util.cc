/**
*    @author: Thorben Quast, CERN Summer Student 2015 (some original
 *        functionality for porting CEDViewer to dd4hep detector geometry).
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
*/
// -- C++ STL headers.

// -- ROOT headers.

// -- LCIO headers.

// -- Marlin headers.
#include "marlin/Exceptions.h"

// -- Includes for detector drawing.
#include "DD4hep/DD4hepUnits.h"

// -- Header for this processor and other project-specific headers.
#include "viewer_util.h"

// -- Using-declarations and global constants.

// ----------------------------------------------------------------------------

// Struct for calculating the track length of a particle.
viewer_util::CalorimeterDrawParams viewer_util::getCalorimeterParameters(
    dd4hep::Detector& the_detector, std::string cal_name, bool self_call) {
  viewer_util::CalorimeterDrawParams params;
  if (self_call) cal_name[1] = tolower(cal_name[1]);
  dd4hep::DetElement calo;
  const std::vector<dd4hep::DetElement>& calorimeters =
      the_detector.detectors("calorimeter");
  for (size_t i = 0, n = calorimeters.size(); i < n; ++i) {
    if ((std::string) calorimeters[i].name() == cal_name) {
      calo = calorimeters[i];
      break;
    }
  }
  dd4hep::rec::LayeredCalorimeterData* calo_geo;
  try {
    calo_geo = calo.extension<dd4hep::rec::LayeredCalorimeterData>();
  } catch (std::runtime_error& e) {
    if (!self_call) {
        return getCalorimeterParameters(the_detector, cal_name, true);
    } else {
      streamlog_out(MESSAGE) << "MC Particles for " << cal_name
          << " cannot be drawn" << std::endl;
      params.delta_z = -1;  // No spatial extension --> No drawing.
      return params;
    }
  }
  params.r_inner = calo_geo->extent[0] / dd4hep::mm;
  params.delta_r = calo_geo->extent[1] / dd4hep::mm - params.r_inner;
  if (calo_geo->layoutType ==
  dd4hep::rec::LayeredCalorimeterData::BarrelLayout) {
    params.z_0 = 0.;
    params.delta_z = calo_geo->extent[3] / dd4hep::mm;
  } else {
    params.z_0 = calo_geo->extent[2] / dd4hep::mm;
    params.delta_z = (calo_geo->extent[3] / dd4hep::mm - params.z_0);
    // -z_0: We are interested in the full length, but CEDGeoTube only
    // requires half length as an argument.
  }
  return params;
}

// ----------------------------------------------------------------------------
// dd4hep draw helpers.

// Get the outer extents of the tracker.
double* viewer_util::getTrackerExtent(dd4hep::Detector& the_detector) {
  double* extent = new double[2];
  extent[0] = the_detector.constant<double>("tracker_region_rmax") / dd4hep::mm;
  extent[1] = the_detector.constant<double>("tracker_region_zmax") / dd4hep::mm;
  return extent;
}

// Get the outer extents of the yoke.
double* viewer_util::getYokeExtent(dd4hep::Detector& the_detector) {
  double* extent = new double[2];
  dd4hep::DetElement yoke;
  const std::vector< dd4hep::DetElement>& calorimeters =
      the_detector.detectors("calorimeter");
  for (size_t i = 0, n = calorimeters.size(); i < n; ++i) {
    std::string calo_name = calorimeters[i].name();
    yoke = calorimeters[i];
    bool is_yoke_barrel = (calo_name == "YokeBarrel");
    bool is_yoke_endcap = (calo_name == "YokeEndcap");
    if (!is_yoke_barrel && !is_yoke_endcap) continue;
    dd4hep::rec::LayeredCalorimeterData* yoke_geo;
    try {
      yoke_geo = yoke.extension<dd4hep::rec::LayeredCalorimeterData>();
    } catch (std::runtime_error& e) {
      streamlog_out(MESSAGE) << "MC Particles for " << calo_name
          << " cannot be drawn" << std::endl;
      extent[0] = extent[1] = 0;
      return extent;
    }
    if (is_yoke_barrel) extent[0] = yoke_geo->extent[1] / dd4hep::mm;
    if (is_yoke_endcap) extent[1] = yoke_geo->extent[3] / dd4hep::mm;
  }
  return extent;
}

// It suffices to perform the calculations in the first quadrant due to the
// detector's symmetry.
// The signs of the tracks' directions are ultimately determined by the momenta.
double viewer_util::calculateTrackLength(std::string type,
    dd4hep::Detector& the_detector,
    double x, double y, double z, double px, double py, double pz) {
  // X0 radiation length in ECAL.
  double rel_x0 = 0.5;
  viewer_util::CalorimeterDrawParams barrel, endcap;
  if (type == "ecal") {
    barrel = getCalorimeterParameters(the_detector, "ECalBarrel");
    endcap = getCalorimeterParameters(the_detector, "ECalEndcap");
  } else if(type == "hcal") {
    barrel = getCalorimeterParameters(the_detector, "HCalBarrel");
    endcap = getCalorimeterParameters(the_detector, "HCalEndcap");
  } else {
    barrel = getCalorimeterParameters(the_detector, "ECalBarrel");
    endcap = getCalorimeterParameters(the_detector, "ECalEndcap");
    rel_x0 = 0.;
  }
  // The case if the parameters could not be loaded properly.
  if (barrel.delta_z == -1 || endcap.delta_z == -1) return 0;
  double pt = sqrt(px*px + py*py);
  double pt_over_pz = pt / fabs(pz);
  double r = sqrt(x*x+y*y);
  if (r > barrel.r_inner || r > endcap.r_inner) return 0;
  double p = 2 * (px * x + py * y) / pt;
  double q_barrel = r*r - barrel.r_inner*barrel.r_inner;
  double dist_to_barrel_r = -p/2 + sqrt(p*p/4 - q_barrel);
  double q_endcap = r*r - endcap.r_inner * endcap.r_inner;
  double dist_to_endcap_r = -p/2 + sqrt(p*p/4 - q_endcap);
  double sign_pz = (pz >= 0) ? 1. : -1.;
  double dist_to_barrel_z = barrel.delta_z - sign_pz*z;
  double dist_to_endcap_z = endcap.z_0 - sign_pz*z;
  enum TrackThroughCalos {
    kEndsInBarrel = 0,
    kBarrelAndEndcap,
    kTroughCalGap,
    kEndcapOnly,
    kThroughEndcapReachesYoke,
  };
  double length;
  TrackThroughCalos track_through;
  // Direction: Traverses only barrel.
  if (pt_over_pz > (dist_to_barrel_r + barrel.delta_r) / dist_to_barrel_z) {
    track_through = kEndsInBarrel;
  // Direction: Touches both barrel and endcap.
  } else if (pt_over_pz > dist_to_barrel_r / dist_to_barrel_z) {
    // x_trav: path traversed in barrel, rotation symmetry is still assumed at
    // this point which is a valid approximation most of the times.
    double x_trav = (dist_to_barrel_z - dist_to_barrel_r / pt_over_pz)
            * sqrt(1 + pow(pt_over_pz, 2));
    // Traversed path in the barrel is larger than defined interaction path.
    if (x_trav > rel_x0 * barrel.delta_r) {
      track_through = kEndsInBarrel;
    // Particle is not absorbed in barrel but reaches the endcap.
    } else {
      length = (rel_x0 - x_trav / barrel.delta_r) * endcap.delta_z
             + dist_to_endcap_z * sqrt(1. + pow(pt_over_pz, 2));
      double length_r = length * pt_over_pz / sqrt(1. + pow(pt_over_pz, 2));
      if (length_r < (dist_to_endcap_r + endcap.delta_r)) {
        track_through = kBarrelAndEndcap; // Keep length calculated above.
      // Distance from z-axis exceeds endcap extension (e.g. if particle travels
      // through gap).
      } else {
        track_through = kTroughCalGap;
      }
    }
  } else if (pt_over_pz > dist_to_endcap_r / dist_to_endcap_z) {
    track_through = kEndcapOnly;
  } else if (pt_over_pz > dist_to_endcap_r /
  (dist_to_endcap_z + endcap.delta_z)) {
    double x_trav = (dist_to_endcap_r/pt_over_pz - dist_to_endcap_z)
                * sqrt(1+pow(pt_over_pz, 2));
    // Traversed path in endcap is larger than defined interaction path.
    if (x_trav > rel_x0*endcap.delta_z) {
      track_through = kEndcapOnly;
    } else {
      track_through = kThroughEndcapReachesYoke;
    }
  }
  switch (track_through) {
    case kEndsInBarrel:
      return barrel.delta_r * rel_x0
           + dist_to_barrel_r * sqrt(1. + pow(1./pt_over_pz, 2));
    case kEndcapOnly:
      return rel_x0 * endcap.delta_z
           + dist_to_endcap_z * sqrt(1. + pow(pt_over_pz,2));
    case kBarrelAndEndcap:
      return length;
    case kTroughCalGap:
      return (dist_to_endcap_r + endcap.delta_r)
            / pt_over_pz * sqrt(1. + pow(pt_over_pz, 2));
    case kThroughEndcapReachesYoke: // Same formula as the default case.
    // If none of the above cases applied, the track must have been to forward
    // to reach any calorimeter.
    default:
      return (dist_to_endcap_z + endcap.delta_z) * sqrt(1. + pow(pt_over_pz,2));
  }
}
// ----------------------------------------------------------------------------
// Functionality related to color.

int viewer_util::fromRGBAToInt(float* rgba_color) {
  return int(rgba_color[2]*(15*16+15))
       + int(rgba_color[1]*(15*16+15))*16*16
       + int(rgba_color[0]*(15*16+15))*16*16*16*16;
}

float* viewer_util::fromIntToRGBA(int hex_value) {
  float rgba[4];
  rgba[0] = ((hex_value >> 32) & 0xFF) / 255.0;  // Extract the transparency.
  rgba[1] = ((hex_value >> 16) & 0xFF) / 255.0;  // Extract the RR byte.
  rgba[2] = ((hex_value >> 8)  & 0xFF) / 255.0;  // Extract the GG byte.
  rgba[3] = ((hex_value)       & 0xFF) / 255.0;  // Extract the BB byte.
  return rgba;
}

// Convert the value laying inside the old scale to a value in a new scale,
// possibly with a non-linear mapping.
// Useful e.g. for conversion of an energy value to a color scale.
double viewer_util::convertScales(double value, double old_max, double old_min,
    double new_max, double new_min, viewer_util::ScaleMapping conv) {
  // Sanity checks for the input.
  if (new_max > new_min) {
    double temp = new_max;
    new_max = new_min;
    new_min = temp;
  }
  if (old_max > old_min) {
    double temp = old_max;
    old_max = old_min;
    old_min = temp;
  }
  if (old_min == old_max || new_min == new_max) {
    streamlog_out(ERROR) << "One of the scales has length 0!" << std::endl;
    return new_max;
  }
  // If the given value lies outside of the old scales boundary, we will just
  // set the converted value to the respective boundary value of the new scale.
  if (value > old_max) {
    return new_max;
  } else if (value < old_min) {
    return new_min;
  }
  if (value < 0 && conv == kLog) {
    streamlog_out(ERROR) << "When choosing a scale with logarithmic conversion "
        "the old minimum and maximum, as well as the old given value must be "
        "positive! Here we encountered" << std::endl
        << old_min << old_max << value << std::endl
        << "respectively. To not break, the highest value of the new scale is "
        "returned. Please reconsider your input to "
        "`convertScales`!" << std::endl;
    return new_max;
  }
  // Here the actual conversion takes place.
  switch (conv) {
    case kLog: default:
      // For this case, all the old-scale values are replaced by their
      // log versions.
      old_max = std::log(old_max + 1);
      old_min = std::log(old_min + 1);
      value   = std::log(value   + 1);
      break;
    case kLinear:
      // The old values are already in the form that is required.
      break;
  }
  return new_min + (value - old_min) / (old_max - old_min)
                                     * (new_max - new_min);
}
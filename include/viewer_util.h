/**
 *  Namespace collection functions that might in the future be used by multiple
 *  Viewer. Right now it only serves DDDSTViewer.
 *
 *  As this is a collection of methods, there is no global explanation. Instead,
 *  each method should be explained individually, or be clear from its name.
 *  It is advisable to check from time to time wether the methods described here
 *  are still in use in the project.
 *
 *    @author: Thorben Quast, CERN Summer Student 2015 (some original
 *        functionality for porting CEDViewer to dd4hep detector geometry).
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
 */

#ifndef _VIEWER_UTIL_H_
#define _VIEWER_UTIL_H_

// -- Includes for detector drawing.
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"

namespace viewer_util {

// Struct for calculating the track length of a particle.
struct CalorimeterDrawParams {
  double r_inner;
  double delta_r;
  double z_0;
  double delta_z;
};

// ----------------------------------------------------------------------------
// dd4hep draw helpers.

// Get the outer extents of the tracker.
double* getTrackerExtent(dd4hep::Detector& the_detector);

// Get the outer extents of the yoke.
double* getYokeExtent(dd4hep::Detector& the_detector);

// Get the relevant calorimeter parameters for track length calculations.
CalorimeterDrawParams getCalorimeterParameters(dd4hep::Detector& the_detector,
    std::string name, bool self_call = false);

// It suffices to perform the calculations in the first quadrant due to the
// detector's symmetry.
// The signs of the tracks' directions are ultimately determined by the momenta.
double calculateTrackLength(std::string type, dd4hep::Detector& the_detector,
    double x, double y, double z, double px, double py, double pz);


// ----------------------------------------------------------------------------
// Functions additional to those from DDCEDViewer.
int fromRGBAToInt(float* rgba_color);
}  // namespace viewer_util

#endif
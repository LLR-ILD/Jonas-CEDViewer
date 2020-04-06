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
// -- C++ STL headers.

// -- ROOT headers.

// -- LCIO headers.
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"

// -- Marlin headers.

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

// Get the relevant calorimeter parameters for track length calculations.
CalorimeterDrawParams getCalorimeterParameters(dd4hep::Detector& the_detector,
    std::string name, bool self_call = false);
// ----------------------------------------------------------------------------
// dd4hep draw helpers.

// Get the outer extents of the tracker.
double* getTrackerExtent(dd4hep::Detector& the_detector);

// Get the outer extents of the yoke.
double* getYokeExtent(dd4hep::Detector& the_detector);

// It suffices to perform the calculations in the first quadrant due to the
// detector's symmetry.
// The signs of the tracks' directions are ultimately determined by the momenta.
double calculateTrackLength(std::string type, dd4hep::Detector& the_detector,
    double x, double y, double z, double px, double py, double pz);
// ----------------------------------------------------------------------------
// Functions additional to those from DDCEDViewer.

// Facilitate simple particle arithmetics.
struct AnyParticle {
  double E = 0;
  double x = 0;
  double y = 0;
  double z = 0;
  int charge = 0;
  AnyParticle() {E = x = y = z = charge = 0;}
  AnyParticle(EVENT::ReconstructedParticle* rp) {
    const double* p = rp->getMomentum();
    x = p[0]; y = p[1], z = p[2];
    E = rp->getEnergy();
    charge = rp->getCharge();
  }
  AnyParticle(EVENT::MCParticle* rp) {
    const double* p = rp->getMomentum();
    x = p[0]; y = p[1], z = p[2];
    E = rp->getEnergy();
    charge = rp->getCharge();
  }
  AnyParticle operator+=(const AnyParticle rhs) {
    this->E += rhs.E;
    this->x += rhs.x;
    this->y += rhs.y;
    this->z += rhs.z;
    this->charge += rhs.charge;
    return *this;
  }
  double getP() {
    return sqrt(x*x + y*y + z*z);
  }
  double getTheta() {
    return atan(sqrt(x*x + y*y) / z);
  }
  double getPhi() {
    return atan2(y, x);
  }
};
// ----------------------------------------------------------------------------
// Functionality related to color.

int fromRGBAToInt(float* rgba_color);
void fromIntToRGBA(int hex_value, float (&rgba)[4]);

enum ScaleMapping {
  kLinear = 'b',
  kLog = 'a',
};

// Convert the value laying inside the old scale to a value in a new scale,
// possibly with a non-linear mapping.
// Useful e.g. for conversion of an energy value to a color scale.
double convertScales(double value, double old_max, double old_min,
    double new_max, double new_min=0, ScaleMapping conv=kLog);
}  // namespace viewer_util

#endif
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
    if (x == 0. && y == 0. && z == 0.) return 0;
    return atan2(sqrt(x*x + y*y), z);
  }
  double getPhi() {
    if (x == 0. && y == 0) return 0;
    return atan2(y, x);
  }
};
// ----------------------------------------------------------------------------
// Functionality related to color.

const int kHexPosAlpha = 16*16*16*16*16*16;
const int kHexPosR     = 16*16*16*16;
const int kHexPosG     = 16*16;
const int kHexPosB     = 1;
const int kByte = 16*16;
/**
 *  Utility class for color and conversion between:
 *    - hex (0x01122): rgb-like.
 *    - hexa (0x00112233): Including the alpha transparency-channel.
 *    - rgb ([0,1,2]): Three integers in [0, 255].
 *    - rgba ([9,1,2,3]) Four integers in [0, 255]. First entry is alpha.
 *    - rgb_float: Three floats in [0,1).
 *    - rgba_float: Four floats in [0,1). First entry is alpha.
 **/
class DSTColor {
 private:
  unsigned int hexa_val = 0x00555555;
 public:
  // Constructors.
  DSTColor() {this->hexa_val = 0x00555555;}
  DSTColor(int hex_in) {this->hexa_val = hex_in;}
  // Setters.
  void setHex(int hex_in) {this->hexa_val = hex_in % kHexPosAlpha;}
  void setHexa(int hexa_in) {this->hexa_val = hexa_in;}
  void setRGBFloat(float* rgb);
  void setRGBFloat(float r, float g, float b);
  void setRGBAFloat(float* rgba);
  void setRGBAFloat(float a, float r, float g, float b);
  void setRGB(unsigned int* rgb);
  void setRGB(unsigned int r, unsigned int g, unsigned int b);
  void setRGBA(unsigned int* rgba);
  void setRGBA(unsigned int a, unsigned int r, unsigned int g, unsigned int b);
  void setRGB(int* rgb) {setRGB((unsigned int*) rgb);}
  void setRGB(int r, int g, int b) {setRGB(r, g, b);}
  void setRGBA(int* rgba) {setRGBA((unsigned int*) rgba);}
  void setRGBA(int a, int r, int g, int b) {setRGBA(a, r, g, b);}
  // Getters.
  unsigned int hex() {return hexa_val % kHexPosAlpha;}
  unsigned int hexa() {return hexa_val;}
  void update_rgb(unsigned int* rgb);
  void update_rgba(unsigned int* rgba);
  void update_rgb(float* rgb);
  void update_rgba(float* rgba);
  // Transparency functionality.
  void addAlpha(unsigned int alpha) {this->hexa_val += alpha * kHexPosAlpha;}
  void addAlpha(int alpha) {addAlpha((unsigned int) alpha);}
  void addAlpha(float alpha) {addAlpha((unsigned int) alpha * kByte);}
  void setAlpha(unsigned int alpha) {
    this->hexa_val =  this->hexa_val % kHexPosAlpha + alpha * kHexPosAlpha;
  }
  void setAlpha(int alpha) {setAlpha((unsigned int) alpha);}
  void setAlpha(float alpha) {
    if (alpha <= 1) alpha *= kByte;
    setAlpha((unsigned int) alpha);
    }
};

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
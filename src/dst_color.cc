/**
*    Define the methods for the DSTColor helper class.
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
*/
#include "viewer_util.h"

// Setter.

void viewer_util::DSTColor::setRGBFloat(float* rgb) {
  unsigned int as_hex = int(rgb[0] * (kByte-1)) * kHexPosR
                      + int(rgb[1] * (kByte-1)) * kHexPosG
                      + int(rgb[2] * (kByte-1)) * kHexPosB;
  this->hexa_val = as_hex;
}

void viewer_util::DSTColor::setRGBFloat(float r, float g, float b) {
  unsigned int as_hex = int(r * (kByte-1)) * kHexPosR
                      + int(g * (kByte-1)) * kHexPosG
                      + int(b * (kByte-1)) * kHexPosB;
  this->hexa_val = as_hex;
}

void viewer_util::DSTColor::setRGBAFloat(float* rgba) {
  unsigned int as_hex = int(rgba[0] * (kByte-1)) * kHexPosR
                      + int(rgba[1] * (kByte-1)) * kHexPosG
                      + int(rgba[2] * (kByte-1)) * kHexPosB
                      + int(rgba[3] * (kByte-1)) * kHexPosAlpha;
  this->hexa_val = as_hex;
}

void viewer_util::DSTColor::setRGBAFloat(float a, float r, float g, float b) {
  unsigned int as_hex = int(r * (kByte-1)) * kHexPosR
                      + int(g * (kByte-1)) * kHexPosG
                      + int(b * (kByte-1)) * kHexPosB
                      + int(a * (kByte-1)) * kHexPosAlpha;
  this->hexa_val = as_hex;
}

void viewer_util::DSTColor::setRGB(unsigned int* rgb) {
  unsigned int as_hex = rgb[0] * kHexPosR
                      + rgb[1] * kHexPosG
                      + rgb[2] * kHexPosB;
  this->hexa_val = as_hex;
}

void viewer_util::DSTColor::setRGBA(unsigned int* rgba) {
  unsigned int as_hex = rgba[0] * kHexPosR
                      + rgba[1] * kHexPosG
                      + rgba[2] * kHexPosB
                      + rgba[3] * kHexPosAlpha;
  this->hexa_val = as_hex;
}

void viewer_util::DSTColor::setRGB(unsigned int r, unsigned int g,
    unsigned int b) {
  unsigned int as_hex = r * kHexPosR
                      + g * kHexPosG
                      + b * kHexPosB;
  this->hexa_val = as_hex;
}

void viewer_util::DSTColor::setRGBA(unsigned int a, unsigned int r,
    unsigned int g, unsigned int b) {
  unsigned int as_hex = r * kHexPosR
                      + g * kHexPosG
                      + b * kHexPosB
                      + a * kHexPosAlpha;
  this->hexa_val = as_hex;
}

// Getter.

void viewer_util::DSTColor::update_rgba(unsigned int* rgba) {
  // The old school method.
  ////rgba[3] = ((this->hexa_val >> 24) & 0xFF);  // Extract the alpha byte.
  ////rgba[0] = ((this->hexa_val >> 16) & 0xFF);  // Extract the RR byte.
  ////rgba[1] = ((this->hexa_val >> 8)  & 0xFF);  // Extract the GG byte.
  ////rgba[2] = ((this->hexa_val)       & 0xFF);  // Extract the BB byte.
  rgba[3] = (this->hexa_val / kHexPosAlpha) % kByte;
  rgba[0] = (this->hexa_val / kHexPosR    ) % kByte;
  rgba[1] = (this->hexa_val / kHexPosG    ) % kByte;
  rgba[2] = (this->hexa_val / kHexPosB    ) % kByte;
  // If no transparency was set, choose the fully opaque version of the color.
  if (rgba[3] < 1) rgba[3] = 0xFF;
}

void viewer_util::DSTColor::update_rgb(unsigned int* rgb) {
  rgb[0] = (this->hexa_val / kHexPosR    ) % kByte;
  rgb[1] = (this->hexa_val / kHexPosG    ) % kByte;
  rgb[2] = (this->hexa_val / kHexPosB    ) % kByte;
}

void viewer_util::DSTColor::update_rgb(float* rgb) {
  rgb[0] = float((this->hexa_val / kHexPosR    ) % kByte) / (kByte-1);
  rgb[1] = float((this->hexa_val / kHexPosG    ) % kByte) / (kByte-1);
  rgb[2] = float((this->hexa_val / kHexPosB    ) % kByte) / (kByte-1);
}

void viewer_util::DSTColor::update_rgba(float* rgba) {
  rgba[3] = float((this->hexa_val / kHexPosAlpha) % kByte) / (kByte-1);
  rgba[0] = float((this->hexa_val / kHexPosR    ) % kByte) / (kByte-1);
  rgba[1] = float((this->hexa_val / kHexPosG    ) % kByte) / (kByte-1);
  rgba[2] = float((this->hexa_val / kHexPosB    ) % kByte) / (kByte-1);
  // If no transparency was set, choose the fully opaque version of the color.
  if (rgba[3] < 1./(kByte-1)) rgba[3] = 1;
}
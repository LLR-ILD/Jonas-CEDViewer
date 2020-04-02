/**
 *  Event display for DST files.
 *
 *  Advanced functionality over the original DSTViewer.
 *  Amongst others, the detector geometry is ported to DD4hep.
 *
 *    @author: Szymon Daraszewicz, DESY summerstudent (original).
 *    @author:  iLCSoft/CEDViewer author list (DSTViewer).
 *    @author Jonas Kunath, LLR, CNRS, École Polytechnique, IPP.
 **/
// TODO: Clean up includes.
#ifndef DD_DST_VIEWER_H
#define DD_DST_VIEWER_H
// -- C++ STL headers.

// -- ROOT headers.

// -- LCIO headers.
#include "IMPL/LCCollectionVec.h"

// -- Marlin headers.
#include "marlin/Processor.h"

// -- Header for this processor and other project-specific headers.

class DDDSTViewer : public Processor {
 public:
  virtual Processor*  newProcessor() { return new DDDSTViewer(); }
  DDDSTViewer();
  virtual void init();
  virtual void processRunHeader(EVENT::LCRunHeader* run);
  virtual void processEvent(EVENT::LCEvent* event);

 protected:
  // -- Parameters registered in steering file.
  std::string rp_col_name_{};
  std::string mc_col_name_{};
  StringVec  jet_col_names_{};
  int wait_for_keyboard_{};
  float e_draw_cut_{};
  // -- Additional constants.
  int n_run_;
  int n_event_;
  // -- Individual member functions.
  void writeLayerDescription(void);


  int returnTrackColor(int type);
  int returnClusterColor(float eneCluster, float cutoff_min, float cutoff_max);
  void showLegendSpectrum(const unsigned int color_steps, char scale, int colorMap, float ene_min, float ene_max, unsigned int ticks);
  int returnRGBClusterColor(float eneCluster, float cutoff_min, float cutoff_max, int color_steps, char scale, int colorMap);
  int returnClusterSize(float eneCluster, float cutoff_min, float cutoff_max);
  int returnTrackSize(float type);
  int returnJetLayer(std::string jetColName);
  int returnIpLayer(std::string jetColName);
  int returnJetColor(std::string jetColName, int colNumber);
  int addAlphaChannelToColor(int color, int alphaChannel);
  float * returnConeColor(std::string jetColName);
};
#endif




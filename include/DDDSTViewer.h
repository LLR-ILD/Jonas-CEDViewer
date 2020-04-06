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
#ifndef DD_DST_VIEWER_H
#define DD_DST_VIEWER_H
// -- C++ STL headers.

// -- ROOT headers.

// -- LCIO headers.

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
  double e_draw_cut_{};
  StringVec detailled_drawn_detector_surfaces_{};
  bool is_drawn_surfaces_ = false;
  // -- Additional constants.
  int n_run_ = 0;
  int n_event_ = 0;
  // -- Additional member functions.
  void writeLayerDescription(void);
  int returnTrackColor(int particle_type);
  void showLegendSpectrum(char scale, int color_map,
      double ene_min, double ene_max, unsigned int ticks);
  int returnRGBClusterColor(double energy,
      double cutoff_min, double cutoff_max, char scale, int color_map);
  double returnClusterSize(double en_cluster,
      double cutoff_min, double cutoff_max);
  int returnJetLayer(std::string jet_col_name);
  int returnIpLayer(std::string jet_col_name);
  int returnJetColor(int col_number);
  int addTransparencyToColor(int color, int transparency);
};
#endif
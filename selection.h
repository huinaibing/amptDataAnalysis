#ifndef AMPT_DATA_ANALYSIS_SELECTION_H
#define AMPT_DATA_ANALYSIS_SELECTION_H

#include "analysisUtils.h"

// Compatibility wrappers used by older ROOT macros.
inline const ampt_analysis::AnalysisConfig &inclusiveSelectionConfig() {
  static const ampt_analysis::AnalysisConfig config = [] {
    auto value = ampt_analysis::defaultConfig();
    value.strictPtBounds = false;
    value.strictEtaBounds = false;
    return value;
  }();
  return config;
}

inline bool particleSelected(const Track &track) {
  return ampt_analysis::acceptsFlowTrack(track, inclusiveSelectionConfig());
}

inline bool cut4Pt(const Track &track) {
  return ampt_analysis::acceptsMeanPtTrack(track, 0,
                                           inclusiveSelectionConfig());
}

#endif // AMPT_DATA_ANALYSIS_SELECTION_H

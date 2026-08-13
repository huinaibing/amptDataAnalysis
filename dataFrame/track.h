#ifndef AMPT_DATA_ANALYSIS_TRACK_H
#define AMPT_DATA_ANALYSIS_TRACK_H

#include <cmath>
#include <limits>

struct Track {
  int eventID = 0;
  int nParticles = 0;
  double imp = 0.;

  int pdgPid = 0;
  double p_x = 0.;
  double p_y = 0.;
  double p_z = 0.;
  int sourceFile = -1;

  double GetPt() const { return std::sqrt(p_x * p_x + p_y * p_y); }

  double GetEta() const {
    const double p = std::sqrt(p_x * p_x + p_y * p_y + p_z * p_z);
    if (p == std::abs(p_z)) {
      return std::numeric_limits<double>::infinity();
    }
    return 0.5 * std::log((p + p_z) / (p - p_z));
  }

  double GetPhi() const { return std::atan2(p_y, p_x); }
};

#endif // AMPT_DATA_ANALYSIS_TRACK_H

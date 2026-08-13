#ifndef AMPT_DATA_ANALYSIS_EVENT_H
#define AMPT_DATA_ANALYSIS_EVENT_H

#include "track.h"

#include <cstddef>
#include <vector>

struct Event {
  struct PtMoments {
    std::size_t count = 0;
    double sum = 0.;
    double sumSquares = 0.;

    double mean() const {
      return count > 0 ? sum / static_cast<double>(count) : 0.;
    }

    double distinctPairMean() const {
      if (count < 2) {
        return 0.;
      }
      const double nPairs =
          static_cast<double>(count) * static_cast<double>(count - 1);
      return (sum * sum - sumSquares) / nPairs;
    }
  };

  int eventID = 0;
  double imp = 0.;

  std::vector<Track> particles;

  template <typename Predicate>
  PtMoments GetPtMoments(const Predicate &accept) const {
    PtMoments moments;
    for (const auto &track : particles) {
      if (!accept(track)) {
        continue;
      }
      const double pt = track.GetPt();
      ++moments.count;
      moments.sum += pt;
      moments.sumSquares += pt * pt;
    }
    return moments;
  }

  template <typename Predicate>
  double GetMeanPt(const Predicate &accept) const {
    return GetPtMoments(accept).mean();
  }

  template <typename Predicate>
  std::size_t nParticlesAfterCut(const Predicate &accept) const {
    return GetPtMoments(accept).count;
  }

  template <typename Predicate>
  double GetPtSquareAve(const Predicate &accept) const {
    return GetPtMoments(accept).distinctPairMean();
  }
};

#endif // AMPT_DATA_ANALYSIS_EVENT_H

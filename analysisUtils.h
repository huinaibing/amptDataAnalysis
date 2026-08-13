#ifndef AMPT_DATA_ANALYSIS_UTILS_H
#define AMPT_DATA_ANALYSIS_UTILS_H

#include "analysisConfig.h"
#include "dataFrame/fileName.h"
#include "eventManager.h"

#include "PWGCF/GenericFramework/Core/GFW.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <string>

namespace ampt_analysis {
struct CorrelationResult {
  double numerator = 0.;
  double pairs = 0.;

  bool hasPairs() const { return pairs != 0.; }
  double value() const { return hasPairs() ? numerator / pairs : 0.; }
  bool isPhysical() const { return hasPairs() && std::abs(value()) < 1.; }
};

struct EventSamples {
  Event::PtMoments charged;
  std::array<Event::PtMoments, 3> pid;

  const Event::PtMoments &forSpecies(Species species) const {
    return pid.at(static_cast<std::size_t>(species));
  }
};

inline void addPt(Event::PtMoments &moments, double pt) {
  ++moments.count;
  moments.sum += pt;
  moments.sumSquares += pt * pt;
}

inline CorrelationResult calculateCorrelation(GFW &gfw,
                                              const GFW::CorrConfig &config) {
  return {gfw.Calculate(config, 0, false).real(),
          gfw.Calculate(config, 0, true).real()};
}

inline CorrelationResult calculateCorrelation(GFW &gfw,
                                              const CorrConfigManager &manager,
                                              CorrType type) {
  return calculateCorrelation(gfw, manager.Get(type));
}

inline CorrelationResult
calculateCombinedCorrelation(GFW &gfw, const CorrConfigManager &manager,
                             CorrType typeA, CorrType typeB) {
  const CorrelationResult a = calculateCorrelation(gfw, manager, typeA);
  const CorrelationResult b = calculateCorrelation(gfw, manager, typeB);
  return {a.numerator + b.numerator, a.pairs + b.pairs};
}

inline double centralityFromImpactParameter(double impactParameter,
                                            const AnalysisConfig &config) {
  return config.centralityPi * impactParameter * impactParameter /
         config.centralityReferenceRadiusSquared;
}

inline bool acceptsFlowTrack(const Track &track, const AnalysisConfig &config) {
  const bool withinEta = config.strictEtaBounds
                             ? std::abs(track.GetEta()) < config.flowEtaMax
                             : std::abs(track.GetEta()) <= config.flowEtaMax;
  return withinEta &&
         config.chargedPt.contains(track.GetPt(), config.strictPtBounds);
}

inline bool acceptsMeanPtTrack(const Track &track, int absPdg,
                               const AnalysisConfig &config) {
  const bool withinEta = config.strictEtaBounds
                             ? std::abs(track.GetEta()) < config.meanPtEtaMax
                             : std::abs(track.GetEta()) <= config.meanPtEtaMax;
  if (!withinEta || !ptRangeForPdg(absPdg, config)
                         .contains(track.GetPt(), config.strictPtBounds)) {
    return false;
  }
  return absPdg == 0 || std::abs(track.pdgPid) == absPdg;
}

inline Event::PtMoments getMeanPtMoments(const Event &event, int absPdg,
                                         const AnalysisConfig &config) {
  return event.GetPtMoments([&](const Track &track) {
    return acceptsMeanPtTrack(track, absPdg, config);
  });
}

/**
 * Build all event-wise mean-pT samples and GFW Q-vectors in one particle pass.
 * The flow and mean-pT eta intervals remain independent.
 */
inline EventSamples fillGfwAndCollectSamples(GFW &gfw, const Event &event,
                                             const AnalysisConfig &config) {
  EventSamples samples;
  for (const auto &track : event.particles) {
    const double pt = track.GetPt();
    const double eta = track.GetEta();
    const int absPdg = std::abs(track.pdgPid);
    const SpeciesDefinition *species = findSpecies(absPdg);

    const bool withinMeanPtEta = config.strictEtaBounds
                                     ? std::abs(eta) < config.meanPtEtaMax
                                     : std::abs(eta) <= config.meanPtEtaMax;
    if (withinMeanPtEta) {
      if (config.chargedPt.contains(pt, config.strictPtBounds)) {
        addPt(samples.charged, pt);
      }
      if (species &&
          ptRangeForPdg(absPdg, config).contains(pt, config.strictPtBounds)) {
        addPt(samples.pid.at(static_cast<std::size_t>(species->species)), pt);
      }
    }

    const bool withinFlowEta = config.strictEtaBounds
                                   ? std::abs(eta) < config.flowEtaMax
                                   : std::abs(eta) <= config.flowEtaMax;
    if (!withinFlowEta) {
      continue;
    }

    if (config.chargedPt.contains(pt, config.strictPtBounds)) {
      gfw.Fill(eta, 0, track.GetPhi(), 1., Mask::kRef);
    }
    if (species &&
        ptRangeForPdg(absPdg, config).contains(pt, config.strictPtBounds)) {
      gfw.Fill(eta, 0, track.GetPhi(), 1.,
               species->mask | species->overlapMask);
    }
  }
  return samples;
}

inline std::size_t limitedConfigCount(std::size_t available, int maxConfigs) {
  if (maxConfigs < 0) {
    return available;
  }
  return std::min(available, static_cast<std::size_t>(maxConfigs));
}

inline int limitedFileCount(int configured, int maxFilesPerConfig) {
  return maxFilesPerConfig < 0 ? configured
                               : std::min(configured, maxFilesPerConfig);
}

template <typename Callback>
Long64_t forEachConfiguredEvent(const std::string &configFile,
                                int maxFilesPerConfig, int maxConfigs,
                                Callback &&callback) {
  if (maxFilesPerConfig == 0 || maxConfigs == 0) {
    return 0;
  }

  const auto configurations = CentConfigReader::load(configFile);
  const std::size_t nConfigurations =
      limitedConfigCount(configurations.size(), maxConfigs);
  Long64_t processedEvents = 0;

  for (std::size_t i = 0; i < nConfigurations; ++i) {
    const auto &configuration = configurations[i];
    const int nFiles =
        limitedFileCount(configuration.n_files, maxFilesPerConfig);
    if (nFiles <= 0) {
      continue;
    }

    AMPTEventReader reader(configuration.path, nFiles);
    for (const auto &event : reader) {
      callback(event, configuration);
      ++processedEvents;
    }
  }
  return processedEvents;
}
} // namespace ampt_analysis

#endif // AMPT_DATA_ANALYSIS_UTILS_H

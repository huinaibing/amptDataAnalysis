#ifndef AMPT_DATA_ANALYSIS_UTILS_H
#define AMPT_DATA_ANALYSIS_UTILS_H

#include "analysisConfig.h"
#include "dataFrame/fileName.h"
#include "eventManager.h"

#include "PWGCF/GenericFramework/Core/GFW.h"

#include "TDatabasePDG.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>

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

/**
 * Require a known PDG particle with non-zero electric charge. Unknown codes
 * are rejected so neutral or unidentified entries cannot enter charged sums.
 */
inline bool isChargedPdg(int pdgPid) {
  if (pdgPid == 0) {
    return false;
  }

  static std::unordered_map<int, bool> chargeCache;
  const auto cached = chargeCache.find(pdgPid);
  if (cached != chargeCache.end()) {
    return cached->second;
  }

  const TParticlePDG *particle =
      TDatabasePDG::Instance()->GetParticle(pdgPid);
  const bool isCharged = particle != nullptr && particle->Charge() != 0.;
  chargeCache.emplace(pdgPid, isCharged);
  return isCharged;
}

inline bool acceptsFlowTrack(const Track &track, const AnalysisConfig &config) {
  const bool withinEta = config.strictEtaBounds
                             ? std::abs(track.GetEta()) < config.flowEtaMax
                             : std::abs(track.GetEta()) <= config.flowEtaMax;
  return isChargedPdg(track.pdgPid) && withinEta &&
         config.chargedPt.contains(track.GetPt(), config.strictPtBounds);
}

inline bool acceptsMeanPtTrack(const Track &track, int absPdg,
                               const AnalysisConfig &config) {
  const bool withinEta = config.strictEtaBounds
                             ? std::abs(track.GetEta()) < config.meanPtEtaMax
                             : std::abs(track.GetEta()) <= config.meanPtEtaMax;
  if (!isChargedPdg(track.pdgPid) || !withinEta ||
      !ptRangeForPdg(absPdg, config)
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
    const bool isCharged = isChargedPdg(track.pdgPid);
    const double pt = track.GetPt();
    const double eta = track.GetEta();
    const int absPdg = std::abs(track.pdgPid);
    const SpeciesDefinition *species = findSpecies(absPdg);

    const bool withinMeanPtEta = config.strictEtaBounds
                                     ? std::abs(eta) < config.meanPtEtaMax
                                     : std::abs(eta) <= config.meanPtEtaMax;
    if (withinMeanPtEta) {
      if (isCharged &&
          config.chargedPt.contains(pt, config.strictPtBounds)) {
        addPt(samples.charged, pt);
      }
      if (isCharged && species &&
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

    if (isCharged && config.chargedPt.contains(pt, config.strictPtBounds)) {
      gfw.Fill(eta, 0, track.GetPhi(), 1., Mask::kRef);
    }
    if (isCharged && species &&
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

struct EventIdentity {
  std::size_t configIndex;
  int fileNumber;
  std::uint64_t eventIndex;
};

inline std::uint64_t mixEventKey(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31);
}

inline double bootstrapRandomValue(unsigned int seed,
                                   const EventIdentity &identity) {
  std::uint64_t key = mixEventKey(seed);
  key = mixEventKey(key ^ identity.configIndex);
  key = mixEventKey(key ^ static_cast<std::uint64_t>(identity.fileNumber));
  key = mixEventKey(key ^ identity.eventIndex);
  return static_cast<double>(key >> 11) * (1.0 / 9007199254740992.0);
}

template <typename Callback>
Long64_t forEachConfiguredEvent(const std::string &configFile,
                                int maxFilesPerConfig, int maxConfigs,
                                int shardIndex, int shardCount,
                                Callback &&callback) {
  if (shardCount <= 0 || shardIndex < 0 || shardIndex >= shardCount) {
    throw std::invalid_argument("shardIndex must be within shardCount");
  }
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
    if (nFiles <= shardIndex) {
      continue;
    }

    AMPTEventReader reader(configuration.path, nFiles, shardIndex, shardCount);
    int currentFile = -1;
    std::uint64_t eventIndex = 0;
    for (const auto &event : reader) {
      if (event.sourceFile != currentFile) {
        currentFile = event.sourceFile;
        eventIndex = 0;
      }
      callback(event, configuration, EventIdentity{i, currentFile, eventIndex++});
      ++processedEvents;
    }
  }
  return processedEvents;
}

template <typename Callback>
Long64_t forEachConfiguredEvent(const std::string &configFile,
                                int maxFilesPerConfig, int maxConfigs,
                                Callback &&callback) {
  return forEachConfiguredEvent(
      configFile, maxFilesPerConfig, maxConfigs, 0, 1,
      [&](const Event &event, const CentralityConfig &configuration,
          const EventIdentity &) { callback(event, configuration); });
}
} // namespace ampt_analysis

#endif // AMPT_DATA_ANALYSIS_UTILS_H

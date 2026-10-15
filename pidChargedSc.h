#ifndef AMPT_PID_CHARGED_SC_H
#define AMPT_PID_CHARGED_SC_H

#include "analysisUtils.h"

#include "TH1D.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include <array>
#include <cmath>
#include <vector>

namespace ampt_analysis {
namespace pid_charged_sc {

inline void fillGfw(GFW &gfw, const Event &event,
                    const AnalysisConfig &config) {
  gfw.Clear();
  const auto &settings = config.pidChargedScOutput;
  for (const auto &track : event.particles) {
    const double pt = track.GetPt();
    const double eta = track.GetEta();
    const double phi = track.GetPhi();
    if (!std::isfinite(pt) || !std::isfinite(eta) || !std::isfinite(phi) ||
        !isChargedPdg(track.pdgPid)) {
      continue;
    }
    const double absEta = std::abs(eta);
    const bool acceptedFullEta = settings.strictEtaBounds
                                    ? absEta < settings.etaMax
                                    : absEta <= settings.etaMax;
    if (!acceptedFullEta) {
      continue;
    }
    const bool acceptedGapEta = eta != 0. &&
        (settings.strictEtaBounds ? absEta > settings.innerAbsEta
                                  : absEta >= settings.innerAbsEta);
    const bool charged = config.chargedPt.contains(pt, settings.strictPtBounds);
    int mask = charged ? PidChargedScCorrConfigManager::kFullChargedMask : 0;
    if (charged && acceptedGapEta) {
      mask |= Mask::kRef;
    }
    const int absPdg = std::abs(track.pdgPid);
    const auto *species = findSpecies(absPdg);
    if (species && ptRangeForPdg(absPdg, config).contains(
                       pt, settings.strictPtBounds)) {
      mask |= PidChargedScCorrConfigManager::kFullPidMasks.at(
          static_cast<std::size_t>(species->species));
      if (acceptedGapEta) {
        mask |= species->mask;
        if (charged) {
          mask |= species->overlapMask;
        }
      }
    }
    if (mask != 0) {
      // Unit track weights for AMPT truth; event weights are the distinct
      // tuple denominators returned by GFW::Calculate(config, 0, true).
      gfw.Fill(eta, 0, phi, 1., mask);
    }
  }
}

struct Moments {
  // CorrelationResult::pairs stores the distinct tuple count (2 or 4 tracks).
  CorrelationResult pidTwo;
  CorrelationResult chargedTwo;
  CorrelationResult mixedFour;

  bool valid() const {
    for (const auto &moment : {pidTwo, chargedTwo, mixedFour}) {
      if (!(moment.pairs > 0.) || !std::isfinite(moment.pairs) ||
          !std::isfinite(moment.numerator)) {
        return false;
      }
    }
    return true;
  }
};

inline Moments calculateMoments(GFW &gfw,
                                const PidChargedScCorrConfigManager &manager,
                                std::size_t species) {
  return {calculateCorrelation(gfw, manager.PidTwo(species)),
          calculateCorrelation(gfw, manager.ChargedTwo()),
          calculateCorrelation(gfw, manager.MixedFour(species))};
}

struct ChargedMoments {
  CorrelationResult two;
  CorrelationResult four;

  bool valid() const {
    for (const auto &moment : {two, four}) {
      if (!(moment.pairs > 0.) || !std::isfinite(moment.pairs) ||
          !std::isfinite(moment.numerator)) {
        return false;
      }
    }
    return true;
  }
};

inline ChargedMoments calculateChargedMoments(
    GFW &gfw, const PidChargedScCorrConfigManager &manager) {
  return {calculateCorrelation(gfw, manager.ChargedTwo()),
          calculateCorrelation(gfw, manager.ChargedFour())};
}

struct FlowVarianceMoments {
  CorrelationResult gapTwo;
  CorrelationResult fullTwo;
  CorrelationResult fullFour;

  bool valid() const {
    for (const auto &moment : {gapTwo, fullTwo, fullFour}) {
      if (!(moment.pairs > 0.) || !std::isfinite(moment.pairs) ||
          !std::isfinite(moment.numerator)) {
        return false;
      }
    }
    return true;
  }
};

inline FlowVarianceMoments calculateFlowVarianceMoments(
    GFW &gfw, const PidChargedScCorrConfigManager &manager,
    const CorrelationResult &gapTwo, std::size_t species) {
  return {gapTwo, calculateCorrelation(gfw, manager.FullTwo(species)),
          calculateCorrelation(gfw, manager.FullFour(species))};
}

// A separate copy of the existing gap correlator keeps all three variance
// ingredients on one species-specific event sample, without imposing a
// four-POI requirement on the existing SC profiles. Store raw <4>, not c2{4}.
struct FlowVarianceProfiles {
  TProfile gapTwo;
  TProfile fullTwo;
  TProfile fullFour;
  TProfile2D gapTwoByGroup;
  TProfile2D fullTwoByGroup;
  TProfile2D fullFourByGroup;
  TH1D validEvents;

  FlowVarianceProfiles(const std::vector<double> &cent,
                       const std::vector<double> &groups)
      : gapTwo("pTwoGapForVariance", ";Centrality (%);<2>_{gap}",
               static_cast<int>(cent.size()) - 1, cent.data()),
        fullTwo("pTwoFull", ";Centrality (%);<2>_{full}",
                static_cast<int>(cent.size()) - 1, cent.data()),
        fullFour("pFourFull", ";Centrality (%);<4>_{full}",
                 static_cast<int>(cent.size()) - 1, cent.data()),
        gapTwoByGroup("pTwoGapForVarianceByGroup",
                      ";Centrality (%);Subsample;<2>_{gap}",
                      static_cast<int>(cent.size()) - 1, cent.data(),
                      static_cast<int>(groups.size()) - 1, groups.data()),
        fullTwoByGroup("pTwoFullByGroup", ";Centrality (%);Subsample;<2>_{full}",
                       static_cast<int>(cent.size()) - 1, cent.data(),
                       static_cast<int>(groups.size()) - 1, groups.data()),
        fullFourByGroup("pFourFullByGroup", ";Centrality (%);Subsample;<4>_{full}",
                        static_cast<int>(cent.size()) - 1, cent.data(),
                        static_cast<int>(groups.size()) - 1, groups.data()),
        validEvents("hFlowVarianceValidEvents", ";Centrality (%);Valid events",
                    static_cast<int>(cent.size()) - 1, cent.data()) {
    for (TH1 *histogram : std::array<TH1 *, 7>{
             &gapTwo, &fullTwo, &fullFour, &gapTwoByGroup, &fullTwoByGroup,
             &fullFourByGroup, &validEvents}) {
      histogram->SetDirectory(nullptr);
      histogram->Sumw2();
    }
  }

  void fill(double centrality, double group, const FlowVarianceMoments &moments,
            const PidChargedScOutputConfig &settings) {
    if (!moments.valid()) {
      return;
    }
    const std::array<CorrelationResult, 3> values{
        moments.gapTwo, moments.fullTwo, moments.fullFour};
    const std::array<TProfile *, 3> nominal{&gapTwo, &fullTwo, &fullFour};
    const std::array<TProfile2D *, 3> grouped{
        &gapTwoByGroup, &fullTwoByGroup, &fullFourByGroup};
    for (std::size_t i = 0; i < values.size(); ++i) {
      const double weight = settings.flowVarianceEventWeight == "full_quadruplets"
                                ? moments.fullFour.pairs : values[i].pairs;
      nominal[i]->Fill(centrality, values[i].value(), weight);
      grouped[i]->Fill(centrality, group, values[i].value(), weight);
    }
    validEvents.Fill(centrality);
  }

  void write() {
    gapTwo.Write();
    fullTwo.Write();
    fullFour.Write();
    gapTwoByGroup.Write();
    fullTwoByGroup.Write();
    fullFourByGroup.Write();
    validEvents.Write();
  }
};

// Independent charged baseline: both moments use events with at least two
// charged tracks in each eta side. PID availability does not restrict it.
struct ChargedProfiles {
  TProfile two;
  TProfile four;
  TProfile2D twoByGroup;
  TProfile2D fourByGroup;
  TH1D validEvents;
  TH1D sumPairs;
  TH1D sumQuadruplets;

  ChargedProfiles(const std::vector<double> &cent,
                  const std::vector<double> &groups)
      : two("pTwoCharged", ";Centrality (%);<2>_{charged}",
            static_cast<int>(cent.size()) - 1, cent.data()),
        four("pFourCharged", ";Centrality (%);<4>_{charged}",
             static_cast<int>(cent.size()) - 1, cent.data()),
        twoByGroup("pTwoChargedByGroup",
                   ";Centrality (%);Subsample;<2>_{charged}",
                   static_cast<int>(cent.size()) - 1, cent.data(),
                   static_cast<int>(groups.size()) - 1, groups.data()),
        fourByGroup("pFourChargedByGroup",
                    ";Centrality (%);Subsample;<4>_{charged}",
                    static_cast<int>(cent.size()) - 1, cent.data(),
                    static_cast<int>(groups.size()) - 1, groups.data()),
        validEvents("hValidEvents", ";Centrality (%);Valid events",
                    static_cast<int>(cent.size()) - 1, cent.data()),
        sumPairs("hSumChargedPairs", ";Centrality (%);Distinct charged pairs",
                 static_cast<int>(cent.size()) - 1, cent.data()),
        sumQuadruplets("hSumChargedQuadruplets",
                       ";Centrality (%);Distinct charged quadruplets",
                       static_cast<int>(cent.size()) - 1, cent.data()) {
    for (TH1 *histogram : std::array<TH1 *, 7>{
             &two, &four, &twoByGroup, &fourByGroup, &validEvents,
             &sumPairs, &sumQuadruplets}) {
      histogram->SetDirectory(nullptr);
      histogram->Sumw2();
    }
  }

  void fill(double centrality, double group, const ChargedMoments &moments,
            const PidChargedScOutputConfig &settings) {
    if (!moments.valid()) {
      return;
    }
    const double fourWeight = moments.four.pairs;
    const double twoWeight = settings.eventWeight == "mixed_quadruplets"
                                 ? fourWeight : moments.two.pairs;
    two.Fill(centrality, moments.two.value(), twoWeight);
    four.Fill(centrality, moments.four.value(), fourWeight);
    twoByGroup.Fill(centrality, group, moments.two.value(), twoWeight);
    fourByGroup.Fill(centrality, group, moments.four.value(), fourWeight);
    validEvents.Fill(centrality);
    sumPairs.Fill(centrality, moments.two.pairs);
    sumQuadruplets.Fill(centrality, moments.four.pairs);
  }

  void write() {
    two.Write();
    four.Write();
    twoByGroup.Write();
    fourByGroup.Write();
    validEvents.Write();
    sumPairs.Write();
    sumQuadruplets.Write();
  }
};

// ROOT profiles preserve weighted sums for merging; nonlinear SC/NSC and
// bootstrap uncertainties are constructed only after all shards are merged.
struct Profiles {
  TProfile pidTwo;
  TProfile chargedTwo;
  TProfile mixedFour;
  TProfile2D pidTwoByGroup;
  TProfile2D chargedTwoByGroup;
  TProfile2D mixedFourByGroup;
  TH1D validEvents;
  TH1D sumPidPairs;
  TH1D sumChargedPairs;
  TH1D sumMixedQuadruplets;

  Profiles(const std::vector<double> &cent, const std::vector<double> &groups)
      : pidTwo("pTwoPID", ";Centrality (%);<2>_{PID}",
               static_cast<int>(cent.size()) - 1, cent.data()),
        chargedTwo("pTwoCharged", ";Centrality (%);<2>_{charged}",
                   static_cast<int>(cent.size()) - 1, cent.data()),
        mixedFour("pFourPIDCharged", ";Centrality (%);<4>_{PID,charged}",
                  static_cast<int>(cent.size()) - 1, cent.data()),
        pidTwoByGroup("pTwoPIDByGroup", ";Centrality (%);Subsample;<2>_{PID}",
                      static_cast<int>(cent.size()) - 1, cent.data(),
                      static_cast<int>(groups.size()) - 1, groups.data()),
        chargedTwoByGroup("pTwoChargedByGroup",
                          ";Centrality (%);Subsample;<2>_{charged}",
                          static_cast<int>(cent.size()) - 1, cent.data(),
                          static_cast<int>(groups.size()) - 1, groups.data()),
        mixedFourByGroup("pFourPIDChargedByGroup",
                         ";Centrality (%);Subsample;<4>_{PID,charged}",
                         static_cast<int>(cent.size()) - 1, cent.data(),
                         static_cast<int>(groups.size()) - 1, groups.data()),
        validEvents("hValidEvents", ";Centrality (%);Valid events",
                    static_cast<int>(cent.size()) - 1, cent.data()),
        sumPidPairs("hSumPIDPairs", ";Centrality (%);Distinct PID pairs",
                    static_cast<int>(cent.size()) - 1, cent.data()),
        sumChargedPairs("hSumChargedPairs",
                        ";Centrality (%);Distinct charged pairs",
                        static_cast<int>(cent.size()) - 1, cent.data()),
        sumMixedQuadruplets("hSumMixedQuadruplets",
                            ";Centrality (%);Distinct mixed quadruplets",
                            static_cast<int>(cent.size()) - 1, cent.data()) {
    const std::array<TH1 *, 10> histograms{
        &pidTwo, &chargedTwo, &mixedFour, &pidTwoByGroup, &chargedTwoByGroup,
        &mixedFourByGroup, &validEvents, &sumPidPairs, &sumChargedPairs,
        &sumMixedQuadruplets};
    for (TH1 *histogram : histograms) {
      histogram->SetDirectory(nullptr);
      histogram->Sumw2();
    }
  }

  void fill(double centrality, double group, const Moments &moments,
            const PidChargedScOutputConfig &settings) {
    if (!moments.valid()) {
      return;
    }
    const bool commonWeight = settings.eventWeight == "mixed_quadruplets";
    const double fourWeight = moments.mixedFour.pairs;
    const double pidWeight = commonWeight ? fourWeight : moments.pidTwo.pairs;
    const double chargedWeight =
        commonWeight ? fourWeight : moments.chargedTwo.pairs;
    pidTwo.Fill(centrality, moments.pidTwo.value(), pidWeight);
    chargedTwo.Fill(centrality, moments.chargedTwo.value(), chargedWeight);
    mixedFour.Fill(centrality, moments.mixedFour.value(), fourWeight);
    pidTwoByGroup.Fill(centrality, group, moments.pidTwo.value(), pidWeight);
    chargedTwoByGroup.Fill(centrality, group, moments.chargedTwo.value(),
                           chargedWeight);
    mixedFourByGroup.Fill(centrality, group, moments.mixedFour.value(), fourWeight);
    validEvents.Fill(centrality);
    sumPidPairs.Fill(centrality, moments.pidTwo.pairs);
    sumChargedPairs.Fill(centrality, moments.chargedTwo.pairs);
    sumMixedQuadruplets.Fill(centrality, moments.mixedFour.pairs);
  }

  void write() {
    pidTwo.Write();
    chargedTwo.Write();
    mixedFour.Write();
    pidTwoByGroup.Write();
    chargedTwoByGroup.Write();
    mixedFourByGroup.Write();
    validEvents.Write();
    sumPidPairs.Write();
    sumChargedPairs.Write();
    sumMixedQuadruplets.Write();
  }
};

} // namespace pid_charged_sc
} // namespace ampt_analysis

#endif // AMPT_PID_CHARGED_SC_H

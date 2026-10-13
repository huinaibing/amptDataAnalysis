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
        eta == 0. || !isChargedPdg(track.pdgPid)) {
      continue;
    }
    const double absEta = std::abs(eta);
    const bool acceptedEta = settings.strictEtaBounds
                                 ? absEta > settings.innerAbsEta &&
                                       absEta < settings.etaMax
                                 : absEta >= settings.innerAbsEta &&
                                       absEta <= settings.etaMax;
    if (!acceptedEta) {
      continue;
    }
    const bool charged = config.chargedPt.contains(pt, settings.strictPtBounds);
    int mask = charged ? Mask::kRef : 0;
    const int absPdg = std::abs(track.pdgPid);
    const auto *species = findSpecies(absPdg);
    if (species && ptRangeForPdg(absPdg, config).contains(
                       pt, settings.strictPtBounds)) {
      mask |= species->mask;
      if (charged) {
        mask |= species->overlapMask;
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

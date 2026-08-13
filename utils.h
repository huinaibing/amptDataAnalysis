#ifndef AMPT_DATA_ANALYSIS_FLOW_CONTAINER_UTILS_H
#define AMPT_DATA_ANALYSIS_FLOW_CONTAINER_UTILS_H

#include "analysisUtils.h"

#include "PWGCF/GenericFramework/Core/FlowContainer.h"
#include "TProfile3D.h"

inline constexpr int kCompatibilityBootstrapSubsamples = 30;

namespace ampt_analysis {
inline void fillFlowProfile(FlowContainer &output, const char *profileName,
                            double centrality,
                            const CorrelationResult &correlation,
                            double random) {
  if (correlation.hasPairs()) {
    output.FillProfile(profileName, centrality, correlation.value(),
                       correlation.pairs, random);
  }
}

inline void fillFlowMeanPtProduct(FlowContainer &output,
                                  const char *profileName, double centrality,
                                  const CorrelationResult &correlation,
                                  const Event::PtMoments &moments,
                                  double random) {
  if (correlation.hasPairs() && moments.count > 0) {
    output.FillProfile(
        profileName, centrality, correlation.value() * moments.mean(),
        correlation.pairs * static_cast<double>(moments.count), random);
  }
}

inline void fillTrackWeightedFlow(FlowContainer &output,
                                  const char *profileName, double centrality,
                                  const CorrelationResult &correlation,
                                  const Event::PtMoments &moments,
                                  double random) {
  if (correlation.hasPairs() && moments.count > 0) {
    output.FillProfile(profileName, centrality, correlation.value(),
                       correlation.pairs * static_cast<double>(moments.count),
                       random);
  }
}

inline void fillMeanPtMoments(FlowContainer &output, double centrality,
                              const Event::PtMoments &moments, double random) {
  if (moments.count < 2) {
    return;
  }

  const double count = static_cast<double>(moments.count);
  const double orderedPairs = count * (count - 1.);
  output.FillProfile("hMeanPt", centrality, moments.mean(), count, random);
  output.FillProfile("ptAve", centrality, moments.mean(), orderedPairs, random);
  output.FillProfile("ptSquareAve", centrality, moments.distinctPairMean(),
                     orderedPairs, random);
}

/**
 * Fill the processData-style meanptCentNbs profiles from correlations that
 * were already calculated in the species loop.
 */
inline void fillV2PtRhoMeanPtProfiles(
    double centrality, double bootstrap, const Event::PtMoments &pidMoments,
    const CorrelationResult &refRef,
    const CorrelationResult &poiRef, const CorrelationResult &pure,
    TProfile3D &poiRefProfile, TProfile3D &refRefProfile,
    TProfile3D &pureProfile, TProfile3D &meanPtProfile) {
  if (pidMoments.count == 0 || !refRef.isPhysical() || !poiRef.hasPairs() ||
      poiRef.value() == 0.) {
    return;
  }

  const double meanPt = pidMoments.mean();
  const double count = static_cast<double>(pidMoments.count);
  poiRefProfile.Fill(meanPt, centrality, bootstrap, poiRef.value(),
                     count * poiRef.pairs);
  refRefProfile.Fill(meanPt, centrality, bootstrap, refRef.value(),
                     count * refRef.pairs);

  // POI-ref/ref-ref are valid independently of the POI-POI result.
  if (!pure.isPhysical()) {
    return;
  }
  pureProfile.Fill(meanPt, centrality, bootstrap, pure.value(), pure.pairs);
  meanPtProfile.Fill(meanPt, centrality, bootstrap, meanPt, count);
}
} // namespace ampt_analysis

// Compatibility wrappers for existing user macros outside this repository.
inline void CalculateAndFill(GFW *gfw, const CorrConfigManager &manager,
                             CorrType type, FlowContainer *output,
                             const char *profileName, double centrality,
                             float random) {
  const auto correlation =
      ampt_analysis::calculateCorrelation(*gfw, manager, type);
  ampt_analysis::fillFlowProfile(*output, profileName, centrality, correlation,
                                 random);
}

template <typename Predicate>
inline void CalculateCovV2ChargedPt(GFW *gfw, const CorrConfigManager &manager,
                                    CorrType type, FlowContainer *output,
                                    const Event &event, double centrality,
                                    const Predicate &accept, float random,
                                    const char *profileName = "covV2Pt") {
  const auto correlation =
      ampt_analysis::calculateCorrelation(*gfw, manager, type);
  const auto moments = event.GetPtMoments(accept);
  ampt_analysis::fillFlowMeanPtProduct(*output, profileName, centrality,
                                       correlation, moments, random);
}

template <typename Predicate>
inline void
CalculateC22TrackWeight(GFW *gfw, const CorrConfigManager &manager,
                        CorrType type, FlowContainer *output,
                        const Event &event, double centrality,
                        const Predicate &accept, float random,
                        const char *profileName = "c22TrackWeight") {
  const auto correlation =
      ampt_analysis::calculateCorrelation(*gfw, manager, type);
  const auto moments = event.GetPtMoments(accept);
  ampt_analysis::fillTrackWeightedFlow(*output, profileName, centrality,
                                       correlation, moments, random);
}

inline double getPidC22InOneEvent(GFW *gfw, const GFW::CorrConfig &configA,
                                  const GFW::CorrConfig &configB) {
  const auto a = ampt_analysis::calculateCorrelation(*gfw, configA);
  const auto b = ampt_analysis::calculateCorrelation(*gfw, configB);
  const ampt_analysis::CorrelationResult combined{a.numerator + b.numerator,
                                                  a.pairs + b.pairs};
  return combined.value();
}

inline void
FillMeanptCentBSProfile(const double &centrality, const double &random,
                        const double &pidMeanPt, const double &nPid, GFW *gfw,
                        const CorrConfigManager &manager, CorrType typeA,
                        CorrType typeB, CorrType pureType,
                        TProfile3D *poiRefProfile, TProfile3D *refRefProfile,
                        TProfile3D *pureProfile, TProfile3D *meanPtProfile) {
  if (!gfw || !poiRefProfile || !refRefProfile || !pureProfile ||
      !meanPtProfile || nPid <= 0.) {
    return;
  }

  Event::PtMoments pidMoments;
  pidMoments.count = static_cast<std::size_t>(nPid);
  pidMoments.sum = pidMeanPt * nPid;

  const auto refRef =
      ampt_analysis::calculateCorrelation(*gfw, manager, CorrType::Ref08Gap22);
  const auto poiRef =
      ampt_analysis::calculateCombinedCorrelation(*gfw, manager, typeA, typeB);
  const auto pure =
      ampt_analysis::calculateCorrelation(*gfw, manager, pureType);
  const auto bootstrapEdges = ampt_analysis::makeUniformEdges(
      kCompatibilityBootstrapSubsamples, 0.,
      static_cast<double>(kCompatibilityBootstrapSubsamples));
  const double bootstrap =
      ampt_analysis::sampleAxisCoordinate(bootstrapEdges, random);
  ampt_analysis::fillV2PtRhoMeanPtProfiles(
      centrality, bootstrap, pidMoments, refRef, poiRef, pure,
      *poiRefProfile, *refRefProfile, *pureProfile, *meanPtProfile);
}

#endif // AMPT_DATA_ANALYSIS_FLOW_CONTAINER_UTILS_H

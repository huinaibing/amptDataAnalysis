#include "../pidChargedSc.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace ampt_pid_charged_sc_test {

void near(double actual, double expected, const std::string &label) {
  if (!std::isfinite(actual) ||
      std::abs(actual - expected) > 1.e-10 * std::max(1., std::abs(expected))) {
    throw std::runtime_error(label + ": got " + std::to_string(actual) +
                             ", expected " + std::to_string(expected));
  }
}

Track track(int pdg, double pt, double eta, double phi) {
  Track result;
  result.pdgPid = pdg;
  result.p_x = pt * std::cos(phi);
  result.p_y = pt * std::sin(phi);
  result.p_z = pt * std::sinh(eta);
  return result;
}

// Independent explicit enumeration of particle tuples. This also tests that
// PID tracks outside the charged pT window remain in the PID sample, and that
// only the actual intersection is subtracted as overlap.
ampt_analysis::pid_charged_sc::Moments bruteForce(
    const Event &event, const ampt_analysis::AnalysisConfig &config,
    std::size_t species) {
  using namespace ampt_analysis;
  std::array<std::vector<std::size_t>, 2> pid, charged;
  const auto &settings = config.pidChargedScOutput;
  // species=3 replaces the PID leg by charged in the independent tuple
  // enumeration, yielding the pure charged four-track baseline.
  const auto *definition = species < 3 ? &speciesDefinitions().at(species) : nullptr;
  for (std::size_t i = 0; i < event.particles.size(); ++i) {
    const auto &particle = event.particles[i];
    const double eta = particle.GetEta();
    const double pt = particle.GetPt();
    if (!std::isfinite(eta) || !std::isfinite(pt) ||
        !std::isfinite(particle.GetPhi()) || eta == 0.) {
      continue;
    }
    const double absEta = std::abs(eta);
    if (settings.strictEtaBounds
            ? !(absEta > settings.innerAbsEta && absEta < settings.etaMax)
            : !(absEta >= settings.innerAbsEta && absEta <= settings.etaMax)) {
      continue;
    }
    const std::size_t side = eta < 0. ? 0 : 1;
    const bool selectedCharged = isChargedPdg(particle.pdgPid) &&
        config.chargedPt.contains(pt, settings.strictPtBounds);
    if (selectedCharged) {
      charged[side].push_back(i);
    }
    if (definition ? (std::abs(particle.pdgPid) == definition->absPdg &&
                      ptRangeForPdg(definition->absPdg, config).contains(
                          pt, settings.strictPtBounds)) : selectedCharged) {
      pid[side].push_back(i);
    }
  }
  pid_charged_sc::Moments result;
  const auto phi = [&](std::size_t i) { return event.particles[i].GetPhi(); };
  for (auto a : pid[0]) {
    for (auto b : pid[1]) {
      result.pidTwo.numerator += std::cos(2. * (phi(a) - phi(b)));
      ++result.pidTwo.pairs;
    }
  }
  for (auto a : charged[0]) {
    for (auto b : charged[1]) {
      result.chargedTwo.numerator += std::cos(2. * (phi(a) - phi(b)));
      ++result.chargedTwo.pairs;
    }
  }
  for (auto pa : pid[0]) {
    for (auto ca : charged[0]) {
      if (pa == ca) {
        continue;
      }
      for (auto pb : pid[1]) {
        for (auto cb : charged[1]) {
          if (pb == cb) {
            continue;
          }
          result.mixedFour.numerator +=
              std::cos(2. * (phi(pa) + phi(ca) - phi(pb) - phi(cb)));
          ++result.mixedFour.pairs;
        }
      }
    }
  }
  return result;
}

ampt_analysis::pid_charged_sc::ChargedMoments bruteForceFull(
    const Event &event, const ampt_analysis::AnalysisConfig &config,
    std::size_t species) {
  using namespace ampt_analysis;
  const auto &settings = config.pidChargedScOutput;
  std::vector<std::size_t> selected;
  for (std::size_t i = 0; i < event.particles.size(); ++i) {
    const auto &p = event.particles[i];
    const double pt = p.GetPt(), eta = p.GetEta();
    if (!std::isfinite(pt) || !std::isfinite(eta) ||
        !std::isfinite(p.GetPhi()) || !isChargedPdg(p.pdgPid) ||
        (settings.strictEtaBounds ? std::abs(eta) >= settings.etaMax
                                  : std::abs(eta) > settings.etaMax)) {
      continue;
    }
    const bool charged = species == 3;
    if ((!charged && std::abs(p.pdgPid) != speciesDefinitions().at(species).absPdg) ||
        !(charged ? config.chargedPt : ptRangeForPdg(std::abs(p.pdgPid), config))
             .contains(pt, settings.strictPtBounds)) {
      continue;
    }
    selected.push_back(i);
  }
  pid_charged_sc::ChargedMoments result;
  const auto phi = [&](std::size_t i) { return event.particles[i].GetPhi(); };
  for (auto a : selected) for (auto b : selected) {
    if (a == b) continue;
    result.two.numerator += std::cos(2. * (phi(a) - phi(b)));
    ++result.two.pairs;
    for (auto c : selected) for (auto d : selected) {
      if (c == a || c == b || d == a || d == b || c == d) continue;
      result.four.numerator +=
          std::cos(2. * (phi(a) + phi(b) - phi(c) - phi(d)));
      ++result.four.pairs;
    }
  }
  return result;
}

} // namespace ampt_pid_charged_sc_test

void test_pid_charged_sc(
    const char *analysisConfigFile = "../config/calculate_pid_charged_sc.json") {
  using namespace ampt_analysis;
  using namespace ampt_analysis::pid_charged_sc;
  using namespace ampt_pid_charged_sc_test;
  AnalysisConfig config = loadAnalysisConfig(
      analysisConfigFile, "pid_charged_sc_output");
  config.chargedPt = {0.4, 1.5};
  config.pionPt = {0.2, 2.};
  config.kaonPt = {0.3, 1.};
  config.protonPt = {0.6, 2.5};
  // Use reconstructed eta values as boundaries to exercise exact equality
  // with both open and closed selection, despite GFW's open region bounds.
  config.pidChargedScOutput.innerAbsEta = track(211, 1., 0.4, 0.).GetEta();
  config.pidChargedScOutput.etaMax = track(211, 1., 0.8, 0.).GetEta();
  GFW gfw;
  const PidChargedScCorrConfigManager correlations(
      gfw, config.pidChargedScOutput.innerAbsEta,
      config.pidChargedScOutput.etaMax);
  const std::vector<double> cent{0., 90.}, groups{0., 2., 5.};
  Profiles profiles(cent, groups);
  Profiles commonProfiles(cent, groups);
  ChargedProfiles chargedProfiles(cent, groups), commonChargedProfiles(cent, groups);
  std::array<std::unique_ptr<FlowVarianceProfiles>, 4> variance, commonVariance;
  for (std::size_t s = 0; s < variance.size(); ++s) {
    variance[s] = std::make_unique<FlowVarianceProfiles>(cent, groups);
    commonVariance[s] = std::make_unique<FlowVarianceProfiles>(cent, groups);
  }
  std::array<std::array<double, 3>, 4> varianceSums{}, varianceWeights{},
      commonVarianceSums{};
  std::array<double, 4> varianceFourWeights{};
  std::array<int, 4> varianceEvents{};
  std::array<double, 2> chargedSums{}, chargedWeights{}, commonChargedSums{};
  double chargedFourWeightSum = 0.;
  int chargedValidEvents = 0;
  std::array<double, 3> sums{}, weights{}, commonSums{};
  double fourWeightSum = 0.;
  int validEvents = 0;
  std::mt19937 engine(7231);
  std::uniform_real_distribution<double> pt(0.1, 2.7), phi(-3.14, 3.14);
  const std::array<int, 8> pdgs{211, -211, 321, -321, 2212, -2212, 11, 111};
  for (int sample = 0; sample < 120; ++sample) {
    config.pidChargedScOutput.strictPtBounds = sample % 2 == 0;
    config.pidChargedScOutput.strictEtaBounds = sample % 3 == 0;
    Event event;
    for (int i = 0; i < 18; ++i) {
      const double eta = (i % 2 ? 1. : -1.) * (i % 5 == 0 ? 0.1 : 0.6);
      event.particles.push_back(track(pdgs[(i + sample) % pdgs.size()],
                                      pt(engine), eta, phi(engine)));
    }
    // Deterministic pT endpoints, full-overlap low multiplicity, and eta=0.
    if (sample < 4) {
      event.particles = {track(211, 0.4, -0.6, 0.), track(-211, 1.5, -0.6, 0.),
                         track(211, 0.4, 0.6, 0.), track(-211, 1.5, 0.6, 0.),
                         track(211, 1., 0., 0.)};
    } else if (sample < 8) {
      event.particles = {track(211, 1., -0.4, 0.), track(-211, 1., -0.8, 0.),
                         track(211, 1., 0.4, 0.), track(-211, 1., 0.8, 0.),
                         track(211, 1., 0., 0.)};
    } else if (sample == 8) {
      event.particles.clear();
      for (int pdg : {211, 321, 2212}) {
        for (double eta : {-0.6, -0.1, 0., 0.6}) {
          event.particles.push_back(track(pdg, 0.8, eta, phi(engine)));
        }
      }
    }
    fillGfw(gfw, event, config);
    for (std::size_t s = 0; s < variance.size(); ++s) {
      const auto gap = bruteForce(event, config, s).pidTwo;
      const auto full = bruteForceFull(event, config, s);
      const auto actual = calculateFlowVarianceMoments(
          gfw, correlations,
          calculateCorrelation(gfw, s == 3 ? correlations.ChargedTwo()
                                           : correlations.PidTwo(s)), s);
      const FlowVarianceMoments expected{gap, full.two, full.four};
      const std::array<CorrelationResult, 3> a{
          actual.gapTwo, actual.fullTwo, actual.fullFour};
      const std::array<CorrelationResult, 3> e{gap, full.two, full.four};
      for (std::size_t i = 0; i < e.size(); ++i) {
        const std::string label = "variance sample " + std::to_string(sample) +
            " species " + std::to_string(s) + " moment " + std::to_string(i);
        near(a[i].pairs, e[i].pairs, label + " count");
        near(a[i].numerator, e[i].numerator, label + " numerator");
      }
      const double group = sampleAxisCoordinate(groups, (sample + 0.5) / 120.);
      config.pidChargedScOutput.flowVarianceEventWeight = "combinations";
      variance[s]->fill(10., group, actual, config.pidChargedScOutput);
      config.pidChargedScOutput.flowVarianceEventWeight = "full_quadruplets";
      commonVariance[s]->fill(10., group, actual, config.pidChargedScOutput);
      if (expected.valid()) {
        ++varianceEvents[s];
        varianceFourWeights[s] += full.four.pairs;
        for (std::size_t i = 0; i < e.size(); ++i) {
          varianceSums[s][i] += e[i].numerator;
          varianceWeights[s][i] += e[i].pairs;
          commonVarianceSums[s][i] += e[i].value() * full.four.pairs;
        }
      }
    }
    const auto chargedActual = calculateChargedMoments(gfw, correlations);
    const auto chargedExpected = bruteForce(event, config, 3);
    near(chargedActual.two.numerator, chargedExpected.chargedTwo.numerator,
         "charged pair numerator");
    near(chargedActual.two.pairs, chargedExpected.chargedTwo.pairs,
         "charged pair count");
    near(chargedActual.four.numerator, chargedExpected.mixedFour.numerator,
         "charged four numerator");
    near(chargedActual.four.pairs, chargedExpected.mixedFour.pairs,
         "charged four count");
    const double chargedGroup = sampleAxisCoordinate(groups, (sample + 0.5) / 120.);
    config.pidChargedScOutput.eventWeight = "combinations";
    chargedProfiles.fill(10., chargedGroup, chargedActual, config.pidChargedScOutput);
    config.pidChargedScOutput.eventWeight = "mixed_quadruplets";
    commonChargedProfiles.fill(10., chargedGroup, chargedActual, config.pidChargedScOutput);
    if (chargedExpected.valid()) {
      ++chargedValidEvents;
      chargedFourWeightSum += chargedExpected.mixedFour.pairs;
      const std::array<CorrelationResult, 2> e{
          chargedExpected.chargedTwo, chargedExpected.mixedFour};
      for (std::size_t i=0; i<e.size(); ++i) {
        chargedSums[i] += e[i].numerator;
        chargedWeights[i] += e[i].pairs;
        commonChargedSums[i] += e[i].value()*chargedExpected.mixedFour.pairs;
      }
    }
    for (std::size_t species = 0; species < 3; ++species) {
      const auto actual = calculateMoments(gfw, correlations, species);
      const auto expected = bruteForce(event, config, species);
      const std::array<CorrelationResult, 3> a{
          actual.pidTwo, actual.chargedTwo, actual.mixedFour};
      const std::array<CorrelationResult, 3> e{
          expected.pidTwo, expected.chargedTwo, expected.mixedFour};
      for (std::size_t i = 0; i < a.size(); ++i) {
        near(a[i].numerator, e[i].numerator, "tuple numerator");
        near(a[i].pairs, e[i].pairs, "tuple count");
      }
      if (species != 0) {
        continue;
      }
      const double group = sampleAxisCoordinate(groups, (sample + 0.5) / 120.);
      config.pidChargedScOutput.eventWeight = "combinations";
      profiles.fill(10., group, actual, config.pidChargedScOutput);
      config.pidChargedScOutput.eventWeight = "mixed_quadruplets";
      commonProfiles.fill(10., group, actual, config.pidChargedScOutput);
      if (!expected.valid()) {
        continue;
      }
      ++validEvents;
      fourWeightSum += expected.mixedFour.pairs;
      for (std::size_t i = 0; i < e.size(); ++i) {
        sums[i] += e[i].numerator;
        weights[i] += e[i].pairs;
        commonSums[i] += e[i].value() * expected.mixedFour.pairs;
      }
    }
  }
  if (validEvents < 2) {
    throw std::runtime_error("Insufficient valid test events");
  }
  const std::array<TProfile *, 3> p{
      &profiles.pidTwo, &profiles.chargedTwo, &profiles.mixedFour};
  const std::array<TProfile *, 3> cp{
      &commonProfiles.pidTwo, &commonProfiles.chargedTwo, &commonProfiles.mixedFour};
  const std::array<TProfile2D *, 3> pg{
      &profiles.pidTwoByGroup, &profiles.chargedTwoByGroup, &profiles.mixedFourByGroup};
  for (std::size_t i = 0; i < p.size(); ++i) {
    near(p[i]->GetBinEntries(1), weights[i], "native profile weight");
    near(p[i]->GetBinContent(1), sums[i] / weights[i], "native profile mean");
    near(cp[i]->GetBinEntries(1), fourWeightSum, "common profile weight");
    near(cp[i]->GetBinContent(1), commonSums[i] / fourWeightSum,
         "common profile mean");
    double groupWeight = 0., groupNumerator = 0.;
    for (int group = 1; group <= 2; ++group) {
      const int bin = pg[i]->GetBin(1, group);
      groupWeight += pg[i]->GetBinEntries(bin);
      groupNumerator += pg[i]->GetBinEntries(bin) * pg[i]->GetBinContent(bin);
    }
    near(groupWeight, weights[i], "subsample weights");
    near(groupNumerator, sums[i], "subsample numerators");
  }
  near(profiles.validEvents.GetBinContent(1), validEvents, "valid event count");
  const std::array<TProfile *, 2> chargedP{&chargedProfiles.two, &chargedProfiles.four};
  const std::array<TProfile *, 2> chargedCP{
      &commonChargedProfiles.two, &commonChargedProfiles.four};
  const std::array<TProfile2D *, 2> chargedPG{
      &chargedProfiles.twoByGroup, &chargedProfiles.fourByGroup};
  for (std::size_t i=0; i<chargedP.size(); ++i) {
    near(chargedP[i]->GetBinEntries(1), chargedWeights[i], "charged native weight");
    near(chargedP[i]->GetBinContent(1), chargedSums[i]/chargedWeights[i], "charged native mean");
    near(chargedCP[i]->GetBinEntries(1), chargedFourWeightSum, "charged common weight");
    near(chargedCP[i]->GetBinContent(1), commonChargedSums[i]/chargedFourWeightSum,
         "charged common mean");
    double w=0., sum=0.;
    for (int g=1; g<=2; ++g) {
      const int b=chargedPG[i]->GetBin(1,g);
      w+=chargedPG[i]->GetBinEntries(b);
      sum+=chargedPG[i]->GetBinEntries(b)*chargedPG[i]->GetBinContent(b);
    }
    near(w, chargedWeights[i], "charged group weights");
    near(sum, chargedSums[i], "charged group numerator");
  }
  near(chargedProfiles.validEvents.GetBinContent(1), chargedValidEvents, "charged valid events");
  near(chargedProfiles.sumPairs.GetBinContent(1), chargedWeights[0], "charged pair sum");
  near(chargedProfiles.sumQuadruplets.GetBinContent(1), chargedWeights[1], "charged four sum");
  for (std::size_t s = 0; s < variance.size(); ++s) {
    const std::array<TProfile *, 3> p{
        &variance[s]->gapTwo, &variance[s]->fullTwo, &variance[s]->fullFour};
    const std::array<TProfile *, 3> cp{
        &commonVariance[s]->gapTwo, &commonVariance[s]->fullTwo,
        &commonVariance[s]->fullFour};
    const std::array<TProfile2D *, 3> pg{
        &variance[s]->gapTwoByGroup, &variance[s]->fullTwoByGroup,
        &variance[s]->fullFourByGroup};
    if (varianceEvents[s] == 0) throw std::runtime_error("No variance test events");
    near(variance[s]->validEvents.GetBinContent(1), varianceEvents[s],
         "variance valid events");
    for (std::size_t i = 0; i < p.size(); ++i) {
      near(p[i]->GetBinEntries(1), varianceWeights[s][i], "variance weight");
      near(p[i]->GetBinContent(1), varianceSums[s][i] / varianceWeights[s][i],
           "variance mean");
      near(cp[i]->GetBinEntries(1), varianceFourWeights[s], "variance common weight");
      near(cp[i]->GetBinContent(1), commonVarianceSums[s][i] / varianceFourWeights[s],
           "variance common mean");
      double sum = 0., weight = 0.;
      for (int g = 1; g <= 2; ++g) {
        const int bin = pg[i]->GetBin(1, g);
        weight += pg[i]->GetBinEntries(bin);
        sum += pg[i]->GetBinEntries(bin) * pg[i]->GetBinContent(bin);
      }
      near(weight, varianceWeights[s][i], "variance group weights");
      near(sum, varianceSums[s][i], "variance group numerator");
    }
  }
  // Full-only tracks must be included even at eta=0; their absence from the
  // gap makes the matched variance triplet invalid, not the full correlator.
  Event centralOnly;
  for (double eta : {-0.2, 0., 0.1, 0.2})
    centralOnly.particles.push_back(track(211, 1., eta, 0.));
  fillGfw(gfw, centralOnly, config);
  near(calculateCorrelation(gfw, correlations.FullTwo(0)).pairs, 12., "central full pairs");
  near(calculateCorrelation(gfw, correlations.FullFour(0)).pairs, 24., "central full fours");
  near(calculateMoments(gfw, correlations, 0).pidTwo.pairs, 0., "central gap empty");
  Event lowPid;
  for (double eta : {-0.6, 0.6}) {
    lowPid.particles.push_back(track(211, 1., eta, 0.));
    lowPid.particles.push_back(track(11, 1., eta, 0.));
  }
  fillGfw(gfw, lowPid, config);
  const auto lowPidSc = calculateMoments(gfw, correlations, 0);
  if (!lowPidSc.valid() ||
      calculateFlowVarianceMoments(gfw, correlations, lowPidSc.pidTwo, 0).valid())
    throw std::runtime_error("Four-POI variance cut changed the SC sample");
  // Charged baseline must remain valid when pion/kaon/proton are absent.
  Event electrons;
  for (double eta : {-0.6,0.6}) for (double angle : {0.,0.3})
    electrons.particles.push_back(track(11,1.,eta,angle));
  fillGfw(gfw, electrons, config);
  if (!calculateChargedMoments(gfw, correlations).valid() ||
      calculateMoments(gfw, correlations, 0).valid())
    throw std::runtime_error("Charged baseline unexpectedly depends on PID");
  fillGfw(gfw, Event{}, config);
  const auto emptyCharged = calculateChargedMoments(gfw, correlations);
  near(emptyCharged.two.pairs, 0., "empty charged pair");
  near(emptyCharged.four.pairs, 0., "empty charged four");
  if (emptyCharged.valid()) throw std::runtime_error("Empty charged baseline valid");
  for (std::size_t species = 0; species < 3; ++species) {
    const auto empty = calculateMoments(gfw, correlations, species);
    near(empty.pidTwo.pairs, 0., "empty event clears PID");
    near(empty.chargedTwo.pairs, 0., "empty event clears charged");
    near(empty.mixedFour.pairs, 0., "empty event clears mixed moment");
    if (empty.valid()) {
      throw std::runtime_error("Empty event unexpectedly valid");
    }
  }
  for (std::size_t s = 0; s < variance.size(); ++s) {
    near(calculateCorrelation(gfw, correlations.FullTwo(s)).pairs, 0., "empty full pairs");
    near(calculateCorrelation(gfw, correlations.FullFour(s)).pairs, 0., "empty full fours");
  }
  std::cout << "Validated GFW against distinct PID/charged tuples for 120 events, "
               "partial overlap, pure charged baseline, four-species full-eta variance "
               "components, multiplicity weight modes and subsamples.\n";
}

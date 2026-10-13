#include "../pidChargedSc.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
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
  const auto &definition = speciesDefinitions().at(species);
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
    if (isChargedPdg(particle.pdgPid) &&
        config.chargedPt.contains(pt, settings.strictPtBounds)) {
      charged[side].push_back(i);
    }
    if (std::abs(particle.pdgPid) == definition.absPdg &&
        ptRangeForPdg(definition.absPdg, config).contains(
            pt, settings.strictPtBounds)) {
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

} // namespace ampt_pid_charged_sc_test

void test_pid_charged_sc(
    const char *analysisConfigFile = "../config/config.json") {
  using namespace ampt_analysis;
  using namespace ampt_analysis::pid_charged_sc;
  using namespace ampt_pid_charged_sc_test;
  AnalysisConfig config = loadAnalysisConfig(analysisConfigFile);
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
    }
    fillGfw(gfw, event, config);
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
  fillGfw(gfw, Event{}, config);
  for (std::size_t species = 0; species < 3; ++species) {
    const auto empty = calculateMoments(gfw, correlations, species);
    near(empty.pidTwo.pairs, 0., "empty event clears PID");
    near(empty.chargedTwo.pairs, 0., "empty event clears charged");
    near(empty.mixedFour.pairs, 0., "empty event clears mixed moment");
    if (empty.valid()) {
      throw std::runtime_error("Empty event unexpectedly valid");
    }
  }
  std::cout << "Validated GFW against distinct PID/charged tuples for 120 events, "
               "partial overlap, both multiplicity weight modes and subsamples.\n";
}

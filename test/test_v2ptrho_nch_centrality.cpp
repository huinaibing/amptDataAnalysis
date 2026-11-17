#include "../analysisUtils.h"

#include <fstream>
#include <iostream>
#include <stdexcept>

namespace ampt_nch_centrality_test {
void require(bool condition, const char *message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

Track track(int pdg, double pt, double eta = 0.) {
  Track result;
  result.pdgPid = pdg;
  result.p_x = pt;
  result.p_z = pt * std::sinh(eta);
  return result;
}
} // namespace ampt_nch_centrality_test

void test_v2ptrho_nch_centrality(
    const char *configFile = "../config/calculate_v2ptrho.json") {
  using namespace ampt_analysis;
  using namespace ampt_nch_centrality_test;
  const auto config = loadAnalysisConfig(configFile, "v2_pt_rho_output");
  const auto &calibration = config.v2PtRhoNchCentrality;
  require(calibration.enabled, "The supplied config must enable the fitted curve");

  const int anchors[]{2321, 1884, 1250, 800, 482, 268, 134, 59, 23, 8};
  const double percentiles[]{5, 10, 20, 30, 40, 50, 60, 70, 80, 90};
  for (int i = 0; i < 10; ++i) {
    require(std::abs(centralityFromNch(anchors[i], calibration) - percentiles[i]) < 0.3,
            "Fit no longer agrees with supplied percentile anchors");
    require(centralityFromNch(anchors[i] + 1, calibration) <
                centralityFromNch(anchors[i], calibration),
            "Curve must vary within old threshold classes");
  }
  double previous = 100.;
  for (int n = 0; n <= 10000; ++n) {
    const double c = centralityFromNch(n, calibration);
    require(std::isfinite(c) && c >= 0. && c <= 100. && c <= previous,
            "Curve must be finite, clipped and decreasing");
    previous = c;
  }
  require(centralityFromNch(0, calibration) == 100., "Peripheral tail clipping");
  require(centralityFromNch(10000, calibration) == 0., "Central tail clipping");
  require(axisBinCount(config.v2PtRhoOutput.centralityAxis) == 90,
          "Analysis must retain 1% bins over 0-90%");
  for (int pdg : {211, -211, 321, -321, 2212, -2212, 213, 3334, 2}) {
    require(acceptsNchCentralityTrack(track(pdg, 1.), calibration),
            "Xinli charged-table entry rejected");
  }
  for (int pdg : {0, 111, 2112, -2112, 11, -11, 13, 999999}) {
    require(!acceptsNchCentralityTrack(track(pdg, 1.), calibration),
            "Particle absent from Xinli charge table accepted");
  }
  for (double pt : {0., 0.2, 3., 4.}) {
    require(!acceptsNchCentralityTrack(track(211, pt), calibration),
            "Strict pT boundary violated");
  }
  for (double eta : {-0.8001, 0.8001}) {
    require(!acceptsNchCentralityTrack(track(211, 1., eta), calibration),
            "Eta acceptance violated");
  }
  auto boundary = calibration;
  const auto atEtaEdge = track(211, 1., 0.8);
  boundary.absEtaMax = std::abs(atEtaEdge.GetEta());
  require(!acceptsNchCentralityTrack(atEtaEdge, boundary),
          "Strict eta equality boundary violated");

  std::ifstream input(configFile);
  nlohmann::json document;
  input >> document;
  auto legacy = document.at("v2_pt_rho_output");
  legacy.erase("nch_centrality");
  legacy["axes"].erase("nch_centrality");
  legacy["qa_axes"].erase("nch_centrality");
  require(!readNchCentralityConfig(legacy).enabled,
          "Old configs must keep the new mode disabled");

  auto custom = document.at("v2_pt_rho_output");
  custom["nch_centrality"]["curve"]["coefficients"][0] = 116.;
  const auto customConfig = readNchCentralityConfig(custom);
  require(customConfig.enabled && customConfig.coefficients[0] == 116.,
          "Custom fit coefficients ignored");
  const std::vector<std::vector<double>> invalidCoefficients = {
      {100, -1, 0}, {100, -1, 0, 1}, {100, 1, 0, -1}};
  for (const auto &p : invalidCoefficients) {
    auto invalid = custom;
    invalid["nch_centrality"]["curve"]["coefficients"] = p;
    bool rejected = false;
    try {
      readNchCentralityConfig(invalid);
    } catch (const std::runtime_error &) {
      rejected = true;
    }
    require(rejected, "Invalid or nonmonotonic fit accepted");
  }
  auto thresholdOnly = custom;
  thresholdOnly["nch_centrality"].erase("curve");
  thresholdOnly["nch_centrality"]["thresholds"] = {2321, 1884};
  bool rejected = false;
  try {
    readNchCentralityConfig(thresholdOnly);
  } catch (const std::runtime_error &) {
    rejected = true;
  }
  require(rejected, "Enabled threshold-only mode must not silently replace the fit");
  std::cout << "PASS: curve anchors, monotonicity, clipping, fine bins, strict "
               "acceptance, legacy configs and invalid-fit rejection.\n";
}

#ifndef AMPT_DATA_ANALYSIS_CONFIG_H
#define AMPT_DATA_ANALYSIS_CONFIG_H

#include "corrConfigManager.h"
#include "dataFrame/json.hpp"

#include "TPDGCode.h"

#include <algorithm>
#include <array>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ampt_analysis {
struct PtRange {
  double min = 0.2;
  double max = 3.0;

  bool contains(double pt, bool strict) const {
    return strict ? (pt > min && pt < max) : (pt >= min && pt <= max);
  }
};

enum class AxisBinning { Uniform, Variable };

// Every histogram axis supports uniform bins or explicit variable-width bins.
struct AxisConfig {
  AxisBinning binning = AxisBinning::Uniform;
  int bins = 0;
  double min = 0.;
  double max = 0.;
  std::vector<double> edges;
};

struct OutputConfig {
  bool strictPtBounds = true;
  bool strictEtaBounds = true;
  AxisConfig centralityAxis;
  AxisConfig meanPtAxis;
  AxisConfig bootstrapAxis;
};

struct AnalysisConfig {
  // Flow particles cover [-flowEtaMax, -flowEtaGap] and
  // [flowEtaGap, flowEtaMax]. Mean-pT particles use an independent interval.
  double flowEtaMax = 0.8;
  double flowEtaGap = 0.4;
  double meanPtEtaMax = 0.4;

  PtRange chargedPt;
  PtRange pionPt;
  PtRange kaonPt;
  PtRange protonPt;

  unsigned int randomSeed = 0;
  double centralityPi = 3.14159265358979323846;
  double centralityReferenceRadiusSquared = 7.84;

  OutputConfig v2PtRhoOutput{
      false,
      false,
      {AxisBinning::Uniform, 90, 0., 90., {}},
      {AxisBinning::Uniform, 1000, 0., 3., {}},
      {AxisBinning::Uniform, 30, 0., 30., {}}};
  OutputConfig c22DeltaPtOutput{
      true,
      true,
      {AxisBinning::Variable,
       0,
       0.,
       0.,
       {0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.,  10., 12., 14.,
        16., 18., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70.,
        80., 90.}},
      {AxisBinning::Uniform, 300, 0., 3., {}},
      {AxisBinning::Uniform, 30, 0., 30., {}}};
  bool c22UsePure = false;

  // These working values are set from the selected output mode by each macro.
  bool strictPtBounds = true;
  bool strictEtaBounds = true;
};

inline const AnalysisConfig &defaultConfig() {
  static const AnalysisConfig config;
  return config;
}

inline std::vector<double> makeUniformEdges(int nBins, double min, double max) {
  if (nBins <= 0 || max <= min) {
    throw std::invalid_argument(
        "Uniform axis requires bins > 0 and max > min");
  }
  std::vector<double> edges(static_cast<std::size_t>(nBins) + 1);
  const double width = (max - min) / static_cast<double>(nBins);
  for (int i = 0; i <= nBins; ++i) {
    edges[static_cast<std::size_t>(i)] = min + static_cast<double>(i) * width;
  }
  return edges;
}

inline std::vector<double> makeAxisEdges(const AxisConfig &axis) {
  if (axis.binning == AxisBinning::Variable) {
    if (axis.edges.size() < 2 ||
        !std::is_sorted(axis.edges.begin(), axis.edges.end()) ||
        std::adjacent_find(axis.edges.begin(), axis.edges.end()) !=
            axis.edges.end()) {
      throw std::invalid_argument(
          "Explicit axis edges must contain at least two increasing values");
    }
    return axis.edges;
  }
  return makeUniformEdges(axis.bins, axis.min, axis.max);
}

inline int axisBinCount(const AxisConfig &axis) {
  return static_cast<int>(makeAxisEdges(axis).size()) - 1;
}

// Map a uniform random number to an axis with equal probability per bin.
// This also works for variable-width bootstrap axes without width bias.
inline double sampleAxisCoordinate(const std::vector<double> &edges,
                                   double random) {
  if (edges.size() < 2 || random < 0. || random >= 1.) {
    throw std::invalid_argument("Cannot sample invalid axis or random value");
  }
  const std::size_t nBins = edges.size() - 1;
  const double scaled = random * static_cast<double>(nBins);
  const std::size_t bin =
      std::min(static_cast<std::size_t>(scaled), nBins - 1);
  const double fraction = scaled - static_cast<double>(bin);
  return edges[bin] + fraction * (edges[bin + 1] - edges[bin]);
}

inline PtRange readPtRange(const nlohmann::json &document,
                           const std::string &name) {
  PtRange range{document.at("min").get<double>(),
                document.at("max").get<double>()};
  if (range.max <= range.min) {
    throw std::runtime_error("Invalid pT range for " + name +
                             ": max must be greater than min");
  }
  return range;
}

inline AxisConfig readAxisConfig(const nlohmann::json &document,
                                 const std::string &name) {
  AxisConfig axis;
  const std::string binning = document.at("binning").get<std::string>();
  if (binning == "variable") {
    axis.binning = AxisBinning::Variable;
    axis.edges = document.at("edges").get<std::vector<double>>();
  } else if (binning == "uniform") {
    axis.binning = AxisBinning::Uniform;
    axis.bins = document.at("bins").get<int>();
    axis.min = document.at("min").get<double>();
    axis.max = document.at("max").get<double>();
  } else {
    throw std::runtime_error("Invalid axis '" + name +
                             "': binning must be 'uniform' or 'variable'");
  }
  try {
    makeAxisEdges(axis);
  } catch (const std::exception &error) {
    throw std::runtime_error("Invalid axis '" + name + "': " + error.what());
  }
  return axis;
}

inline OutputConfig readOutputConfig(const nlohmann::json &document,
                                     const std::string &name) {
  const auto &axes = document.at("axes");
  return {document.at("strict_pt_bounds").get<bool>(),
          document.at("strict_eta_bounds").get<bool>(),
          readAxisConfig(axes.at("centrality"), name + ".centrality"),
          readAxisConfig(axes.at("mean_pt"), name + ".mean_pt"),
          readAxisConfig(axes.at("bootstrap"), name + ".bootstrap")};
}

inline AnalysisConfig loadAnalysisConfig(const std::string &jsonPath) {
  std::ifstream input(jsonPath);
  if (!input.is_open()) {
    throw std::runtime_error("Cannot open analysis configuration: " +
                             jsonPath);
  }

  nlohmann::json document;
  input >> document;

  AnalysisConfig config;
  const auto &selection = document.at("selection");
  const auto &flowEta = selection.at("flow_subevents");
  config.flowEtaMax = flowEta.at("eta_max").get<double>();
  config.flowEtaGap = flowEta.at("inner_abs_eta").get<double>();
  config.meanPtEtaMax =
      selection.at("mean_pt_abs_eta_max").get<double>();
  const auto &ptRanges = selection.at("pt_ranges");
  config.chargedPt = readPtRange(ptRanges.at("charged"), "charged");
  config.pionPt = readPtRange(ptRanges.at("pion"), "pion");
  config.kaonPt = readPtRange(ptRanges.at("kaon"), "kaon");
  config.protonPt = readPtRange(ptRanges.at("proton"), "proton");

  const auto &centrality = document.at("centrality_from_impact_parameter");
  config.centralityPi = centrality.at("pi").get<double>();
  config.centralityReferenceRadiusSquared =
      centrality.at("reference_radius_squared").get<double>();

  const auto &bootstrap = document.at("bootstrap");
  config.randomSeed = bootstrap.at("random_seed").get<unsigned int>();

  config.v2PtRhoOutput =
      readOutputConfig(document.at("v2_pt_rho_output"), "v2_pt_rho_output");
  const auto &c22Output = document.at("c22_delta_pt_output");
  config.c22DeltaPtOutput =
      readOutputConfig(c22Output, "c22_delta_pt_output");
  config.c22UsePure = c22Output.at("use_pure").get<bool>();

  if (config.flowEtaGap < 0. || config.flowEtaMax <= config.flowEtaGap) {
    throw std::runtime_error(
        "selection.flow_subevents requires 0 <= inner_abs_eta < eta_max");
  }
  if (config.meanPtEtaMax <= 0.) {
    throw std::runtime_error(
        "selection.mean_pt_abs_eta_max must be positive");
  }
  if (config.centralityPi <= 0. ||
      config.centralityReferenceRadiusSquared <= 0.) {
    throw std::runtime_error(
        "Centrality conversion constants must be positive");
  }
  return config;
}

enum class Species { Pion, Kaon, Proton };

struct SpeciesDefinition {
  Species species;
  int absPdg;
  int mask;
  int overlapMask;
  CorrType poiRefA;
  CorrType poiRefB;
  CorrType pure;
  CorrType fullTwoA;
  CorrType fullTwoB;
  CorrType fourA;
  CorrType fourB;
  const char *flowContainerName;
  const char *speciesName;
  const char *axisLabel;
};

inline const std::array<SpeciesDefinition, 3> &speciesDefinitions() {
  static const std::array<SpeciesDefinition, 3> definitions{
      {{Species::Pion, PDG_t::kPiPlus, Mask::kPion, Mask::kPionOverlap,
        CorrType::Pion08Gap22a, CorrType::Pion08Gap22b, CorrType::PiPi08Gap22,
        CorrType::Pion0Gap22a_Full, CorrType::Pion0Gap22b_Full,
        CorrType::Pion0Gap24a, CorrType::Pion0Gap24b, "FlowContainerPi", "Pion",
        "#pi"},
       {Species::Kaon, PDG_t::kKPlus, Mask::kKaon, Mask::kKaonOverlap,
        CorrType::Kaon08Gap22a, CorrType::Kaon08Gap22b, CorrType::KaKa08Gap22,
        CorrType::Kaon0Gap22a_Full, CorrType::Kaon0Gap22b_Full,
        CorrType::Kaon0Gap24a, CorrType::Kaon0Gap24b, "FlowContainerKa", "Kaon",
        "K"},
       {Species::Proton, PDG_t::kProton, Mask::kProton, Mask::kProtonOverlap,
        CorrType::Prot08Gap22a, CorrType::Prot08Gap22b, CorrType::PrPr08Gap22,
        CorrType::Prot0Gap22a_Full, CorrType::Prot0Gap22b_Full,
        CorrType::Prot0Gap24a, CorrType::Prot0Gap24b, "FlowContainerPr",
        "Proton", "p"}}};
  return definitions;
}

inline const PtRange &ptRangeForPdg(int absPdg, const AnalysisConfig &config) {
  switch (absPdg) {
  case PDG_t::kPiPlus:
    return config.pionPt;
  case PDG_t::kKPlus:
    return config.kaonPt;
  case PDG_t::kProton:
    return config.protonPt;
  default:
    return config.chargedPt;
  }
}

inline const SpeciesDefinition *findSpecies(int absPdg) {
  for (const auto &definition : speciesDefinitions()) {
    if (definition.absPdg == absPdg) {
      return &definition;
    }
  }
  return nullptr;
}
} // namespace ampt_analysis

#endif // AMPT_DATA_ANALYSIS_CONFIG_H

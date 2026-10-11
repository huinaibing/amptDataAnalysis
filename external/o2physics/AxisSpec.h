// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// Minimal ROOT-only subset of O2 Framework/HistogramSpec.h (AxisSpec).
// See README.md in this directory for source revision and licensing.
#ifndef AMPT_STANDALONE_AXIS_SPEC_H
#define AMPT_STANDALONE_AXIS_SPEC_H

#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace o2::framework {
struct AxisSpec {
  AxisSpec(std::vector<double> edges,
           std::optional<std::string> axisTitle = std::nullopt,
           std::optional<std::string> axisName = std::nullopt)
      : binEdges(std::move(edges)), title(std::move(axisTitle)),
        name(std::move(axisName)) {}

  AxisSpec(int bins, double minimum, double maximum,
           std::optional<std::string> axisTitle = std::nullopt,
           std::optional<std::string> axisName = std::nullopt)
      : nBins(bins), binEdges{minimum, maximum}, title(std::move(axisTitle)),
        name(std::move(axisName)) {
    if (minimum > maximum) {
      throw std::invalid_argument("Defined ill-defined axis");
    }
  }

  std::optional<int> nBins{};
  std::vector<double> binEdges{};
  std::optional<std::string> title{};
  std::optional<std::string> name{};
};
} // namespace o2::framework

#endif

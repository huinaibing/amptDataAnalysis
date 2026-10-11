#include "../analysisConfig.h"
#include "../analysisUtils.h"

#include "TFile.h"
#include "TH1D.h"
#include "TNamed.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ampt_pid_v2_rp_gap_variance_macro {

struct Subevent {
  int count = 0;
  double sumCos2Phi = 0.;
  void add(double phi) {
    ++count;
    sumCos2Phi += std::cos(2. * phi); // fixed AMPT reaction plane = 0
  }
  double v2() const { return sumCos2Phi / count; }
};

struct SpeciesMoments {
  Subevent negative;
  Subevent positive;
};

// The two first moments use the multiplicity of their own eta subevent as the
// event weight. The cross moment uses the product of the two subevent
// multiplicities. TProfile/TProfile2D keep mergeable weighted sums for ROOT
// hadd.
struct Profiles {
  TProfile negative;
  TProfile positive;
  TProfile product;
  TH1D validEvents;
  TProfile2D negativeByGroup;
  TProfile2D positiveByGroup;
  TProfile2D productByGroup;

  Profiles(const std::string &name, const std::vector<double> &edges, int nGroups)
      : negative(("pV2RPMinus" + name).c_str(),
                 ";Centrality (%);<v_{2,-}^{RP}>",
                 static_cast<int>(edges.size()) - 1, edges.data()),
        positive(("pV2RPPlus" + name).c_str(),
                 ";Centrality (%);<v_{2,+}^{RP}>",
                 static_cast<int>(edges.size()) - 1, edges.data()),
        product(("pV2RPMinusTimesPlus" + name).c_str(),
                ";Centrality (%);<v_{2,-}^{RP} v_{2,+}^{RP}>",
                static_cast<int>(edges.size()) - 1, edges.data()),
        validEvents(("hValidEvents" + name).c_str(),
                    ";Centrality (%);Events with both PID subevents",
                    static_cast<int>(edges.size()) - 1, edges.data()),
        negativeByGroup(("pV2RPMinusByGroup" + name).c_str(),
                        ";Centrality (%);Subsample;v_{2,-}^{RP}",
                        static_cast<int>(edges.size()) - 1, edges.data(),
                        nGroups, 0., static_cast<double>(nGroups)),
        positiveByGroup(("pV2RPPlusByGroup" + name).c_str(),
                        ";Centrality (%);Subsample;v_{2,+}^{RP}",
                        static_cast<int>(edges.size()) - 1, edges.data(),
                        nGroups, 0., static_cast<double>(nGroups)),
        productByGroup(("pV2RPMinusTimesPlusByGroup" + name).c_str(),
                       ";Centrality (%);Subsample;v_{2,-}^{RP}v_{2,+}^{RP}",
                       static_cast<int>(edges.size()) - 1, edges.data(),
                       nGroups, 0., static_cast<double>(nGroups)) {}

  void fill(double centrality, int group, double a, double b,
            int countNegative, int countPositive) {
    const double negativeWeight = static_cast<double>(countNegative);
    const double positiveWeight = static_cast<double>(countPositive);
    const double productWeight = negativeWeight * positiveWeight;
    negative.Fill(centrality, a, negativeWeight);
    positive.Fill(centrality, b, positiveWeight);
    product.Fill(centrality, a * b, productWeight);
    validEvents.Fill(centrality);
    const double groupCoordinate = group + 0.5;
    negativeByGroup.Fill(centrality, groupCoordinate, a, negativeWeight);
    positiveByGroup.Fill(centrality, groupCoordinate, b, positiveWeight);
    productByGroup.Fill(centrality, groupCoordinate, a * b, productWeight);
  }

  void write() {
    negative.Write();
    positive.Write();
    product.Write();
    validEvents.Write();
    negativeByGroup.Write();
    positiveByGroup.Write();
    productByGroup.Write();
  }
};

} // namespace ampt_pid_v2_rp_gap_variance_macro

/**
 * Save mergeable, multiplicity-weighted RP-projected PID v2 moments for two
 * eta subevents. The A and B first moments use N_A and N_B, respectively;
 * the cross moment uses N_A*N_B. Both sides must have at least one selected
 * PID particle.
 * No variance (or its uncertainty) is constructed in this upstream macro.
 *
 * Run directly in a ROOT environment from script/, or use
 * run_parallel.sh calculate_pid_v2_rp_gap_variance <workers> ...
 */
void calculate_pid_v2_rp_gap_variance(
    const char *inputConfigFile = "../config/cent_cfg_ampt.json",
    const char *outputFile = "pid_v2_rp_gap_moments.root",
    int maxFilesPerConfig = -1, int maxConfigs = -1,
    const char *analysisConfigFile = "../config/config.json",
    int shardIndex = 0, int shardCount = 1)
{
  using namespace ampt_analysis;
  using namespace ampt_pid_v2_rp_gap_variance_macro;
  const AnalysisConfig config = loadAnalysisConfig(analysisConfigFile);
  const auto &settings = config.pidV2RPGapVarianceOutput;
  const auto edges = makeAxisEdges(settings.centralityAxis);
  const int nGroups = axisBinCount(settings.bootstrapAxis);
  const std::array<const char *, 3> names{"Pion", "Kaon", "Proton"};
  std::array<std::unique_ptr<Profiles>, 3> profiles;
  for (int species = 0; species < 3; ++species) {
    profiles[species] = std::make_unique<Profiles>(names[species], edges, nGroups);
  }

  const Long64_t processedEvents = forEachConfiguredEvent(
      inputConfigFile, maxFilesPerConfig, maxConfigs, shardIndex, shardCount,
      [&](const Event &event, const CentralityConfig &,
          const EventIdentity &identity) {
        const double centrality = centralityFromImpactParameter(event.imp, config);
        if (!std::isfinite(centrality) || centrality < edges.front() ||
            centrality >= edges.back()) {
          return;
        }
        std::array<SpeciesMoments, 3> eventMoments{};
        for (const auto &track : event.particles) {
          if (!isChargedPdg(track.pdgPid)) {
            continue;
          }
          const int absPdg = std::abs(track.pdgPid);
          const SpeciesDefinition *species = findSpecies(absPdg);
          if (!species || !ptRangeForPdg(absPdg, config).contains(
                              track.GetPt(), settings.strictPtBounds)) {
            continue;
          }
          const double eta = track.GetEta();
          const double absEta = std::abs(eta);
          const bool accepted = settings.strictEtaBounds
                                    ? absEta > settings.innerAbsEta &&
                                          absEta < settings.etaMax
                                    : absEta >= settings.innerAbsEta &&
                                          absEta <= settings.etaMax;
          if (!accepted || eta == 0.) {
            continue;
          }
          auto &moments = eventMoments[static_cast<int>(species->species)];
          (eta < 0. ? moments.negative : moments.positive).add(track.GetPhi());
        }
        const double random = bootstrapRandomValue(config.randomSeed, identity);
        const int group = std::min(static_cast<int>(random * nGroups), nGroups - 1);
        for (int species = 0; species < 3; ++species) {
          const auto &sample = eventMoments[species];
          if (sample.negative.count == 0 || sample.positive.count == 0) {
            continue;
          }
          profiles[species]->fill(centrality, group,
                                  sample.negative.v2(), sample.positive.v2(),
                                  sample.negative.count, sample.positive.count);
        }
      });

  TFile output(outputFile, "RECREATE");
  if (output.IsZombie()) {
    throw std::runtime_error(std::string("Cannot create output file: ") + outputFile);
  }
  auto *directory = output.mkdir("pid-v2-rp-gap-moments");
  directory->cd();
  for (auto &species : profiles) {
    species->write();
  }
  // hadd sums TParameter values; TNamed preserves these fixed settings.
  const auto writeSetting = [](const char *name, double value) {
    std::ostringstream text;
    text << std::setprecision(17) << value;
    TNamed(name, text.str().c_str()).Write();
  };
  writeSetting("reactionPlaneAngle", 0.);
  writeSetting("innerAbsEta", settings.innerAbsEta);
  writeSetting("etaMax", settings.etaMax);
  TNamed("subsampleGroups", std::to_string(nGroups).c_str()).Write();
  output.Close();
  std::cout << "Processed " << processedEvents << " AMPT events; wrote "
            << outputFile << "." << std::endl;
}

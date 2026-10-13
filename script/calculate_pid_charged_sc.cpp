#include "../pidChargedSc.h"

#include "TFile.h"
#include "TH1D.h"
#include "TNamed.h"

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

// Save the two pure pair moments and the distinct two-PID/two-charged
// four-track moment needed by SC(v2,PID^2, v2,ch^2). No SC is computed here.
void calculate_pid_charged_sc(
    const char *inputConfigFile = "../config/cent_cfg_ampt.json",
    const char *outputFile = "pid_charged_sc_moments.root",
    int maxFilesPerConfig = -1, int maxConfigs = -1,
    const char *analysisConfigFile = "../config/config.json",
    int shardIndex = 0, int shardCount = 1) {
  using namespace ampt_analysis;
  const AnalysisConfig config = loadAnalysisConfig(analysisConfigFile);
  const auto &settings = config.pidChargedScOutput;
  const auto centralityEdges = makeAxisEdges(settings.centralityAxis);
  const auto groupEdges = makeAxisEdges(settings.bootstrapAxis);
  GFW gfw;
  const PidChargedScCorrConfigManager correlations(
      gfw, settings.innerAbsEta, settings.etaMax);
  std::array<std::unique_ptr<pid_charged_sc::Profiles>, 3> profiles;
  for (auto &profile : profiles) {
    profile = std::make_unique<pid_charged_sc::Profiles>(centralityEdges,
                                                       groupEdges);
  }
  TH1D eventCount("hEventCount", ";Selection;Events", 2, 0., 2.);
  eventCount.SetDirectory(nullptr);
  eventCount.GetXaxis()->SetBinLabel(1, "Read");
  eventCount.GetXaxis()->SetBinLabel(2, "Centrality accepted");
  TH1D centralityCount("hCent", ";Centrality (%);Events",
                       static_cast<int>(centralityEdges.size()) - 1,
                       centralityEdges.data());
  centralityCount.SetDirectory(nullptr);

  const Long64_t processedEvents = forEachConfiguredEvent(
      inputConfigFile, maxFilesPerConfig, maxConfigs, shardIndex, shardCount,
      [&](const Event &event, const CentralityConfig &,
          const EventIdentity &identity) {
        eventCount.Fill(0.5);
        const double centrality = centralityFromImpactParameter(event.imp, config);
        if (!std::isfinite(centrality) || centrality < centralityEdges.front() ||
            centrality >= centralityEdges.back()) {
          return;
        }
        eventCount.Fill(1.5);
        centralityCount.Fill(centrality);
        pid_charged_sc::fillGfw(gfw, event, config);
        const double group = sampleAxisCoordinate(
            groupEdges, bootstrapRandomValue(config.randomSeed, identity));
        for (std::size_t species = 0; species < profiles.size(); ++species) {
          profiles[species]->fill(
              centrality, group,
              pid_charged_sc::calculateMoments(gfw, correlations, species), settings);
        }
      });

  TFile output(outputFile, "RECREATE");
  if (output.IsZombie()) {
    throw std::runtime_error(std::string("Cannot create output file: ") + outputFile);
  }
  auto *directory = output.mkdir("pid-charged-sc");
  directory->cd();
  eventCount.Write();
  centralityCount.Write();
  TNamed("eventWeight", settings.eventWeight.c_str()).Write();
  TNamed("correlatorEngine", "GFW").Write();
  TNamed("harmonic", "2").Write();
  TNamed("subeventLayout", "A: negative eta; B: positive eta; +2,+2 in A; -2,-2 in B").Write();
  TNamed("eventSample", "Each species: all three moments require a valid mixed quadruplet").Write();
  TNamed("subsampleGroups",
         std::to_string(groupEdges.size() - 1).c_str()).Write();
  // Fixed metadata use TNamed: TParameter would be added during shard merges.
  std::ifstream input(analysisConfigFile);
  nlohmann::json configuration;
  input >> configuration;
  TNamed("analysisConfig", configuration.dump().c_str()).Write();
  for (const auto &species : speciesDefinitions()) {
    auto *speciesDirectory = directory->mkdir(species.speciesName);
    speciesDirectory->cd();
    profiles[static_cast<std::size_t>(species.species)]->write();
  }
  output.Close();
  std::cout << "Processed " << processedEvents << " AMPT events; wrote "
            << outputFile << "." << std::endl;
}

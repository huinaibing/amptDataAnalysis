#include "../analysisConfig.h"
#include "../analysisUtils.h"

#include "PWGCF/GenericFramework/Core/GFW.h"
#include "TFile.h"
#include "TH1D.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TRandom3.h"

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace ampt_c22_delta_pt_macro {
using ampt_analysis::CorrelationResult;
using ampt_analysis::SpeciesDefinition;

std::unique_ptr<TProfile3D>
makeC22Profile(const char *name, const char *title,
               const std::vector<double> &centralityEdges,
               const std::vector<double> &meanPtEdges,
               const std::vector<double> &bootstrapEdges) {
  return std::make_unique<TProfile3D>(
      name, title, static_cast<int>(centralityEdges.size()) - 1,
      centralityEdges.data(), static_cast<int>(meanPtEdges.size()) - 1,
      meanPtEdges.data(), static_cast<int>(bootstrapEdges.size()) - 1,
      bootstrapEdges.data());
}

std::unique_ptr<TProfile2D>
makeMeanPtProfile(const char *name, const char *title,
                  const std::vector<double> &centralityEdges,
                  const std::vector<double> &bootstrapEdges) {
  return std::make_unique<TProfile2D>(
      name, title, static_cast<int>(centralityEdges.size()) - 1,
      centralityEdges.data(), static_cast<int>(bootstrapEdges.size()) - 1,
      bootstrapEdges.data());
}

void fillC22Profile(TProfile3D &profile, double centrality, double meanPt,
                    double bootstrap, const CorrelationResult &correlation) {
  if (correlation.isPhysical()) {
    // Mean pT is a coordinate; only the GFW pair count enters the weight.
    profile.Fill(centrality, meanPt, bootstrap, correlation.value(),
                 correlation.pairs);
  }
}

struct C22SpeciesOutput {
  const SpeciesDefinition *definition = nullptr;
  std::unique_ptr<TProfile2D> meanPt;
  std::unique_ptr<TProfile3D> refRef;
  std::unique_ptr<TProfile3D> signal;
};

C22SpeciesOutput makeSpeciesOutput(const SpeciesDefinition &definition,
                                   bool usePure,
                                   const std::vector<double> &centralityEdges,
                                   const std::vector<double> &meanPtEdges,
                                   const std::vector<double> &bootstrapEdges) {
  const std::string species = definition.speciesName;
  const std::string axisLabel = definition.axisLabel;

  C22SpeciesOutput output;
  output.definition = &definition;
  output.meanPt = makeMeanPtProfile(
      ("hMeanPt" + species).c_str(),
      (species + " mean p_{T};Centrality (%);Bootstrap subsample;[p_{T}]_{" +
       axisLabel + "} (GeV/#it{c})")
          .c_str(),
      centralityEdges, bootstrapEdges);
  output.refRef =
      makeC22Profile(("c22dmeanpt" + species + "RefRef").c_str(),
                     ("Ref-ref c_{2}{2} versus " + species +
                      " event [p_{T}];Centrality (%);[p_{T}]_{" + axisLabel +
                      "} (GeV/#it{c});Bootstrap subsample;c_{2}{2}")
                         .c_str(),
                     centralityEdges, meanPtEdges, bootstrapEdges);

  const std::string mode = usePure ? "Pure" : "POIRef";
  const std::string correlationLabel = usePure ? "POI-POI" : "POI-ref";
  output.signal =
      makeC22Profile(("c22dmeanpt" + species + mode).c_str(),
                     (species + " " + correlationLabel + " c_{2}{2} versus " +
                      species + " event [p_{T}];Centrality (%);[p_{T}]_{" +
                      axisLabel + "} (GeV/#it{c});Bootstrap subsample;c_{2}{2}")
                         .c_str(),
                     centralityEdges, meanPtEdges, bootstrapEdges);
  return output;
}

void writeObject(TDirectory &directory, TObject &object) {
  directory.WriteTObject(&object, object.GetName());
}
} // namespace ampt_c22_delta_pt_macro

/**
 * Produce the same c22DeltaPt histogram layout as processDataC22DeltaPt.
 *
 * Test example (one ROOT file from the first centrality configuration):
 * root -l -b -q 'calculate_c22deltapt.cpp("../config/cent_cfg.json",
 * "/tmp/amptC22.root", 1, 1, false)'
 *
 * @param inputConfigFile JSON input configuration shared with
 * calculate_v2ptrho.cpp.
 * @param outputFile AnalysisResults-style output ROOT file.
 * @param maxFilesPerConfig Negative means use every configured input file.
 * @param maxConfigs Negative means process every centrality configuration.
 * @param usePure -1 uses JSON, 1 writes PID POI-POI, and 0 writes PID POI-ref.
 * @param analysisConfigFile Physics selections and histogram axes.
 */
void calculate_c22deltapt(
    const char *inputConfigFile = "../config/cent_cfg.json",
    const char *outputFile = "myAnalysisResultC22DeltaPt.root",
    int maxFilesPerConfig = -1, int maxConfigs = -1, int usePure = -1,
    const char *analysisConfigFile = "../config/config.json") {
  using namespace ampt_analysis;
  using namespace ampt_c22_delta_pt_macro;

  AnalysisConfig config = loadAnalysisConfig(analysisConfigFile);
  config.strictPtBounds = config.c22DeltaPtOutput.strictPtBounds;
  config.strictEtaBounds = config.c22DeltaPtOutput.strictEtaBounds;
  const bool usePureMode =
      usePure < 0 ? config.c22UsePure : static_cast<bool>(usePure);
  const auto centralityEdges =
      makeAxisEdges(config.c22DeltaPtOutput.centralityAxis);
  const auto meanPtEdges = makeAxisEdges(config.c22DeltaPtOutput.meanPtAxis);
  const auto bootstrapEdges =
      makeAxisEdges(config.c22DeltaPtOutput.bootstrapAxis);

  auto chargedC22 = makeC22Profile(
      "c22dmeanptCharged",
      "Charged c_{2}{2} versus charged event [p_{T}];Centrality "
      "(%);[p_{T}]_{ch} (GeV/#it{c});Bootstrap subsample;c_{2}{2}",
      centralityEdges, meanPtEdges, bootstrapEdges);
  auto chargedMeanPt =
      makeMeanPtProfile("hMeanPtCharged",
                        "Charged mean p_{T};Centrality (%);Bootstrap "
                        "subsample;[p_{T}]_{ch} (GeV/#it{c})",
                        centralityEdges, bootstrapEdges);

  std::vector<C22SpeciesOutput> speciesOutputs;
  speciesOutputs.reserve(speciesDefinitions().size());
  for (const auto &definition : speciesDefinitions()) {
    speciesOutputs.emplace_back(makeSpeciesOutput(
        definition, usePureMode, centralityEdges, meanPtEdges, bootstrapEdges));
  }

  TH1D eventCount("processDataC22DeltaPt", "", 14, 0., 14.);
  eventCount.GetXaxis()->SetBinLabel(1, "Filtered event");
  eventCount.GetXaxis()->SetBinLabel(2, "after sel8");

  TRandom3 random(config.randomSeed);
  GFW gfw;
  CorrConfigManager manager(&gfw, config.flowEtaGap, config.flowEtaMax);

  const Long64_t processedEvents = forEachConfiguredEvent(
      inputConfigFile, maxFilesPerConfig, maxConfigs,
      [&](const Event &event, const CentralityConfig &) {
        eventCount.Fill(0.5);
        eventCount.Fill(
            1.5); // AMPT has no reconstructed-event sel8 equivalent.

        gfw.Clear();
        const EventSamples samples =
            fillGfwAndCollectSamples(gfw, event, config);

        const double centrality =
            centralityFromImpactParameter(event.imp, config);
        const double bootstrap =
            sampleAxisCoordinate(bootstrapEdges, random.Rndm());
        const Event::PtMoments &chargedMoments = samples.charged;
        if (chargedMoments.count == 0) {
          return;
        }

        const CorrelationResult refRef =
            calculateCorrelation(gfw, manager, CorrType::Ref08Gap22);
        chargedMeanPt->Fill(centrality, bootstrap, chargedMoments.mean(),
                            static_cast<double>(chargedMoments.count));
        fillC22Profile(*chargedC22, centrality, chargedMoments.mean(),
                       bootstrap, refRef);

        for (auto &output : speciesOutputs) {
          const SpeciesDefinition &definition = *output.definition;
          const Event::PtMoments &pidMoments =
              samples.forSpecies(definition.species);
          if (pidMoments.count == 0) {
            continue;
          }

          output.meanPt->Fill(centrality, bootstrap, pidMoments.mean(),
                              static_cast<double>(pidMoments.count));
          fillC22Profile(*output.refRef, centrality, pidMoments.mean(),
                         bootstrap, refRef);

          const CorrelationResult signal =
              usePureMode
                  ? calculateCorrelation(gfw, manager, definition.pure)
                  : calculateCombinedCorrelation(gfw, manager,
                                                 definition.poiRefA,
                                                 definition.poiRefB);
          fillC22Profile(*output.signal, centrality, pidMoments.mean(),
                         bootstrap, signal);
        }
      });

  TFile outputFileHandle(outputFile, "RECREATE");
  if (outputFileHandle.IsZombie()) {
    throw std::runtime_error(std::string("Cannot create output file: ") +
                             outputFile);
  }

  TDirectory *taskDirectory = outputFileHandle.mkdir("pid-flow-pt-corr");
  TDirectory *eventCountDirectory = taskDirectory->mkdir("hEventCount");
  writeObject(*eventCountDirectory, eventCount);

  TDirectory *c22Directory = taskDirectory->mkdir("c22DeltaPt");
  writeObject(*c22Directory, *chargedC22);
  writeObject(*c22Directory, *chargedMeanPt);
  for (auto &output : speciesOutputs) {
    writeObject(*c22Directory, *output.meanPt);
  }
  for (auto &output : speciesOutputs) {
    writeObject(*c22Directory, *output.refRef);
  }
  for (auto &output : speciesOutputs) {
    writeObject(*c22Directory, *output.signal);
  }

  outputFileHandle.Close();
  std::cout << "Processed " << processedEvents << " AMPT events; wrote "
            << outputFile << " in " << (usePureMode ? "Pure" : "POI-ref")
            << " mode." << std::endl;
}

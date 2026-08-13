#include "../analysisConfig.h"
#include "../analysisUtils.h"
#include "../utils.h"

#include "Framework/HistogramSpec.h"
#include "PWGCF/GenericFramework/Core/FlowContainer.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TProfile3D.h"
#include "TRandom3.h"

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace ampt_v2_pt_rho_macro {
using ampt_analysis::CorrelationResult;
using ampt_analysis::SpeciesDefinition;

std::unique_ptr<TObjArray> makeChargedProfileNames() {
  auto names = std::make_unique<TObjArray>();
  names->SetOwner(true);
  for (const char *name : {"c22", "c32", "c24", "c34", "c22Full",
                           "c22TrackWeight", "c32TrackWeight", "c24TrackWeight",
                           "c34TrackWeight", "c22FullTrackWeight", "covV2Pt",
                           "covV3Pt", "ptSquareAve", "ptAve", "hMeanPt"}) {
    names->Add(new TNamed(name, name));
  }
  return names;
}

std::unique_ptr<TObjArray> makePidProfileNames(const TObjArray &chargedNames) {
  auto names = std::unique_ptr<TObjArray>(
      static_cast<TObjArray *>(chargedNames.Clone()));
  names->SetOwner(true);
  for (const char *name :
       {"c22pure", "c32pure", "covV2PtPID", "c22TrackWeightPID"}) {
    names->Add(new TNamed(name, name));
  }
  return names;
}

std::unique_ptr<TProfile3D>
makeV2PtRhoProfile(const std::string &name,
                   const std::vector<double> &meanPtEdges,
                   const std::vector<double> &centralityEdges,
                   const std::vector<double> &bootstrapEdges) {
  return std::make_unique<TProfile3D>(
      name.c_str(), name.c_str(), static_cast<int>(meanPtEdges.size()) - 1,
      meanPtEdges.data(), static_cast<int>(centralityEdges.size()) - 1,
      centralityEdges.data(), static_cast<int>(bootstrapEdges.size()) - 1,
      bootstrapEdges.data());
}

struct V2PtRhoSpeciesOutput {
  const SpeciesDefinition *definition = nullptr;
  std::unique_ptr<FlowContainer> flow;
  std::unique_ptr<TProfile3D> poiRef;
  std::unique_ptr<TProfile3D> refRef;
  std::unique_ptr<TProfile3D> pure;
  std::unique_ptr<TProfile3D> meanPt;
};

V2PtRhoSpeciesOutput
makeSpeciesOutput(const SpeciesDefinition &definition,
                  const o2::framework::AxisSpec &centralityAxis,
                  TObjArray &profileNames,
                  const std::vector<double> &meanPtEdges,
                  const std::vector<double> &centralityEdges,
                  const std::vector<double> &bootstrapEdges, int nBootstrap) {
  V2PtRhoSpeciesOutput output;
  output.definition = &definition;
  output.flow = std::make_unique<FlowContainer>(definition.flowContainerName);
  output.flow->Initialize(&profileNames, centralityAxis, nBootstrap);

  const std::string baseName = definition.speciesName;
  output.poiRef = makeV2PtRhoProfile("h" + baseName, meanPtEdges,
                                     centralityEdges, bootstrapEdges);
  output.refRef = makeV2PtRhoProfile(
      "hCharged" + baseName + "Full", meanPtEdges, centralityEdges,
      bootstrapEdges);
  output.pure = makeV2PtRhoProfile("h" + baseName + baseName, meanPtEdges,
                                   centralityEdges, bootstrapEdges);
  output.meanPt = makeV2PtRhoProfile(
      "h" + baseName + "Meanpt", meanPtEdges, centralityEdges,
      bootstrapEdges);
  return output;
}

void writeObject(TDirectory &directory, TObject &object) {
  directory.WriteTObject(&object, object.GetName());
}
} // namespace ampt_v2_pt_rho_macro

/**
 * Produce the processData-compatible v2-pT-rho FlowContainer and
 * meanptCentNbs output.
 *
 * For a quick test, set maxFilesPerConfig=1 and maxConfigs=1. Negative limits
 * preserve the original full-input behavior.
 */
void calculate_v2ptrho(const char *inputConfigFile = "../config/cent_cfg.json",
                       const char *outputFile = "myAnalysisResultV2PtRho.root",
                       int maxFilesPerConfig = -1, int maxConfigs = -1,
                       const char *analysisConfigFile =
                           "../config/config.json") {
  using namespace ampt_analysis;
  using namespace ampt_v2_pt_rho_macro;

  AnalysisConfig config = loadAnalysisConfig(analysisConfigFile);
  config.strictPtBounds = config.v2PtRhoOutput.strictPtBounds;
  config.strictEtaBounds = config.v2PtRhoOutput.strictEtaBounds;
  const auto meanPtEdges = makeAxisEdges(config.v2PtRhoOutput.meanPtAxis);
  const auto bootstrapEdges =
      makeAxisEdges(config.v2PtRhoOutput.bootstrapAxis);
  const auto centralityEdges =
      makeAxisEdges(config.v2PtRhoOutput.centralityAxis);
  const int nBootstrap = axisBinCount(config.v2PtRhoOutput.bootstrapAxis);
  const o2::framework::AxisSpec centralityAxis{centralityEdges,
                                               "Centrality (%)"};

  auto chargedNames = makeChargedProfileNames();
  auto pidNames = makePidProfileNames(*chargedNames);

  FlowContainer chargedFlow("FlowContainerCharged");
  chargedFlow.Initialize(chargedNames.get(), centralityAxis, nBootstrap);

  std::vector<V2PtRhoSpeciesOutput> speciesOutputs;
  speciesOutputs.reserve(speciesDefinitions().size());
  for (const auto &definition : speciesDefinitions()) {
    speciesOutputs.emplace_back(
        makeSpeciesOutput(definition, centralityAxis, *pidNames, meanPtEdges,
                          centralityEdges, bootstrapEdges, nBootstrap));
  }

  TRandom3 random(config.randomSeed);
  GFW gfw;
  CorrConfigManager manager(&gfw, config.flowEtaGap, config.flowEtaMax);

  const Long64_t processedEvents = forEachConfiguredEvent(
      inputConfigFile, maxFilesPerConfig, maxConfigs,
      [&](const Event &event, const CentralityConfig &) {
        gfw.Clear();
        const EventSamples samples =
            fillGfwAndCollectSamples(gfw, event, config);

        const double centrality =
            centralityFromImpactParameter(event.imp, config);
        const double randomValue = random.Rndm();
        const double bootstrap =
            sampleAxisCoordinate(bootstrapEdges, randomValue);
        const Event::PtMoments &chargedMoments = samples.charged;

        const CorrelationResult chargedGap =
            calculateCorrelation(gfw, manager, CorrType::Ref08Gap22);
        fillFlowProfile(chargedFlow, "c22", centrality, chargedGap,
                        randomValue);
        fillFlowProfile(chargedFlow, "c24", centrality,
                        calculateCorrelation(gfw, manager, CorrType::Ref0Gap24),
                        randomValue);
        fillFlowProfile(chargedFlow, "c22Full", centrality,
                        calculateCorrelation(gfw, manager, CorrType::Ref0Gap22),
                        randomValue);
        fillFlowMeanPtProduct(chargedFlow, "covV2Pt", centrality, chargedGap,
                              chargedMoments, randomValue);
        fillTrackWeightedFlow(chargedFlow, "c22TrackWeight", centrality,
                              chargedGap, chargedMoments, randomValue);
        fillMeanPtMoments(chargedFlow, centrality, chargedMoments, randomValue);

        for (auto &output : speciesOutputs) {
          const SpeciesDefinition &definition = *output.definition;
          const Event::PtMoments &pidMoments =
              samples.forSpecies(definition.species);

          // Keep the two eta orientations as separate FlowContainer fills;
          // this preserves the processData profile-error bookkeeping.
          const CorrelationResult poiRefA =
              calculateCorrelation(gfw, manager, definition.poiRefA);
          const CorrelationResult poiRefB =
              calculateCorrelation(gfw, manager, definition.poiRefB);
          fillFlowProfile(*output.flow, "c22", centrality, poiRefA,
                          randomValue);
          fillFlowProfile(*output.flow, "c22", centrality, poiRefB,
                          randomValue);

          const CorrelationResult pure =
              calculateCorrelation(gfw, manager, definition.pure);
          fillFlowProfile(*output.flow, "c22pure", centrality, pure,
                          randomValue);
          fillFlowProfile(
              *output.flow, "c22Full", centrality,
              calculateCorrelation(gfw, manager, definition.fullTwoA),
              randomValue);
          fillFlowProfile(
              *output.flow, "c22Full", centrality,
              calculateCorrelation(gfw, manager, definition.fullTwoB),
              randomValue);
          fillFlowProfile(*output.flow, "c24", centrality,
                          calculateCorrelation(gfw, manager, definition.fourA),
                          randomValue);
          fillFlowProfile(*output.flow, "c24", centrality,
                          calculateCorrelation(gfw, manager, definition.fourB),
                          randomValue);

          fillFlowMeanPtProduct(*output.flow, "covV2Pt", centrality, pure,
                                chargedMoments, randomValue);
          fillTrackWeightedFlow(*output.flow, "c22TrackWeight", centrality,
                                pure, chargedMoments, randomValue);
          fillFlowMeanPtProduct(*output.flow, "covV2PtPID", centrality, pure,
                                pidMoments, randomValue);
          fillTrackWeightedFlow(*output.flow, "c22TrackWeightPID", centrality,
                                pure, pidMoments, randomValue);

          if (pidMoments.count > 1) {
            const CorrelationResult poiRef{poiRefA.numerator +
                                               poiRefB.numerator,
                                           poiRefA.pairs + poiRefB.pairs};
            fillV2PtRhoMeanPtProfiles(
                centrality, bootstrap, pidMoments, chargedGap, poiRef, pure,
                *output.poiRef, *output.refRef, *output.pure, *output.meanPt);
            fillMeanPtMoments(*output.flow, centrality, pidMoments,
                              randomValue);
          }
        }
      });

  TFile outputFileHandle(outputFile, "RECREATE");
  if (outputFileHandle.IsZombie()) {
    throw std::runtime_error(std::string("Cannot create output file: ") +
                             outputFile);
  }

  TDirectory *taskDirectory = outputFileHandle.mkdir("pid-flow-pt-corr");
  writeObject(*taskDirectory, chargedFlow);
  for (auto &output : speciesOutputs) {
    writeObject(*taskDirectory, *output.flow);
  }

  TDirectory *meanPtDirectory = taskDirectory->mkdir("meanptCentNbs");
  for (auto &output : speciesOutputs) {
    writeObject(*meanPtDirectory, *output.poiRef);
    writeObject(*meanPtDirectory, *output.refRef);
    writeObject(*meanPtDirectory, *output.pure);
  }
  for (auto &output : speciesOutputs) {
    writeObject(*meanPtDirectory, *output.meanPt);
  }

  outputFileHandle.Close();
  std::cout << "Processed " << processedEvents << " AMPT events; wrote "
            << outputFile << "." << std::endl;
}

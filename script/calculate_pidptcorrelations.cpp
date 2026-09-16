#include "../analysisConfig.h"
#include "../analysisUtils.h"

#include "Framework/HistogramSpec.h"
#include "PWGCF/GenericFramework/Core/FlowContainer.h"
#include "TFile.h"
#include "TH1D.h"
#include "TObjArray.h"
#include "TRandom3.h"

#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

namespace ampt_pid_pt_correlations_macro
{
  struct PidPtMoments
  {
    double sumWeight = 0.;
    double sumPt = 0.;
    double sumWeight2 = 0.;
    double sumPtWeight2 = 0.;
    double sumPt2Weight2 = 0.;

    void add(double pt, double particleWeight = 1.)
    {
      const double weight2 = particleWeight * particleWeight;
      sumWeight += particleWeight;
      sumPt += particleWeight * pt;
      sumWeight2 += weight2;
      sumPtWeight2 += weight2 * pt;
      sumPt2Weight2 += weight2 * pt * pt;
    }
  };

  std::unique_ptr<TObjArray> makeProfileNames()
  {
    auto names = std::make_unique<TObjArray>();
    names->SetOwner(true);
    for (const char *name : {
             "meanPtPi", "meanPtKa", "meanPtPr",
             "ptProductPiPi", "ptPiInPiPiPairs",
             "ptProductKaKa", "ptKaInKaKaPairs",
             "ptProductPrPr", "ptPrInPrPrPairs",
             "ptProductPiKa", "ptPiInPiKaPairs", "ptKaInPiKaPairs",
             "ptProductPiPr", "ptPiInPiPrPairs", "ptPrInPiPrPairs",
             "ptProductKaPr", "ptKaInKaPrPairs", "ptPrInKaPrPairs"})
    {
      names->Add(new TNamed(name, name));
    }
    return names;
  }

  void setEventCountLabels(TH1D &eventCount)
  {
    const std::array<const char *, 14> labels{
        "Filtered event", "after sel8", "after kTVXinTRD",
        "after kNoTimeFrameBorder", "after kNoITSROFrameBorder",
        "after kDoNoSameBunchPileup", "after kIsGoodZvtxFT0vsPV",
        "after kNoCollInTimeRangeStandard", "after kIsGoodITSLayersAll",
        "after MultPVCut", "after TPC occupancy cut", "after V0AT0Acut",
        "after IRmincut", "after IRmaxcut"};
    for (std::size_t index = 0; index < labels.size(); ++index)
    {
      eventCount.GetXaxis()->SetBinLabel(static_cast<int>(index) + 1,
                                         labels[index]);
    }
  }

  void fillSelfMoments(FlowContainer &output, const char *productName,
                       const char *singlePtName, double centrality,
                       double randomValue, const PidPtMoments &moments)
  {
    const double pairWeight =
        moments.sumWeight * moments.sumWeight - moments.sumWeight2;
    if (pairWeight <= 1.e-3)
    {
      return;
    }
    output.FillProfile(
        productName, centrality,
        (moments.sumPt * moments.sumPt - moments.sumPt2Weight2) / pairWeight,
        pairWeight, randomValue);
    output.FillProfile(
        singlePtName, centrality,
        (moments.sumWeight * moments.sumPt - moments.sumPtWeight2) /
            pairWeight,
        pairWeight, randomValue);
  }

  void fillCrossMoments(FlowContainer &output, const char *productName,
                        const char *alphaPtName, const char *betaPtName,
                        double centrality, double randomValue,
                        const PidPtMoments &alpha,
                        const PidPtMoments &beta)
  {
    if (alpha.sumWeight <= 0. || beta.sumWeight <= 0.)
    {
      return;
    }
    const double pairWeight = alpha.sumWeight * beta.sumWeight;
    output.FillProfile(productName, centrality,
                       alpha.sumPt * beta.sumPt / pairWeight, pairWeight,
                       randomValue);
    output.FillProfile(alphaPtName, centrality,
                       alpha.sumPt / alpha.sumWeight, pairWeight, randomValue);
    output.FillProfile(betaPtName, centrality, beta.sumPt / beta.sumWeight,
                       pairWeight, randomValue);
  }

  void writeObject(TDirectory &directory, TObject &object)
  {
    directory.WriteTObject(&object, object.GetName());
  }
} // namespace ampt_pid_pt_correlations_macro

/**
 * Produce the processPidPtCorrelations-compatible PID mean-pT correlations.
 *
 * AMPT truth PID replaces detector PID and every accepted particle has unit
 * efficiency weight. The full flow-track eta interval is used, matching
 * processPidPtCorrelations rather than the narrower mean-pT interval.
 *
 * Quick test:
 * root -l -b -q 'calculate_pidptcorrelations.cpp(
 * "../config/cent_cfg.json", "/tmp/amptPidPt.root", 1, 1,
 * "../config/config.json")'
 */
void calculate_pidptcorrelations(
    const char *inputConfigFile = "../config/cent_cfg.json",
    const char *outputFile = "myAnalysisResultPidPtCorrelations.root",
    int maxFilesPerConfig = -1, int maxConfigs = -1,
    const char *analysisConfigFile = "../config/config.json")
{
  using namespace ampt_analysis;
  using namespace ampt_pid_pt_correlations_macro;

  const AnalysisConfig config = loadAnalysisConfig(analysisConfigFile);
  const PidPtCorrelationsOutputConfig &outputConfig =
      config.pidPtCorrelationsOutput;
  const auto centralityEdges =
      makeAxisEdges(outputConfig.centralityAxis);
  const int nBootstrap = axisBinCount(outputConfig.bootstrapAxis);
  const o2::framework::AxisSpec centralityAxis{centralityEdges,
                                               "Centrality (%)"};

  auto profileNames = makeProfileNames();
  FlowContainer pidPtCorrelations("FlowContainerPidPtCorr");
  pidPtCorrelations.Initialize(profileNames.get(), centralityAxis, nBootstrap);

  TH1D eventCount("processPidPtCorrelations", "", 14, 0., 14.);
  setEventCountLabels(eventCount);

  TRandom3 random(config.randomSeed);
  const Long64_t processedEvents = forEachConfiguredEvent(
      inputConfigFile, maxFilesPerConfig, maxConfigs,
      [&](const Event &event, const CentralityConfig &)
      {
        eventCount.Fill(0.5);
        if (event.particles.empty())
        {
          return;
        }
        eventCount.Fill(1.5);

        std::array<PidPtMoments, 3> moments;
        for (const auto &track : event.particles)
        {
          if (!isChargedPdg(track.pdgPid))
          {
            continue;
          }
          const bool withinEta =
              outputConfig.strictEtaBounds
                  ? std::abs(track.GetEta()) < config.flowEtaMax
                  : std::abs(track.GetEta()) <= config.flowEtaMax;
          if (!withinEta)
          {
            continue;
          }
          const int absPdg = std::abs(track.pdgPid);
          const SpeciesDefinition *species = findSpecies(absPdg);
          if (!species ||
              !ptRangeForPdg(absPdg, config)
                   .contains(track.GetPt(), outputConfig.strictPtBounds))
          {
            continue;
          }
          moments.at(static_cast<std::size_t>(species->species))
              .add(track.GetPt());
        }

        const double centrality =
            centralityFromImpactParameter(event.imp, config);
        const double randomValue = random.Rndm();
        constexpr std::size_t pion = static_cast<std::size_t>(Species::Pion);
        constexpr std::size_t kaon = static_cast<std::size_t>(Species::Kaon);
        constexpr std::size_t proton =
            static_cast<std::size_t>(Species::Proton);

        if (moments[pion].sumWeight > 0.)
        {
          pidPtCorrelations.FillProfile(
              "meanPtPi", centrality,
              moments[pion].sumPt / moments[pion].sumWeight,
              moments[pion].sumWeight, randomValue);
        }
        if (moments[kaon].sumWeight > 0.)
        {
          pidPtCorrelations.FillProfile(
              "meanPtKa", centrality,
              moments[kaon].sumPt / moments[kaon].sumWeight,
              moments[kaon].sumWeight, randomValue);
        }
        if (moments[proton].sumWeight > 0.)
        {
          pidPtCorrelations.FillProfile(
              "meanPtPr", centrality,
              moments[proton].sumPt / moments[proton].sumWeight,
              moments[proton].sumWeight, randomValue);
        }

        fillSelfMoments(pidPtCorrelations, "ptProductPiPi",
                        "ptPiInPiPiPairs", centrality, randomValue,
                        moments[pion]);
        fillSelfMoments(pidPtCorrelations, "ptProductKaKa",
                        "ptKaInKaKaPairs", centrality, randomValue,
                        moments[kaon]);
        fillSelfMoments(pidPtCorrelations, "ptProductPrPr",
                        "ptPrInPrPrPairs", centrality, randomValue,
                        moments[proton]);

        fillCrossMoments(pidPtCorrelations, "ptProductPiKa",
                         "ptPiInPiKaPairs", "ptKaInPiKaPairs", centrality,
                         randomValue, moments[pion], moments[kaon]);
        fillCrossMoments(pidPtCorrelations, "ptProductPiPr",
                         "ptPiInPiPrPairs", "ptPrInPiPrPairs", centrality,
                         randomValue, moments[pion], moments[proton]);
        fillCrossMoments(pidPtCorrelations, "ptProductKaPr",
                         "ptKaInKaPrPairs", "ptPrInKaPrPairs", centrality,
                         randomValue, moments[kaon], moments[proton]);
      });

  TFile outputFileHandle(outputFile, "RECREATE");
  if (outputFileHandle.IsZombie())
  {
    throw std::runtime_error(std::string("Cannot create output file: ") +
                             outputFile);
  }

  TDirectory *taskDirectory = outputFileHandle.mkdir("pid-flow-pt-corr");
  writeObject(*taskDirectory, pidPtCorrelations);
  TDirectory *eventCountDirectory = taskDirectory->mkdir("hEventCount");
  writeObject(*eventCountDirectory, eventCount);

  outputFileHandle.Close();
  std::cout << "Processed " << processedEvents << " AMPT events; wrote "
            << outputFile << "." << std::endl;
}

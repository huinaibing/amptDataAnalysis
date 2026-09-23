#include "../analysisConfig.h"
#include "../analysisUtils.h"

#include "Framework/HistogramSpec.h"
#include "PWGCF/GenericFramework/Core/FlowContainer.h"
#include "PWGCF/GenericFramework/Core/GFW.h"

#include "TAxis.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TProfile.h"
#include "TRandom3.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ampt_flow_pbpb_pikp_macro
{
  constexpr const char *kTaskDirectory = "flow-pbpb-pikp_TPC_TOF";

  void writeObject(TDirectory &directory, TObject &object)
  {
    directory.WriteTObject(&object, object.GetName());
  }

  void fillProfile(GFW &gfw, const GFW::CorrConfig &correlation,
                   TProfile &profile, double activity)
  {
    if (correlation.pTDif)
    {
      return;
    }
    const double denominator =
        std::real(gfw.Calculate(correlation, 0, true));
    if (denominator == 0.)
    {
      return;
    }
    const double value =
        std::real(gfw.Calculate(correlation, 0, false)) / denominator;
    if (std::fabs(value) < 1.)
    {
      profile.Fill(activity, value, denominator);
    }
  }

  void fillFlowContainer(GFW &gfw, const GFW::CorrConfig &correlation,
                         FlowContainer &output, double activity,
                         double random, int nPtBins)
  {
    if (!correlation.pTDif)
    {
      const double denominator =
          std::real(gfw.Calculate(correlation, 0, true));
      if (denominator == 0.)
      {
        return;
      }
      const double value =
          std::real(gfw.Calculate(correlation, 0, false)) / denominator;
      if (std::fabs(value) < 1.)
      {
        output.FillProfile(correlation.Head.c_str(), activity, value,
                           denominator, random);
      }
      return;
    }

    for (int ptBin = 1; ptBin <= nPtBins; ++ptBin)
    {
      const double denominator =
          std::real(gfw.Calculate(correlation, ptBin - 1, true));
      if (denominator == 0.)
      {
        continue;
      }
      const double value =
          std::real(gfw.Calculate(correlation, ptBin - 1, false)) /
          denominator;
      if (std::fabs(value) < 1.)
      {
        const std::string profileName =
            correlation.Head + "_pt_" + std::to_string(ptBin);
        output.FillProfile(profileName.c_str(), activity, value, denominator,
                           random);
      }
    }
  }

  int pidIndexFromPdg(int pdgCode)
  {
    switch (std::abs(pdgCode))
    {
    case PDG_t::kPiPlus:
      return 1;
    case PDG_t::kKPlus:
      return 2;
    case PDG_t::kProton:
      return 3;
    default:
      return 0;
    }
  }

  double phiInZeroTwoPi(double phi)
  {
    return phi >= 0. ? phi : phi + TMath::TwoPi();
  }
} // namespace ampt_flow_pbpb_pikp_macro

/**
 * Reproduce flowPbpbPikp::process for AMPT truth. Particle identity comes
 * directly from the PDG code, and all GFW particle weights are one.
 */
void calculate_flowPbpbPikp(
    const char *inputConfigFile = "../config/cent_cfg.json",
    const char *outputFile = "myAnalysisResultFlowPbpbPikp.root",
    int maxFilesPerConfig = -1, int maxConfigs = -1,
    const char *analysisConfigFile = "../config/config.json")
{
  using namespace ampt_analysis;
  using namespace ampt_flow_pbpb_pikp_macro;

  const AnalysisConfig analysisConfig =
      loadAnalysisConfig(analysisConfigFile);
  const FlowPbpbPikpOutputConfig &config =
      analysisConfig.flowPbpbPikpOutput;
  const std::vector<double> ptEdges = makeAxisEdges(config.ptAxis);
  const std::vector<double> centralityEdges =
      makeAxisEdges(config.centralityAxis);
  const std::vector<double> nchEdges = makeAxisEdges(config.nchAxis);
  const std::vector<double> qaCentralityEdges =
      makeAxisEdges(config.qaCentralityAxis);
  const std::vector<double> qaMultiplicityEdges =
      makeAxisEdges(config.qaMultiplicityAxis);
  const std::vector<double> qaPhiEdges = makeAxisEdges(config.qaPhiAxis);
  const std::vector<double> qaEtaEdges = makeAxisEdges(config.qaEtaAxis);

  GFW gfw;
  const int nPtBins = static_cast<int>(ptEdges.size()) - 1;
  for (std::size_t i = 0; i < config.regionNames.size(); ++i)
  {
    const int nRegionPtBins = config.regionPtDifferential[i]
                                  ? nPtBins + 1
                                  : 1;
    gfw.AddRegion(config.regionNames[i], config.regionEtaMin[i],
                  config.regionEtaMax[i], nRegionPtBins,
                  config.regionMasks[i]);
  }

  std::vector<GFW::CorrConfig> correlations;
  correlations.reserve(config.correlations.size());
  for (std::size_t i = 0; i < config.correlations.size(); ++i)
  {
    correlations.emplace_back(gfw.GetCorrelatorConfig(
        config.correlations[i], config.correlationHeads[i],
        config.correlationPtDifferential[i] != 0));
  }
  if (correlations.empty())
  {
    throw std::runtime_error("flowPbpbPikp_output has no correlations");
  }
  gfw.CreateRegions();

  TAxis ptAxis(nPtBins, ptEdges.data());
  TObjArray profileNames;
  profileNames.SetOwner(true);
  for (const GFW::CorrConfig &correlation : correlations)
  {
    if (!correlation.pTDif)
    {
      profileNames.Add(
          new TNamed(correlation.Head.c_str(), correlation.Head.c_str()));
      continue;
    }
    for (int ptBin = 1; ptBin <= nPtBins; ++ptBin)
    {
      const std::string name =
          correlation.Head + "_pt_" + std::to_string(ptBin);
      profileNames.Add(new TNamed(name.c_str(),
                                  (correlation.Head + "_ptDiff").c_str()));
    }
  }

  FlowContainer flowContainer("FlowContainer");
  flowContainer.SetXAxis(&ptAxis);
  const std::vector<double> &flowActivityEdges =
      config.useNch ? nchEdges : centralityEdges;
  const o2::framework::AxisSpec flowActivityAxis{flowActivityEdges};
  flowContainer.Initialize(&profileNames, flowActivityAxis,
                           config.nBootstrap);

  TH1D hMult("hMult", "", static_cast<int>(qaMultiplicityEdges.size()) - 1,
             qaMultiplicityEdges.data());
  TH1D hCent("hCent", "", static_cast<int>(qaCentralityEdges.size()) - 1,
             qaCentralityEdges.data());
  TH1D hPhi("hPhi", "", static_cast<int>(qaPhiEdges.size()) - 1,
            qaPhiEdges.data());
  TH1D hEta("hEta", "", static_cast<int>(qaEtaEdges.size()) - 1,
            qaEtaEdges.data());
  TH1D hPt("hPt", "", nPtBins, ptEdges.data());
  TProfile c22FullCharged(
      "c22_full_ch", "", static_cast<int>(centralityEdges.size()) - 1,
      centralityEdges.data());
  TProfile c22FullChargedNch(
      "c22_full_ch_Nch", "", static_cast<int>(nchEdges.size()) - 1,
      nchEdges.data());

  TRandom3 random(analysisConfig.randomSeed);
  const Long64_t processedEvents = forEachConfiguredEvent(
      inputConfigFile, maxFilesPerConfig, maxConfigs,
      [&](const Event &event, const CentralityConfig &)
      {
        std::vector<const Track *> selectedTracks;
        selectedTracks.reserve(event.particles.size());
        for (const Track &track : event.particles)
        {
          const double pt = track.GetPt();
          if (!isChargedPdg(track.pdgPid) ||
              std::abs(track.GetEta()) >= config.cutEta ||
              pt <= config.cutPtPoiMin || pt >= config.cutPtPoiMax)
          {
            continue;
          }
          selectedTracks.push_back(&track);
        }

        const int nTotal = static_cast<int>(selectedTracks.size());
        if (nTotal < 1)
        {
          return;
        }

        const double centrality =
            centralityFromImpactParameter(event.imp, analysisConfig);
        hMult.Fill(nTotal);
        hCent.Fill(centrality);
        gfw.Clear();

        for (const Track *track : selectedTracks)
        {
          const double pt = track->GetPt();
          const double eta = track->GetEta();
          const double phi = phiInZeroTwoPi(track->GetPhi());
          const bool withinReference =
              pt > config.cutPtRefMin && pt < config.cutPtRefMax;
          const int ptBin = ptAxis.FindBin(pt) - 1;
          const int pidIndex = pidIndexFromPdg(track->pdgPid);

          hPhi.Fill(phi);
          hEta.Fill(eta);
          hPt.Fill(pt);

          if (withinReference)
          {
            gfw.Fill(eta, ptBin, phi, 1., 1);
          }
          gfw.Fill(eta, ptBin, phi, 1., 128);
          if (withinReference)
          {
            gfw.Fill(eta, ptBin, phi, 1., 256);
          }

          if (pidIndex != 0)
          {
            gfw.Fill(eta, ptBin, phi, 1., 1 << pidIndex);
            if (withinReference)
            {
              gfw.Fill(eta, ptBin, phi, 1., 1 << (pidIndex + 3));
            }
          }
        }

        fillProfile(gfw, correlations.front(), c22FullCharged, centrality);
        fillProfile(gfw, correlations.front(), c22FullChargedNch, nTotal);
        const double activity = config.useNch ? nTotal : centrality;
        const double randomValue = random.Rndm();
        for (const GFW::CorrConfig &correlation : correlations)
        {
          fillFlowContainer(gfw, correlation, flowContainer, activity,
                            randomValue, nPtBins);
        }
      });

  TFile outputFileHandle(outputFile, "RECREATE");
  if (outputFileHandle.IsZombie())
  {
    throw std::runtime_error(std::string("Cannot create output file: ") +
                             outputFile);
  }
  TDirectory *taskDirectory = outputFileHandle.mkdir(kTaskDirectory);
  writeObject(*taskDirectory, hMult);
  writeObject(*taskDirectory, hCent);
  writeObject(*taskDirectory, hPhi);
  writeObject(*taskDirectory, hEta);
  writeObject(*taskDirectory, hPt);
  writeObject(*taskDirectory, c22FullCharged);
  writeObject(*taskDirectory, c22FullChargedNch);
  writeObject(*taskDirectory, flowContainer);

  outputFileHandle.Close();
  std::cout << "Processed " << processedEvents << " AMPT events; wrote "
            << outputFile << "." << std::endl;
}

#include "../analysisConfig.h"
#include "../analysisUtils.h"

#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TMath.h"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ampt_pt_spectra_macro
{
  constexpr const char *kTaskDirectory = "pt-spectra";

  void writeObject(TDirectory &directory, TObject &object)
  {
    directory.WriteTObject(&object, object.GetName());
  }

  double phiInZeroTwoPi(double phi)
  {
    return phi >= 0. ? phi : phi + TMath::TwoPi();
  }
} // namespace ampt_pt_spectra_macro

/**
 * Fill charged and identified-particle pT spectra from AMPT truth.
 *
 * Event and track iteration follow the same streaming pattern as the other
 * macros in this directory. The selection cuts and all histogram axes are
 * read from the "pt_spectra_output" block in config/config.json.
 */
void calculate_ptspectra(
    const char *inputConfigFile = "../config/cent_cfg.json",
    const char *outputFile = "myAnalysisResultPtSpectra.root",
    int maxFilesPerConfig = -1, int maxConfigs = -1,
    const char *analysisConfigFile = "../config/config.json")
{
  using namespace ampt_analysis;
  using namespace ampt_pt_spectra_macro;

  AnalysisConfig config = loadAnalysisConfig(analysisConfigFile);
  config.strictPtBounds = config.ptSpectraOutput.strictPtBounds;
  config.strictEtaBounds = config.ptSpectraOutput.strictEtaBounds;

  const auto ptEdges = makeAxisEdges(config.ptSpectraOutput.ptAxis);
  const auto centralityEdges =
      makeAxisEdges(config.ptSpectraOutput.centralityAxis);
  const auto phiEdges = makeAxisEdges(config.ptSpectraOutput.phiAxis);
  const auto etaEdges = makeAxisEdges(config.ptSpectraOutput.etaAxis);

  const int nPtBins = static_cast<int>(ptEdges.size()) - 1;
  const int nCentralityBins = static_cast<int>(centralityEdges.size()) - 1;
  const int nEtaBins = static_cast<int>(etaEdges.size()) - 1;

  TH1D hPtCharged("hPtCharged",
                  "Charged particle p_{T} spectrum;p_{T} (GeV/#it{c});"
                  "Tracks",
                  static_cast<int>(ptEdges.size()) - 1, ptEdges.data());
  TH1D hPtPion("hPtPion", "#pi p_{T} spectrum;p_{T} (GeV/#it{c});Tracks",
               static_cast<int>(ptEdges.size()) - 1, ptEdges.data());
  TH1D hPtKaon("hPtKaon", "K p_{T} spectrum;p_{T} (GeV/#it{c});Tracks",
               static_cast<int>(ptEdges.size()) - 1, ptEdges.data());
  TH1D hPtProton("hPtProton", "p p_{T} spectrum;p_{T} (GeV/#it{c});Tracks",
                 static_cast<int>(ptEdges.size()) - 1, ptEdges.data());
  TH1D hCent("hCent",
             "Event centrality distribution (AMPT impact parameter);"
             "Centrality (%);Events",
             static_cast<int>(centralityEdges.size()) - 1,
             centralityEdges.data());
  TH1D hPhi("hPhi", "Track #phi distribution;#phi (rad);Tracks",
            static_cast<int>(phiEdges.size()) - 1, phiEdges.data());
  TH1D hEta("hEta", "Track #eta distribution;#eta;Tracks",
            static_cast<int>(etaEdges.size()) - 1, etaEdges.data());

  TH1D hEventCount("hEventCount", "Events accepted;all events;Events", 1,
                   0., 1.);

  TH2D hPtEtaCharged("hPtEtaCharged",
                     "Charged raw p_{T}-#eta counts;p_{T} (GeV/#it{c});"
                     "#eta",
                     nPtBins, ptEdges.data(), nEtaBins, etaEdges.data());
  TH2D hPtEtaPion("hPtEtaPion", "#pi raw p_{T}-#eta counts;"
                               "p_{T} (GeV/#it{c});#eta",
                  nPtBins, ptEdges.data(), nEtaBins, etaEdges.data());
  TH2D hPtEtaKaon("hPtEtaKaon", "K raw p_{T}-#eta counts;"
                                "p_{T} (GeV/#it{c});#eta",
                  nPtBins, ptEdges.data(), nEtaBins, etaEdges.data());
  TH2D hPtEtaProton("hPtEtaProton", "p raw p_{T}-#eta counts;"
                                    "p_{T} (GeV/#it{c});#eta",
                    nPtBins, ptEdges.data(), nEtaBins, etaEdges.data());

  TH2D hPtCentCharged("hPtCentCharged",
                      "Charged raw p_{T}-centrality counts;"
                      "p_{T} (GeV/#it{c});Centrality (%)",
                      nPtBins, ptEdges.data(), nCentralityBins,
                      centralityEdges.data());
  TH2D hPtCentPion("hPtCentPion", "#pi raw p_{T}-centrality counts;"
                                  "p_{T} (GeV/#it{c});Centrality (%)",
                   nPtBins, ptEdges.data(), nCentralityBins,
                   centralityEdges.data());
  TH2D hPtCentKaon("hPtCentKaon", "K raw p_{T}-centrality counts;"
                                  "p_{T} (GeV/#it{c});Centrality (%)",
                   nPtBins, ptEdges.data(), nCentralityBins,
                   centralityEdges.data());
  TH2D hPtCentProton("hPtCentProton", "p raw p_{T}-centrality counts;"
                                      "p_{T} (GeV/#it{c});Centrality (%)",
                     nPtBins, ptEdges.data(), nCentralityBins,
                     centralityEdges.data());

  TH3D hPtEtaCentCharged("hPtEtaCentCharged",
                         "Charged raw p_{T}-#eta-centrality counts;"
                         "p_{T} (GeV/#it{c});#eta;Centrality (%)",
                         nPtBins, ptEdges.data(), nEtaBins, etaEdges.data(),
                         nCentralityBins, centralityEdges.data());
  TH3D hPtEtaCentPion("hPtEtaCentPion",
                      "#pi raw p_{T}-#eta-centrality counts;"
                      "p_{T} (GeV/#it{c});#eta;Centrality (%)",
                      nPtBins, ptEdges.data(), nEtaBins, etaEdges.data(),
                      nCentralityBins, centralityEdges.data());
  TH3D hPtEtaCentKaon("hPtEtaCentKaon",
                      "K raw p_{T}-#eta-centrality counts;"
                      "p_{T} (GeV/#it{c});#eta;Centrality (%)",
                      nPtBins, ptEdges.data(), nEtaBins, etaEdges.data(),
                      nCentralityBins, centralityEdges.data());
  TH3D hPtEtaCentProton("hPtEtaCentProton",
                        "p raw p_{T}-#eta-centrality counts;"
                        "p_{T} (GeV/#it{c});#eta;Centrality (%)",
                        nPtBins, ptEdges.data(), nEtaBins, etaEdges.data(),
                        nCentralityBins, centralityEdges.data());

  const Long64_t processedEvents = forEachConfiguredEvent(
      inputConfigFile, maxFilesPerConfig, maxConfigs,
      [&](const Event &event, const CentralityConfig &)
      {
        const double centrality =
            centralityFromImpactParameter(event.imp, config);
        hCent.Fill(centrality);
        hEventCount.Fill(0.5);

        for (const auto &track : event.particles)
        {
          const double eta = track.GetEta();
          const bool etaOk =
              config.strictEtaBounds
                  ? eta > config.ptSpectraOutput.etaMin &&
                        eta < config.ptSpectraOutput.etaMax
                  : eta >= config.ptSpectraOutput.etaMin &&
                        eta <= config.ptSpectraOutput.etaMax;
          if (!etaOk || !isChargedPdg(track.pdgPid))
          {
            continue;
          }

          const double pt = track.GetPt();
          const int absPdg = std::abs(track.pdgPid);

          if (config.chargedPt.contains(pt, config.strictPtBounds))
          {
            hPtCharged.Fill(pt);
            hPhi.Fill(phiInZeroTwoPi(track.GetPhi()));
            hEta.Fill(eta);
            hPtEtaCharged.Fill(pt, eta);
            hPtCentCharged.Fill(pt, centrality);
            hPtEtaCentCharged.Fill(pt, eta, centrality);
          }

          const SpeciesDefinition *species = findSpecies(absPdg);
          if (!species ||
              !ptRangeForPdg(absPdg, config)
                   .contains(pt, config.strictPtBounds))
          {
            continue;
          }

          switch (species->species)
          {
          case Species::Pion:
            hPtPion.Fill(pt);
            hPtEtaPion.Fill(pt, eta);
            hPtCentPion.Fill(pt, centrality);
            hPtEtaCentPion.Fill(pt, eta, centrality);
            break;
          case Species::Kaon:
            hPtKaon.Fill(pt);
            hPtEtaKaon.Fill(pt, eta);
            hPtCentKaon.Fill(pt, centrality);
            hPtEtaCentKaon.Fill(pt, eta, centrality);
            break;
          case Species::Proton:
            hPtProton.Fill(pt);
            hPtEtaProton.Fill(pt, eta);
            hPtCentProton.Fill(pt, centrality);
            hPtEtaCentProton.Fill(pt, eta, centrality);
            break;
          }
        }
      });

  TFile outputFileHandle(outputFile, "RECREATE");
  if (outputFileHandle.IsZombie())
  {
    throw std::runtime_error(std::string("Cannot create output file: ") +
                             outputFile);
  }

  TDirectory *taskDirectory = outputFileHandle.mkdir(kTaskDirectory);
  writeObject(*taskDirectory, hEventCount);
  writeObject(*taskDirectory, hPtCharged);
  writeObject(*taskDirectory, hPtPion);
  writeObject(*taskDirectory, hPtKaon);
  writeObject(*taskDirectory, hPtProton);
  writeObject(*taskDirectory, hPtEtaCharged);
  writeObject(*taskDirectory, hPtEtaPion);
  writeObject(*taskDirectory, hPtEtaKaon);
  writeObject(*taskDirectory, hPtEtaProton);
  writeObject(*taskDirectory, hPtCentCharged);
  writeObject(*taskDirectory, hPtCentPion);
  writeObject(*taskDirectory, hPtCentKaon);
  writeObject(*taskDirectory, hPtCentProton);
  writeObject(*taskDirectory, hPtEtaCentCharged);
  writeObject(*taskDirectory, hPtEtaCentPion);
  writeObject(*taskDirectory, hPtEtaCentKaon);
  writeObject(*taskDirectory, hPtEtaCentProton);
  writeObject(*taskDirectory, hCent);
  writeObject(*taskDirectory, hPhi);
  writeObject(*taskDirectory, hEta);

  outputFileHandle.Close();
  std::cout << "Processed " << processedEvents << " AMPT events; wrote "
            << outputFile << "." << std::endl;
}

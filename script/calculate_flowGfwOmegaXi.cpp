#include "../analysisConfig.h"
#include "../analysisUtils.h"

#include "PWGCF/GenericFramework/Core/GFW.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "THnSparse.h"
#include "TAxis.h"
#include "TMath.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace ampt_gfw_omega_xi_macro
{
  using ampt_analysis::AnalysisConfig;
  using ampt_analysis::centralityFromImpactParameter;
  using ampt_analysis::isChargedPdg;

  // Particle masses (GeV/c^2), matching the O2 physics constants used by
  // flowGfwOmegaXi. AMPT truth has no detector PID, so these are the only
  // masses needed for the invariant-mass reconstruction.
  constexpr double kMassPi = 0.1395704;
  constexpr double kMassK = 0.493677;
  constexpr double kMassPr = 0.9382721;
  constexpr double kMassK0Short = 0.497611;
  constexpr double kMassLambda = 1.115683;
  constexpr double kMassXi = 1.32171;
  constexpr double kMassOmega = 1.67245;

  // O2 cfgaxisPt (charged pT, 38 bins).
  const double kChargedPtEdges[] = {0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55,
                                    0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95,
                                    1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70,
                                    1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00,
                                    3.50, 4.00, 4.50, 5.00, 5.50, 6.00, 10.0};
  const int kNChargedPtBins = static_cast<int>(sizeof(kChargedPtEdges) / sizeof(kChargedPtEdges[0])) - 1;

  // O2 cfgaxisPtXi/cfgaxisPtOmega/cfgaxisPtK0s/cfgaxisPtLambda (14 bins).
  const double kSpeciesPtEdges[] = {0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3,
                                    2.5, 2.7, 2.9, 3.9, 4.9, 5.9, 9.9};
  const int kNSpeciesPtBins = static_cast<int>(sizeof(kSpeciesPtEdges) / sizeof(kSpeciesPtEdges[0])) - 1;

  // O2 axisMultiplicity (centrality, 10 bins).
  const double kCentralityEdges[] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90};
  const int kNCentralityBins = static_cast<int>(sizeof(kCentralityEdges) / sizeof(kCentralityEdges[0])) - 1;

  // Number of K0s, Lambda, Xi, Omega bins used by the internal GFW mass axes.
  // Order: {K0s, Lambda, Xi, Omega}; copied from O2 cfgMassBins.
  const int kMassBins[] = {80, 32, 14, 16};

  constexpr double kTwoPi = 2.0 * TMath::Pi();
  constexpr double kCutPtMin = 0.2;
  constexpr double kCutPtMax = 10.0;
  constexpr double kCutEta = 0.8;
  constexpr double kCascRapidity = 0.5;

  constexpr double kV0K0sMassWindow = 0.1;
  constexpr double kV0LambdaMassWindow = 0.04;
  constexpr double kV0CompRejLambda = 0.01;
  constexpr double kV0CompRejK0s = 0.005;
  constexpr double kCascLambdaMassWindow = 0.04;
  constexpr double kCascCompMassRej = 0.008;

  inline double totalMomentum(double px, double py, double pz)
  {
    return std::sqrt(px * px + py * py + pz * pz);
  }

  inline double invariantMass2(double m1, double px1, double py1, double pz1,
                               double m2, double px2, double py2, double pz2)
  {
    const double p1sq = px1 * px1 + py1 * py1 + pz1 * pz1;
    const double p2sq = px2 * px2 + py2 * py2 + pz2 * pz2;
    const double e1 = std::sqrt(m1 * m1 + p1sq);
    const double e2 = std::sqrt(m2 * m2 + p2sq);
    const double dot = px1 * px2 + py1 * py2 + pz1 * pz2;
    return m1 * m1 + m2 * m2 + 2.0 * (e1 * e2 - dot);
  }

  inline double invariantMass(double m1, const Track &a, double m2, const Track &b)
  {
    const double m2val = invariantMass2(m1, a.p_x, a.p_y, a.p_z,
                                        m2, b.p_x, b.p_y, b.p_z);
    return std::sqrt(m2val > 0. ? m2val : 0.);
  }

  inline double threeBodyMass(double m1, const Track &a, double m2, const Track &b,
                              double m3, const Track &c)
  {
    const double px = a.p_x + b.p_x + c.p_x;
    const double py = a.p_y + b.p_y + c.p_y;
    const double pz = a.p_z + b.p_z + c.p_z;
    const double e = std::sqrt(m1 * m1 + a.p_x * a.p_x + a.p_y * a.p_y + a.p_z * a.p_z) +
                     std::sqrt(m2 * m2 + b.p_x * b.p_x + b.p_y * b.p_y + b.p_z * b.p_z) +
                     std::sqrt(m3 * m3 + c.p_x * c.p_x + c.p_y * c.p_y + c.p_z * c.p_z);
    const double m2val = e * e - (px * px + py * py + pz * pz);
    return std::sqrt(m2val > 0. ? m2val : 0.);
  }

  inline double pseudoRapidity(double px, double py, double pz)
  {
    const double p = totalMomentum(px, py, pz);
    if (p == std::abs(pz)) {
      return std::copysign(std::numeric_limits<double>::infinity(), pz);
    }
    return 0.5 * std::log((p + pz) / (p - pz));
  }

  inline double rapidity(double mass, double px, double py, double pz)
  {
    const double p = totalMomentum(px, py, pz);
    const double e = std::sqrt(mass * mass + p * p);
    if (e == std::abs(pz)) {
      return std::copysign(std::numeric_limits<double>::infinity(), pz);
    }
    return 0.5 * std::log((e + pz) / (e - pz));
  }

  inline std::vector<double> uniformEdges(int nBins, double min, double max)
  {
    std::vector<double> edges(static_cast<std::size_t>(nBins) + 1);
    const double width = (max - min) / static_cast<double>(nBins);
    for (int i = 0; i <= nBins; ++i) {
      edges[static_cast<std::size_t>(i)] = min + static_cast<double>(i) * width;
    }
    return edges;
  }

  void writeObject(TDirectory &directory, TObject &object)
  {
    directory.WriteTObject(&object, object.GetName());
  }

  void fillProfile(GFW &gfw, const GFW::CorrConfig &corrconf,
                   TProfile &profile, double cent)
  {
    const double dnx = std::real(gfw.Calculate(corrconf, 0, true));
    if (dnx == 0.) {
      return;
    }
    if (!corrconf.pTDif) {
      const double val = std::real(gfw.Calculate(corrconf, 0, false)) / dnx;
      if (std::fabs(val) < 1.) {
        profile.Fill(cent, val, dnx);
      }
    }
  }

  void fillProfilepT(GFW &gfw, const GFW::CorrConfig &corrconf,
                     TProfile2D &profile, int ptbin, TAxis *ptAxis, double cent)
  {
    const double dnx = std::real(gfw.Calculate(corrconf, ptbin - 1, true));
    if (dnx == 0.) {
      return;
    }
    const double val = std::real(gfw.Calculate(corrconf, ptbin - 1, false)) / dnx;
    if (std::fabs(val) < 1.) {
      profile.Fill(ptAxis->GetBinCenter(ptbin), cent, val, dnx);
    }
  }

  void fillProfilepTMass(GFW &gfw, const GFW::CorrConfig &corrconf,
                         TProfile3D &profile, int ptbin, TAxis *ptAxis,
                         TAxis *massAxis, int nMassBins, int nPtBins, double cent)
  {
    for (int massbin = 1; massbin <= nMassBins; ++massbin) {
      const int ptIndex = (ptbin - 1) + (massbin - 1) * nPtBins;
      const double dnx = std::real(gfw.Calculate(corrconf, ptIndex, true));
      if (dnx == 0.) {
        continue;
      }
      const double val = std::real(gfw.Calculate(corrconf, ptIndex, false)) / dnx;
      if (std::fabs(val) < 1.) {
        profile.Fill(ptAxis->GetBinCenter(ptbin), massAxis->GetBinCenter(massbin),
                     cent, val, dnx);
      }
    }
  }

  struct LambdaCandidate
  {
    int charge = 0; // +1 Lambda (p + pi-), -1 antiLambda (pbar + pi+)
    double px = 0.;
    double py = 0.;
    double pz = 0.;
    double pt = 0.;
    double eta = 0.;
    double phi = 0.;
    double mLambda = 0.;
    double mK0sRej = 0.;
    int protonIndex = -1;
    int pionIndex = -1;
  };

  struct K0sCandidate
  {
    double px = 0.;
    double py = 0.;
    double pz = 0.;
    double pt = 0.;
    double eta = 0.;
    double phi = 0.;
    double mK0s = 0.;
    double mLambdaRej = 0.;
    int posPiIndex = -1;
    int negPiIndex = -1;
  };

  TH1D makeEventCountHistogram()
  {
    TH1D eventCount("hEventCount", "", 14, 0., 14.);
    eventCount.GetXaxis()->SetBinLabel(1, "Filtered event");
    eventCount.GetXaxis()->SetBinLabel(2, "after sel8");
    eventCount.GetXaxis()->SetBinLabel(3, "after kTVXinTRD");
    eventCount.GetXaxis()->SetBinLabel(4, "after kNoTimeFrameBorder");
    eventCount.GetXaxis()->SetBinLabel(5, "after kNoITSROFrameBorder");
    eventCount.GetXaxis()->SetBinLabel(6, "after kDoNoSameBunchPileup");
    eventCount.GetXaxis()->SetBinLabel(7, "after kIsGoodZvtxFT0vsPV");
    eventCount.GetXaxis()->SetBinLabel(8, "after kNoCollInTimeRangeStandard");
    eventCount.GetXaxis()->SetBinLabel(9, "after kIsGoodITSLayersAll");
    eventCount.GetXaxis()->SetBinLabel(10, "after MultPVCut");
    eventCount.GetXaxis()->SetBinLabel(11, "after TPC occupancy cut");
    eventCount.GetXaxis()->SetBinLabel(12, "after V0AT0Acut");
    eventCount.GetXaxis()->SetBinLabel(13, "after IRmincut");
    eventCount.GetXaxis()->SetBinLabel(14, "after IRmaxcut");
    return eventCount;
  }

  std::unique_ptr<THnSparseF> makeInvMass(const char *name, int massBins,
                                          double massMin, double massMax)
  {
    std::vector<TAxis> axes;
    axes.reserve(4);
    axes.emplace_back(TAxis(kNSpeciesPtBins, kSpeciesPtEdges));
    axes.emplace_back(TAxis(massBins, massMin, massMax));
    axes.emplace_back(TAxis(40, -1., 1.));
    axes.emplace_back(TAxis(kNCentralityBins, kCentralityEdges));
    auto hist = std::make_unique<THnSparseF>(name, "", axes);
    return hist;
  }
} // namespace ampt_gfw_omega_xi_macro

/**
 * Produce the processData-compatible flowGfwOmegaXi output from AMPT truth.
 *
 * AMPT truth has no detector track-quality, PID, NUA/NUE, CCDB or local-density
 * inputs, so the corresponding QA and correction histograms are not produced.
 * K0s, Lambda/antiLambda, Xi/antiXi and Omega/antiOmega are reconstructed from
 * daughter PDG codes plus invariant-mass, rapidity and competing-mass cuts.
 */
void calculate_flowGfwOmegaXi(
    const char *inputConfigFile = "../config/cent_cfg.json",
    const char *outputFile = "myAnalysisResultFlowGfwOmegaXi.root",
    int maxFilesPerConfig = -1, int maxConfigs = -1,
    const char *analysisConfigFile = "../config/config.json",
    int shardIndex = 0, int shardCount = 1)
{
  using namespace ampt_analysis;
  using namespace ampt_gfw_omega_xi_macro;

  AnalysisConfig config = loadAnalysisConfig(analysisConfigFile);

  // Internal GFW axes. These match fPtAxis/fXiPtAxis/... in flowGfwOmegaXi.
  std::unique_ptr<TAxis> fPtAxis(new TAxis(kNChargedPtBins, kChargedPtEdges));
  std::unique_ptr<TAxis> fXiPtAxis(new TAxis(kNSpeciesPtBins, kSpeciesPtEdges));
  std::unique_ptr<TAxis> fOmegaPtAxis(new TAxis(kNSpeciesPtBins, kSpeciesPtEdges));
  std::unique_ptr<TAxis> fK0sPtAxis(new TAxis(kNSpeciesPtBins, kSpeciesPtEdges));
  std::unique_ptr<TAxis> fLambdaPtAxis(new TAxis(kNSpeciesPtBins, kSpeciesPtEdges));
  std::unique_ptr<TAxis> fXiMass(new TAxis(kMassBins[2], 1.29, 1.36));
  std::unique_ptr<TAxis> fOmegaMass(new TAxis(kMassBins[3], 1.63, 1.71));
  std::unique_ptr<TAxis> fK0sMass(new TAxis(kMassBins[0], 0.4, 0.6));
  std::unique_ptr<TAxis> fLambdaMass(new TAxis(kMassBins[1], 1.08, 1.16));

  const int nXiPtBins = kNSpeciesPtBins;
  const int nOmegaPtBins = kNSpeciesPtBins;
  const int nK0sPtBins = kNSpeciesPtBins;
  const int nLambdaPtBins = kNSpeciesPtBins;

  const int nXiPtMassBins = nXiPtBins * kMassBins[2];
  const int nOmegaPtMassBins = nXiPtBins * kMassBins[3];
  const int nK0sPtMassBins = nK0sPtBins * kMassBins[0];
  const int nLambdaPtMassBins = nLambdaPtBins * kMassBins[1];

  // GFW regions. Strings and order are copied from flowGfwOmegaXi::init().
  GFW gfw;
  gfw.AddRegion("reffull", -0.8, 0.8, 1, 1);
  gfw.AddRegion("refN10", -0.8, -0.4, 1, 1);
  gfw.AddRegion("refP10", 0.4, 0.8, 1, 1);
  gfw.AddRegion("poiN10dpt", -0.8, -0.4, kNChargedPtBins, 32);
  gfw.AddRegion("poiP10dpt", 0.4, 0.8, kNChargedPtBins, 32);
  gfw.AddRegion("poifulldpt", -0.8, 0.8, kNChargedPtBins, 32);
  gfw.AddRegion("poioldpt", -0.8, 0.8, kNChargedPtBins, 1);
  gfw.AddRegion("poiXiPdpt", 0.4, 0.8, nXiPtMassBins, 2);
  gfw.AddRegion("poiXiNdpt", -0.8, -0.4, nXiPtMassBins, 2);
  gfw.AddRegion("poiXifulldpt", -0.8, 0.8, nXiPtMassBins, 2);
  gfw.AddRegion("poiOmegaPdpt", 0.4, 0.8, nOmegaPtMassBins, 4);
  gfw.AddRegion("poiOmegaNdpt", -0.8, -0.4, nOmegaPtMassBins, 4);
  gfw.AddRegion("poiOmegafulldpt", -0.8, 0.8, nOmegaPtMassBins, 4);
  gfw.AddRegion("poiK0sPdpt", 0.4, 0.8, nK0sPtMassBins, 8);
  gfw.AddRegion("poiK0sNdpt", -0.8, -0.4, nK0sPtMassBins, 8);
  gfw.AddRegion("poiK0sfulldpt", -0.8, 0.8, nK0sPtMassBins, 8);
  gfw.AddRegion("poiLambdaPdpt", 0.4, 0.8, nLambdaPtMassBins, 16);
  gfw.AddRegion("poiLambdaNdpt", -0.8, -0.4, nLambdaPtMassBins, 16);
  gfw.AddRegion("poiLambdafulldpt", -0.8, 0.8, nLambdaPtMassBins, 16);
  // MC regions (only needed so that the full 45-entry corrconfig vector below
  // can be parsed; the MC configs are never calculated in processData).
  gfw.AddRegion("refN10MC", -0.8, -0.4, 1, 64);
  gfw.AddRegion("refP10MC", 0.4, 0.8, 1, 64);
  gfw.AddRegion("poiXiPdptMC", 0.4, 0.8, nXiPtMassBins, 128);
  gfw.AddRegion("poiXiNdptMC", -0.8, -0.4, nXiPtMassBins, 128);
  gfw.AddRegion("poiOmegaPdptMC", 0.4, 0.8, nOmegaPtMassBins, 256);
  gfw.AddRegion("poiOmegaNdptMC", -0.8, -0.4, nOmegaPtMassBins, 256);
  gfw.AddRegion("poiK0sPdptMC", 0.4, 0.8, nK0sPtMassBins, 512);
  gfw.AddRegion("poiK0sNdptMC", -0.8, -0.4, nK0sPtMassBins, 512);
  gfw.AddRegion("poiLambdaPdptMC", 0.4, 0.8, nLambdaPtMassBins, 1024);
  gfw.AddRegion("poiLambdaNdptMC", -0.8, -0.4, nLambdaPtMassBins, 1024);

  std::vector<GFW::CorrConfig> corrconfigs;
  corrconfigs.reserve(45);
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiP10dpt {2} refN10 {-2}", "Poi10Gap22dpta", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiN10dpt {2} refP10 {-2}", "Poi10Gap22dptb", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poifulldpt reffull | poioldpt {2 2 -2 -2}", "Poi10Gap24dpt", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poifulldpt reffull | poioldpt {2 -2}", "PoiFull22dpt", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiXiPdpt {2} refN10 {-2}", "Xi10Gap22a", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiXiNdpt {2} refP10 {-2}", "Xi10Gap22b", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiXifulldpt reffull {2 2 -2 -2}", "Xi10Gap24", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiXifulldpt {2} reffull {-2}", "XiFull22", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiOmegaPdpt {2} refN10 {-2}", "Omega10Gap22a", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiOmegaNdpt {2} refP10 {-2}", "Omega10Gap22b", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiOmegafulldpt reffull {2 2 -2 -2}", "Omega10Gap24", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiOmegafulldpt {2} reffull {-2}", "OmegaFull22", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiK0sPdpt {2} refN10 {-2}", "K0short10Gap22a", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiK0sNdpt {2} refP10 {-2}", "K0short10Gap22b", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiK0sfulldpt reffull {2 2 -2 -2}", "K0short10Gap24", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiK0sfulldpt {2} reffull {-2}", "K0shortFull22", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiLambdaPdpt {2} refN10 {-2}", "Lambda10Gap22a", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiLambdaNdpt {2} refP10 {-2}", "Lambda10Gap22b", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiLambdafulldpt reffull {2 2 -2 -2}", "LambdaFull24", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiLambdafulldpt {2} reffull {-2}", "LambdaFull22", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("refP10 {2} refN10 {-2}", "Ref10Gap22a", false));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("reffull reffull {2 2 -2 -2}", "Ref10Gap24", false));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("reffull reffull {2 -2}", "RefFull22", false));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiXiPdpt {3} refN10 {-3}", "Xi10Gap32a", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiXiNdpt {3} refP10 {-3}", "Xi10Gap32b", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiOmegaPdpt {3} refN10 {-3}", "Omega10Gap32a", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiOmegaNdpt {3} refP10 {-3}", "Omega10Gap32b", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiK0sPdpt {3} refN10 {-3}", "K0short10Gap32a", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiK0sNdpt {3} refP10 {-3}", "K0short10Gap32b", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiLambdaPdpt {3} refN10 {-3}", "Lambda10Gap32a", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiLambdaNdpt {3} refP10 {-3}", "Lambda10Gap32b", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("refP10 {3} refN10 {-3}", "Ref10Gap32a", false));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiXiPdptMC {2} refN10MC {-2}", "MCXi10Gap22a", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiXiNdptMC {2} refP10MC {-2}", "MCXi10Gap22b", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiOmegaPdptMC {2} refN10MC {-2}", "MCOmega10Gap22a", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiOmegaNdptMC {2} refP10MC {-2}", "MCOmega10Gap22b", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiK0sPdptMC {2} refN10MC {-2}", "MCK0s10Gap22a", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiK0sNdptMC {2} refP10MC {-2}", "MCK0s10Gap22b", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiLambdaPdptMC {2} refN10MC {-2}", "MCLambda10Gap22a", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiLambdaNdptMC {2} refP10MC {-2}", "MCLambda10Gap22b", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("refP10MC {2} refN10MC {-2}", "MCRef10Gap22a", false));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiXiPdpt refP10 {2, 2} refN10 {-2 -2}", "Xi10Gap24a", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiXiNdpt refN10 {2, 2} refP10 {-2 -2}", "Xi10Gap24b", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiOmegaPdpt refP10 {2, 2} refN10 {-2 -2}", "Omega10Gap24a", true));
  corrconfigs.push_back(gfw.GetCorrelatorConfig("poiOmegaNdpt refN10 {2, 2} refP10 {-2 -2}", "Omega10Gap24b", true));
  gfw.CreateRegions();

  // Output objects.
  std::unique_ptr<TH1D> hPhi(new TH1D("hPhi", "", 60, 0., kTwoPi));
  std::unique_ptr<TH1D> hPhicorr(new TH1D("hPhicorr", "", 60, 0., kTwoPi));
  std::unique_ptr<TH1D> hPhiK0s(new TH1D("hPhiK0s", "", 60, 0., kTwoPi));
  std::unique_ptr<TH1D> hPhiLambda(new TH1D("hPhiLambda", "", 60, 0., kTwoPi));
  std::unique_ptr<TH1D> hPhiXi(new TH1D("hPhiXi", "", 60, 0., kTwoPi));
  std::unique_ptr<TH1D> hPhiOmega(new TH1D("hPhiOmega", "", 60, 0., kTwoPi));
  std::unique_ptr<TH1D> hPhiK0scorr(new TH1D("hPhiK0scorr", "", 60, 0., kTwoPi));
  std::unique_ptr<TH1D> hPhiLambdacorr(new TH1D("hPhiLambdacorr", "", 60, 0., kTwoPi));
  std::unique_ptr<TH1D> hPhiXicorr(new TH1D("hPhiXicorr", "", 60, 0., kTwoPi));
  std::unique_ptr<TH1D> hPhiOmegacorr(new TH1D("hPhiOmegacorr", "", 60, 0., kTwoPi));
  std::unique_ptr<TH1D> hEta(new TH1D("hEta", "", 40, -1., 1.));
  std::unique_ptr<TH1D> hVtxZ(new TH1D("hVtxZ", "", 20, -10., 10.));
  std::unique_ptr<TH1D> hMult(new TH1D("hMult", "", 3000, 0.5, 3000.5));
  std::unique_ptr<TH1D> hMultTPC(new TH1D("hMultTPC", "", 3000, 0.5, 3000.5));
  std::unique_ptr<TH1D> hCent(new TH1D("hCent", "", 90, 0., 90.));
  std::unique_ptr<TH1D> hPt(new TH1D("hPt", "", kNChargedPtBins, kChargedPtEdges));
  std::unique_ptr<TH2D> hNTracksPVvsCentrality(new TH2D("hNTracksPVvsCentrality", "", 5000, 0., 5000., 100, 0., 100.));
  std::unique_ptr<TH2D> hmultFV0AvsmultFT0A(new TH2D("hmultFV0AvsmultFT0A", "", 5000, 0., 5000., 5000, 0., 5000.));
  TH1D hEventCount = makeEventCountHistogram();
  std::unique_ptr<TH1D> hInteractionRate(new TH1D("hInteractionRate", "", 1000, 0., 1000.));

  // Integrated and pT-differential charged cumulant profiles.
  std::unique_ptr<TProfile> c22(new TProfile("c22", ";Centrality  (%) ; C_{2}{2}", kNCentralityBins, kCentralityEdges));
  std::unique_ptr<TProfile> c32(new TProfile("c32", ";Centrality  (%) ; C_{2}{2}", kNCentralityBins, kCentralityEdges));
  std::unique_ptr<TProfile> c24(new TProfile("c24", ";Centrality  (%) ; C_{2}{2}", kNCentralityBins, kCentralityEdges));
  std::unique_ptr<TProfile> c22Full(new TProfile("c22Full", ";Centrality  (%) ; C_{2}{2}", kNCentralityBins, kCentralityEdges));
  std::unique_ptr<TProfile2D> c22dpt(new TProfile2D("c22dpt", ";Centrality  (%) ; C_{2}{2}", kNChargedPtBins, kChargedPtEdges, kNCentralityBins, kCentralityEdges));
  std::unique_ptr<TProfile2D> c24dpt(new TProfile2D("c24dpt", ";Centrality  (%) ; C_{2}{4}", kNChargedPtBins, kChargedPtEdges, kNCentralityBins, kCentralityEdges));
  std::unique_ptr<TProfile2D> c22Fulldpt(new TProfile2D("c22Fulldpt", ";Centrality  (%) ; C_{2}{2}", kNChargedPtBins, kChargedPtEdges, kNCentralityBins, kCentralityEdges));

  const std::vector<double> xiMassProfileEdges = uniformEdges(14, 1.3, 1.37);
  const std::vector<double> omegaMassProfileEdges = uniformEdges(16, 1.63, 1.71);
  const std::vector<double> k0sMassProfileEdges = uniformEdges(40, 0.4, 0.6);
  const std::vector<double> lambdaMassProfileEdges = uniformEdges(32, 1.08, 1.16);

  auto makeSpeciesProfile3D = [&](const char *name, const std::vector<double> &massEdges) {
    return std::unique_ptr<TProfile3D>(new TProfile3D(
        name, ";pt ; C_{2}{2}", kNSpeciesPtBins, kSpeciesPtEdges,
        static_cast<int>(massEdges.size()) - 1, massEdges.data(),
        kNCentralityBins, kCentralityEdges));
  };

  std::unique_ptr<TProfile3D> Xic22dpt = makeSpeciesProfile3D("Xic22dpt", xiMassProfileEdges);
  std::unique_ptr<TProfile3D> Xic24dpt = makeSpeciesProfile3D("Xic24dpt", xiMassProfileEdges);
  std::unique_ptr<TProfile3D> Xic22Fulldpt = makeSpeciesProfile3D("Xic22Fulldpt", xiMassProfileEdges);
  std::unique_ptr<TProfile3D> Xic24_gapdpt = makeSpeciesProfile3D("Xic24_gapdpt", xiMassProfileEdges);
  std::unique_ptr<TProfile3D> Xic32dpt = makeSpeciesProfile3D("Xic32dpt", xiMassProfileEdges);
  std::unique_ptr<TProfile3D> Omegac22dpt = makeSpeciesProfile3D("Omegac22dpt", omegaMassProfileEdges);
  std::unique_ptr<TProfile3D> Omegac24dpt = makeSpeciesProfile3D("Omegac24dpt", omegaMassProfileEdges);
  std::unique_ptr<TProfile3D> Omegac22Fulldpt = makeSpeciesProfile3D("Omegac22Fulldpt", omegaMassProfileEdges);
  std::unique_ptr<TProfile3D> Omegac24_gapdpt = makeSpeciesProfile3D("Omegac24_gapdpt", omegaMassProfileEdges);
  std::unique_ptr<TProfile3D> Omegac32dpt = makeSpeciesProfile3D("Omegac32dpt", omegaMassProfileEdges);
  std::unique_ptr<TProfile3D> K0sc22dpt = makeSpeciesProfile3D("K0sc22dpt", k0sMassProfileEdges);
  std::unique_ptr<TProfile3D> K0sc24dpt = makeSpeciesProfile3D("K0sc24dpt", k0sMassProfileEdges);
  std::unique_ptr<TProfile3D> K0sc22Fulldpt = makeSpeciesProfile3D("K0sc22Fulldpt", k0sMassProfileEdges);
  std::unique_ptr<TProfile3D> K0sc32dpt = makeSpeciesProfile3D("K0sc32dpt", k0sMassProfileEdges);
  std::unique_ptr<TProfile3D> Lambdac22dpt = makeSpeciesProfile3D("Lambdac22dpt", lambdaMassProfileEdges);
  std::unique_ptr<TProfile3D> Lambdac24dpt = makeSpeciesProfile3D("Lambdac24dpt", lambdaMassProfileEdges);
  std::unique_ptr<TProfile3D> Lambdac22Fulldpt = makeSpeciesProfile3D("Lambdac22Fulldpt", lambdaMassProfileEdges);
  std::unique_ptr<TProfile3D> Lambdac32dpt = makeSpeciesProfile3D("Lambdac32dpt", lambdaMassProfileEdges);

  // Invariant-mass sparse histograms. The mass axes intentionally differ from
  // the GFW/internal mass axes, matching the O2 HistogramRegistry layout.
  std::unique_ptr<THnSparseF> InvMassXi_all = makeInvMass("InvMassXi_all", 80, 1.29, 1.37);
  std::unique_ptr<THnSparseF> InvMassXi = makeInvMass("InvMassXi", 80, 1.29, 1.37);
  std::unique_ptr<THnSparseF> InvMassOmega_all = makeInvMass("InvMassOmega_all", 80, 1.63, 1.71);
  std::unique_ptr<THnSparseF> InvMassOmega = makeInvMass("InvMassOmega", 80, 1.63, 1.71);
  std::unique_ptr<THnSparseF> InvMassK0s_all = makeInvMass("InvMassK0s_all", 400, 0.4, 0.6);
  std::unique_ptr<THnSparseF> InvMassK0s = makeInvMass("InvMassK0s", 400, 0.4, 0.6);
  std::unique_ptr<THnSparseF> InvMassLambda_all = makeInvMass("InvMassLambda_all", 160, 1.08, 1.16);
  std::unique_ptr<THnSparseF> InvMassLambda = makeInvMass("InvMassLambda", 160, 1.08, 1.16);
  std::unique_ptr<THnSparseF> InvMassALambda_all = makeInvMass("InvMassALambda_all", 160, 1.08, 1.16);
  std::unique_ptr<THnSparseF> InvMassALambda = makeInvMass("InvMassALambda", 160, 1.08, 1.16);

  const Long64_t processedEvents = forEachConfiguredEvent(
      inputConfigFile, maxFilesPerConfig, maxConfigs, shardIndex, shardCount,
      [&](const Event &event, const CentralityConfig &,
          const EventIdentity &)
      {
        hEventCount.Fill(0.5);
        if (event.particles.size() < 1) {
          return;
        }
        gfw.Clear();

        const double cent = centralityFromImpactParameter(event.imp, config);
        hEventCount.Fill(1.5);
        hCent->Fill(cent);

        // Index charged tracks used for daughters and the charged flow.
        std::vector<int> posPi;
        std::vector<int> negPi;
        std::vector<int> proton;
        std::vector<int> antiproton;
        std::vector<int> posKaon;
        std::vector<int> negKaon;
        posPi.reserve(event.particles.size());
        negPi.reserve(event.particles.size());
        proton.reserve(event.particles.size());
        antiproton.reserve(event.particles.size());
        posKaon.reserve(event.particles.size());
        negKaon.reserve(event.particles.size());

        int nCharged = 0;
        for (int i = 0; i < static_cast<int>(event.particles.size()); ++i) {
          const Track &track = event.particles[i];
          if (!isChargedPdg(track.pdgPid)) {
            continue;
          }
          const double pt = track.GetPt();
          const double eta = track.GetEta();
          if (pt <= kCutPtMin || pt >= kCutPtMax || std::fabs(eta) >= kCutEta) {
            continue;
          }

          ++nCharged;
          hPhi->Fill(track.GetPhi());
          hPhicorr->Fill(track.GetPhi(), 1.0);
          hEta->Fill(eta);
          hPt->Fill(pt);
          const int ptbin = fPtAxis->FindBin(pt) - 1;
          gfw.Fill(eta, ptbin, track.GetPhi(), 1.0, 1);
          gfw.Fill(eta, ptbin, track.GetPhi(), 1.0, 32);

          switch (track.pdgPid) {
          case 211:
            posPi.push_back(i);
            break;
          case -211:
            negPi.push_back(i);
            break;
          case 2212:
            proton.push_back(i);
            break;
          case -2212:
            antiproton.push_back(i);
            break;
          case 321:
            posKaon.push_back(i);
            break;
          case -321:
            negKaon.push_back(i);
            break;
          default:
            break;
          }
        }
        hMult->Fill(nCharged);

        // ---- V0 reconstruction -------------------------------------------
        std::vector<K0sCandidate> k0sCandidates;
        std::vector<LambdaCandidate> lambdaCandidates;

        // K0s: pi+ pi-
        for (int pi : posPi) {
          const Track &a = event.particles[pi];
          for (int nj : negPi) {
            const Track &b = event.particles[nj];
            const double px = a.p_x + b.p_x;
            const double py = a.p_y + b.p_y;
            const double pz = a.p_z + b.p_z;
            const double pt = std::sqrt(px * px + py * py);
            if (pt <= kCutPtMin || pt >= kCutPtMax) {
              continue;
            }
            const double mK0s = invariantMass(kMassPi, a, kMassPi, b);
            if (std::fabs(mK0s - kMassK0Short) >= kV0K0sMassWindow) {
              continue;
            }
            K0sCandidate cand;
            cand.px = px;
            cand.py = py;
            cand.pz = pz;
            cand.pt = pt;
            cand.eta = pseudoRapidity(px, py, pz);
            cand.phi = std::atan2(py, px);
            cand.mK0s = mK0s;
            cand.mLambdaRej = invariantMass(kMassPr, a, kMassPi, b);
            cand.posPiIndex = pi;
            cand.negPiIndex = nj;
            k0sCandidates.push_back(cand);
          }
        }

        // Lambda: p + pi- ; antiLambda: pbar + pi+
        for (int p : proton) {
          const Track &prot = event.particles[p];
          for (int pi : negPi) {
            const Track &pion = event.particles[pi];
            const double px = prot.p_x + pion.p_x;
            const double py = prot.p_y + pion.p_y;
            const double pz = prot.p_z + pion.p_z;
            const double pt = std::sqrt(px * px + py * py);
            if (pt <= kCutPtMin || pt >= kCutPtMax) {
              continue;
            }
            const double mLambda = invariantMass(kMassPr, prot, kMassPi, pion);
            if (std::fabs(mLambda - kMassLambda) >= kV0LambdaMassWindow) {
              continue;
            }
            LambdaCandidate cand;
            cand.charge = 1;
            cand.px = px;
            cand.py = py;
            cand.pz = pz;
            cand.pt = pt;
            cand.eta = pseudoRapidity(px, py, pz);
            cand.phi = std::atan2(py, px);
            cand.mLambda = mLambda;
            cand.mK0sRej = invariantMass(kMassPi, prot, kMassPi, pion);
            cand.protonIndex = p;
            cand.pionIndex = pi;
            lambdaCandidates.push_back(cand);
          }
        }
        for (int pbar : antiproton) {
          const Track &prot = event.particles[pbar];
          for (int pi : posPi) {
            const Track &pion = event.particles[pi];
            const double px = prot.p_x + pion.p_x;
            const double py = prot.p_y + pion.p_y;
            const double pz = prot.p_z + pion.p_z;
            const double pt = std::sqrt(px * px + py * py);
            if (pt <= kCutPtMin || pt >= kCutPtMax) {
              continue;
            }
            const double mLambda = invariantMass(kMassPr, prot, kMassPi, pion);
            if (std::fabs(mLambda - kMassLambda) >= kV0LambdaMassWindow) {
              continue;
            }
            LambdaCandidate cand;
            cand.charge = -1;
            cand.px = px;
            cand.py = py;
            cand.pz = pz;
            cand.pt = pt;
            cand.eta = pseudoRapidity(px, py, pz);
            cand.phi = std::atan2(py, px);
            cand.mLambda = mLambda;
            cand.mK0sRej = invariantMass(kMassPi, prot, kMassPi, pion);
            cand.protonIndex = pbar;
            cand.pionIndex = pi;
            lambdaCandidates.push_back(cand);
          }
        }

        // V0 _all histograms and final acceptance.
        for (const K0sCandidate &cand : k0sCandidates) {
          double coords[4] = {cand.pt, cand.mK0s, cand.eta, cent};
          InvMassK0s_all->Fill(coords, 1.0);
        }
        for (const LambdaCandidate &cand : lambdaCandidates) {
          double coords[4] = {cand.pt, cand.mLambda, cand.eta, cent};
          if (cand.charge > 0) {
            InvMassLambda_all->Fill(coords, 1.0);
          } else {
            InvMassALambda_all->Fill(coords, 1.0);
          }
        }

        for (const K0sCandidate &cand : k0sCandidates) {
          bool isK0s = true;
          if (std::fabs(cand.mLambdaRej - kMassLambda) < kV0CompRejLambda) {
            isK0s = false;
          }
          if (!isK0s) {
            continue;
          }
          double coords[4] = {cand.pt, cand.mK0s, cand.eta, cent};
          InvMassK0s->Fill(coords, 1.0);
          hPhiK0s->Fill(cand.phi);
          hPhiK0scorr->Fill(cand.phi, 1.0);
          gfw.Fill(cand.eta,
                   fK0sPtAxis->FindBin(cand.pt) - 1 +
                       (fK0sMass->FindBin(cand.mK0s) - 1) * nK0sPtBins,
                   cand.phi, 1.0, 8);
          if (fK0sPtAxis->FindBin(cand.pt) - 1 == 6) {
            gfw.Fill(cand.eta, fK0sMass->FindBin(cand.mK0s) - 1, cand.phi, 1.0, 2048);
          }
        }

        for (const LambdaCandidate &cand : lambdaCandidates) {
          bool isLambda = true;
          if (std::fabs(cand.mK0sRej - kMassK0Short) < kV0CompRejK0s) {
            isLambda = false;
          }
          if (!isLambda) {
            continue;
          }
          hPhiLambda->Fill(cand.phi);
          hPhiLambdacorr->Fill(cand.phi, 1.0);
          double coords[4] = {cand.pt, cand.mLambda, cand.eta, cent};
          if (cand.charge > 0) {
            InvMassLambda->Fill(coords, 1.0);
          } else {
            InvMassALambda->Fill(coords, 1.0);
          }
          gfw.Fill(cand.eta,
                   fLambdaPtAxis->FindBin(cand.pt) - 1 +
                       (fLambdaMass->FindBin(cand.mLambda) - 1) * nLambdaPtBins,
                   cand.phi, 1.0, 16);
        }

        // ---- Cascade reconstruction --------------------------------------
        // Xi- : pi- bachelor + Lambda ; antiXi : pi+ bachelor + antiLambda
        for (const LambdaCandidate &cand : lambdaCandidates) {
          const std::vector<int> &bachelors = cand.charge > 0 ? negPi : posPi;
          for (int bi : bachelors) {
            if (bi == cand.protonIndex || bi == cand.pionIndex) {
              continue;
            }
            const Track &bach = event.particles[bi];
            const double px = cand.px + bach.p_x;
            const double py = cand.py + bach.p_y;
            const double pz = cand.pz + bach.p_z;
            const double pt = std::sqrt(px * px + py * py);
            if (pt <= kCutPtMin || pt >= kCutPtMax) {
              continue;
            }
            const Track &prot = event.particles[cand.protonIndex];
            const Track &pion = event.particles[cand.pionIndex];
            const double mXi = threeBodyMass(kMassPi, bach, kMassPr, prot, kMassPi, pion);
            const double yXi = rapidity(kMassXi, px, py, pz);
            if (std::fabs(yXi) >= kCascRapidity) {
              continue;
            }
            const double eta = pseudoRapidity(px, py, pz);
            const double phi = std::atan2(py, px);
            double coords[4] = {pt, mXi, eta, cent};
            InvMassXi_all->Fill(coords, 1.0);

            const double mOmega = threeBodyMass(kMassK, bach, kMassPr, prot, kMassPi, pion);
            bool isXi = true;
            if (std::fabs(mOmega - kMassOmega) < kCascCompMassRej) {
              isXi = false;
            }
            if (!isXi) {
              continue;
            }
            InvMassXi->Fill(coords, 1.0);
            hPhiXi->Fill(phi);
            hPhiXicorr->Fill(phi, 1.0);
            gfw.Fill(eta,
                     fXiPtAxis->FindBin(pt) - 1 +
                         (fXiMass->FindBin(mXi) - 1) * nXiPtBins,
                     phi, 1.0, 2);
          }
        }

        // Omega- : K- bachelor + Lambda ; antiOmega : K+ bachelor + antiLambda
        for (const LambdaCandidate &cand : lambdaCandidates) {
          const std::vector<int> &bachelors = cand.charge > 0 ? negKaon : posKaon;
          for (int bi : bachelors) {
            if (bi == cand.protonIndex || bi == cand.pionIndex) {
              continue;
            }
            const Track &bach = event.particles[bi];
            const double px = cand.px + bach.p_x;
            const double py = cand.py + bach.p_y;
            const double pz = cand.pz + bach.p_z;
            const double pt = std::sqrt(px * px + py * py);
            if (pt <= kCutPtMin || pt >= kCutPtMax) {
              continue;
            }
            const Track &prot = event.particles[cand.protonIndex];
            const Track &pion = event.particles[cand.pionIndex];
            const double mOmega = threeBodyMass(kMassK, bach, kMassPr, prot, kMassPi, pion);
            const double yOmega = rapidity(kMassOmega, px, py, pz);
            if (std::fabs(yOmega) >= kCascRapidity) {
              continue;
            }
            const double eta = pseudoRapidity(px, py, pz);
            const double phi = std::atan2(py, px);
            double coords[4] = {pt, mOmega, eta, cent};
            InvMassOmega_all->Fill(coords, 1.0);

            const double mXi = threeBodyMass(kMassPi, bach, kMassPr, prot, kMassPi, pion);
            bool isOmega = true;
            if (std::fabs(mXi - kMassXi) < kCascCompMassRej) {
              isOmega = false;
            }
            if (!isOmega) {
              continue;
            }
            InvMassOmega->Fill(coords, 1.0);
            hPhiOmega->Fill(phi);
            hPhiOmegacorr->Fill(phi, 1.0);
            gfw.Fill(eta,
                     fOmegaPtAxis->FindBin(pt) - 1 +
                         (fOmegaMass->FindBin(mOmega) - 1) * nOmegaPtBins,
                     phi, 1.0, 4);
          }
        }

        // ---- Cumulant profiles -------------------------------------------
        fillProfile(gfw, corrconfigs.at(20), *c22, cent);
        fillProfile(gfw, corrconfigs.at(21), *c24, cent);
        fillProfile(gfw, corrconfigs.at(22), *c22Full, cent);
        fillProfile(gfw, corrconfigs.at(31), *c32, cent);

        for (int i = 1; i <= kNChargedPtBins; ++i) {
          fillProfilepT(gfw, corrconfigs.at(0), *c22dpt, i, fPtAxis.get(), cent);
          fillProfilepT(gfw, corrconfigs.at(1), *c22dpt, i, fPtAxis.get(), cent);
          fillProfilepT(gfw, corrconfigs.at(2), *c24dpt, i, fPtAxis.get(), cent);
          fillProfilepT(gfw, corrconfigs.at(3), *c22Fulldpt, i, fPtAxis.get(), cent);
        }

        for (int i = 1; i <= nK0sPtBins; ++i) {
          fillProfilepTMass(gfw, corrconfigs.at(12), *K0sc22dpt, i, fK0sPtAxis.get(), fK0sMass.get(), kMassBins[0], nK0sPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(13), *K0sc22dpt, i, fK0sPtAxis.get(), fK0sMass.get(), kMassBins[0], nK0sPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(14), *K0sc24dpt, i, fK0sPtAxis.get(), fK0sMass.get(), kMassBins[0], nK0sPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(15), *K0sc22Fulldpt, i, fK0sPtAxis.get(), fK0sMass.get(), kMassBins[0], nK0sPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(27), *K0sc32dpt, i, fK0sPtAxis.get(), fK0sMass.get(), kMassBins[0], nK0sPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(28), *K0sc32dpt, i, fK0sPtAxis.get(), fK0sMass.get(), kMassBins[0], nK0sPtBins, cent);
        }
        for (int i = 1; i <= nLambdaPtBins; ++i) {
          fillProfilepTMass(gfw, corrconfigs.at(16), *Lambdac22dpt, i, fLambdaPtAxis.get(), fLambdaMass.get(), kMassBins[1], nLambdaPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(17), *Lambdac22dpt, i, fLambdaPtAxis.get(), fLambdaMass.get(), kMassBins[1], nLambdaPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(18), *Lambdac24dpt, i, fLambdaPtAxis.get(), fLambdaMass.get(), kMassBins[1], nLambdaPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(19), *Lambdac22Fulldpt, i, fLambdaPtAxis.get(), fLambdaMass.get(), kMassBins[1], nLambdaPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(29), *Lambdac32dpt, i, fLambdaPtAxis.get(), fLambdaMass.get(), kMassBins[1], nLambdaPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(30), *Lambdac32dpt, i, fLambdaPtAxis.get(), fLambdaMass.get(), kMassBins[1], nLambdaPtBins, cent);
        }
        for (int i = 1; i <= nXiPtBins; ++i) {
          fillProfilepTMass(gfw, corrconfigs.at(4), *Xic22dpt, i, fXiPtAxis.get(), fXiMass.get(), kMassBins[2], nXiPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(5), *Xic22dpt, i, fXiPtAxis.get(), fXiMass.get(), kMassBins[2], nXiPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(6), *Xic24dpt, i, fXiPtAxis.get(), fXiMass.get(), kMassBins[2], nXiPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(7), *Xic22Fulldpt, i, fXiPtAxis.get(), fXiMass.get(), kMassBins[2], nXiPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(23), *Xic32dpt, i, fXiPtAxis.get(), fXiMass.get(), kMassBins[2], nXiPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(24), *Xic32dpt, i, fXiPtAxis.get(), fXiMass.get(), kMassBins[2], nXiPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(41), *Xic24_gapdpt, i, fXiPtAxis.get(), fXiMass.get(), kMassBins[2], nXiPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(42), *Xic24_gapdpt, i, fXiPtAxis.get(), fXiMass.get(), kMassBins[2], nXiPtBins, cent);
        }
        for (int i = 1; i <= nOmegaPtBins; ++i) {
          fillProfilepTMass(gfw, corrconfigs.at(8), *Omegac22dpt, i, fOmegaPtAxis.get(), fOmegaMass.get(), kMassBins[3], nOmegaPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(9), *Omegac22dpt, i, fOmegaPtAxis.get(), fOmegaMass.get(), kMassBins[3], nOmegaPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(10), *Omegac24dpt, i, fOmegaPtAxis.get(), fOmegaMass.get(), kMassBins[3], nOmegaPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(11), *Omegac22Fulldpt, i, fOmegaPtAxis.get(), fOmegaMass.get(), kMassBins[3], nOmegaPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(25), *Omegac32dpt, i, fOmegaPtAxis.get(), fOmegaMass.get(), kMassBins[3], nOmegaPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(26), *Omegac32dpt, i, fOmegaPtAxis.get(), fOmegaMass.get(), kMassBins[3], nOmegaPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(43), *Omegac24_gapdpt, i, fOmegaPtAxis.get(), fOmegaMass.get(), kMassBins[3], nOmegaPtBins, cent);
          fillProfilepTMass(gfw, corrconfigs.at(44), *Omegac24_gapdpt, i, fOmegaPtAxis.get(), fOmegaMass.get(), kMassBins[3], nOmegaPtBins, cent);
        }
      });

  TFile outputFileHandle(outputFile, "RECREATE");
  if (outputFileHandle.IsZombie()) {
    throw std::runtime_error(std::string("Cannot create output file: ") + outputFile);
  }

  TDirectory *taskDirectory = outputFileHandle.mkdir("flow-gfw-omega-xi");
  writeObject(*taskDirectory, *hPhi);
  writeObject(*taskDirectory, *hPhicorr);
  writeObject(*taskDirectory, *hPhiK0s);
  writeObject(*taskDirectory, *hPhiLambda);
  writeObject(*taskDirectory, *hPhiXi);
  writeObject(*taskDirectory, *hPhiOmega);
  writeObject(*taskDirectory, *hPhiK0scorr);
  writeObject(*taskDirectory, *hPhiLambdacorr);
  writeObject(*taskDirectory, *hPhiXicorr);
  writeObject(*taskDirectory, *hPhiOmegacorr);
  writeObject(*taskDirectory, *hEta);
  writeObject(*taskDirectory, *hVtxZ);
  writeObject(*taskDirectory, *hMult);
  writeObject(*taskDirectory, *hMultTPC);
  writeObject(*taskDirectory, *hCent);
  writeObject(*taskDirectory, *hPt);
  writeObject(*taskDirectory, hEventCount);
  writeObject(*taskDirectory, *hInteractionRate);
  writeObject(*taskDirectory, *hNTracksPVvsCentrality);
  writeObject(*taskDirectory, *hmultFV0AvsmultFT0A);

  writeObject(*taskDirectory, *c22);
  writeObject(*taskDirectory, *c32);
  writeObject(*taskDirectory, *c24);
  writeObject(*taskDirectory, *c22Full);
  writeObject(*taskDirectory, *c22dpt);
  writeObject(*taskDirectory, *c24dpt);
  writeObject(*taskDirectory, *c22Fulldpt);

  writeObject(*taskDirectory, *Xic22dpt);
  writeObject(*taskDirectory, *Xic24dpt);
  writeObject(*taskDirectory, *Xic22Fulldpt);
  writeObject(*taskDirectory, *Xic24_gapdpt);
  writeObject(*taskDirectory, *Xic32dpt);
  writeObject(*taskDirectory, *Omegac22dpt);
  writeObject(*taskDirectory, *Omegac24dpt);
  writeObject(*taskDirectory, *Omegac22Fulldpt);
  writeObject(*taskDirectory, *Omegac24_gapdpt);
  writeObject(*taskDirectory, *Omegac32dpt);
  writeObject(*taskDirectory, *K0sc22dpt);
  writeObject(*taskDirectory, *K0sc24dpt);
  writeObject(*taskDirectory, *K0sc22Fulldpt);
  writeObject(*taskDirectory, *K0sc32dpt);
  writeObject(*taskDirectory, *Lambdac22dpt);
  writeObject(*taskDirectory, *Lambdac24dpt);
  writeObject(*taskDirectory, *Lambdac22Fulldpt);
  writeObject(*taskDirectory, *Lambdac32dpt);

  writeObject(*taskDirectory, *InvMassXi_all);
  writeObject(*taskDirectory, *InvMassXi);
  writeObject(*taskDirectory, *InvMassOmega_all);
  writeObject(*taskDirectory, *InvMassOmega);
  writeObject(*taskDirectory, *InvMassK0s_all);
  writeObject(*taskDirectory, *InvMassK0s);
  writeObject(*taskDirectory, *InvMassLambda_all);
  writeObject(*taskDirectory, *InvMassLambda);
  writeObject(*taskDirectory, *InvMassALambda_all);
  writeObject(*taskDirectory, *InvMassALambda);

  outputFileHandle.Close();
  std::cout << "Processed " << processedEvents << " AMPT events; wrote "
            << outputFile << "." << std::endl;
}

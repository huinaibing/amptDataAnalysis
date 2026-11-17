#include "../script/calculate_v2ptrho.cpp"

#include "TSystem.h"
#include "TTree.h"

#include <array>
#include <fstream>

namespace ampt_rho_event_selection_test {
void require(bool ok, const char *message) {
  if (!ok) {
    throw std::runtime_error(message);
  }
}

void near(double actual, double expected, const char *message) {
  require(std::abs(actual - expected) < 1.e-8 * (1. + std::abs(expected)),
          message);
}

int labelBin(const TProfile2D &profile, const char *label) {
  for (int bin = 1; bin <= profile.GetNbinsY(); ++bin) {
    if (std::string(profile.GetYaxis()->GetBinLabel(bin)) == label) {
      return bin;
    }
  }
  throw std::runtime_error(std::string("Missing profile ") + label);
}
} // namespace ampt_rho_event_selection_test

void test_v2ptrho_event_selection(
    const char *configFile = "../config/calculate_v2ptrho.json") {
  using namespace ampt_analysis;
  using namespace ampt_rho_event_selection_test;

  const auto config = loadAnalysisConfig(configFile, "v2_pt_rho_output");
  require(!config.v2PtRhoNchCentrality.enabled &&
              config.v2PtRhoUseImpactParameterCentrality,
          "Default config must select legacy impact-parameter centrality");
  require(config.v2PtRhoUseXinliEventSelection,
          "Default config must enable the comparison event gate");

  // A, negative B, eta=0 B, positive B, C. Include valid global gates with
  // only one track in one eta half, and events rejected only by each cut.
  const std::vector<std::array<int, 5>> events{
      {1, 1, 0, 0, 1}, {1, 1, 0, 1, 1}, {2, 2, 0, 0, 1},
      {1, 0, 0, 2, 2}, {1, 1, 0, 1, 2}, {0, 2, 0, 1, 2},
      {2, 2, 0, 1, 0}, {2, 1, 0, 1, 2}, {2, 0, 0, 0, 2},
      {0, 1, 0, 0, 1}, {2, 0, 2, 0, 1}, {1, 0, 2, 0, 2}};
  const std::string dir = "/tmp/v2ptrho_event_selection_test";
  gSystem->mkdir(dir.c_str(), true);
  {
    TFile file((dir + "/input_0.root").c_str(), "RECREATE");
    TTree tree("particles", "particles");
    int eventID = 0, pdgPid = 211, nParticles = 1000;
    double imp = 0., px = 0., py = 0., pz = 0.;
    tree.Branch("eventID", &eventID);
    tree.Branch("nParticles", &nParticles);
    tree.Branch("imp", &imp);
    tree.Branch("pdgPid", &pdgPid);
    tree.Branch("p_x", &px);
    tree.Branch("p_y", &py);
    tree.Branch("p_z", &pz);
    const double eta[]{-0.6, -0.2, 0., 0.2, 0.6};
    for (std::size_t event = 0; event < events.size(); ++event) {
      eventID = static_cast<int>(event);
      imp = std::sqrt((event + 1.5) * config.centralityReferenceRadiusSquared /
                      config.centralityPi);
      int track = 0;
      for (int region = 0; region < 5; ++region) {
        for (int i = 0; i < events[event][region]; ++i, ++track) {
          const double pt = 0.5 + 0.1 * track;
          const double phi = 0.27 + 0.59 * track;
          px = pt * std::cos(phi);
          py = pt * std::sin(phi);
          pz = pt * std::sinh(eta[region]);
          tree.Fill();
        }
      }
    }
    tree.Write();
  }
  const nlohmann::json inputConfig = nlohmann::json::array(
      {{{"path", dir + "/input_"}, {"n_files", 1}, {"bin_val", 0}}});
  {
    std::ofstream file(dir + "/input.json");
    file << inputConfig;
  }
  calculate_v2ptrho((dir + "/input.json").c_str(),
                     (dir + "/result.root").c_str(), 1, 1, configFile);

  TFile file((dir + "/result.root").c_str());
  auto *charged = file.Get<FlowContainer>("pid-flow-pt-corr/FlowContainerCharged");
  auto *unit = file.Get<FlowContainer>("pid-flow-pt-corr/FlowContainerChargedUnitEvent");
  auto *pion = file.Get<FlowContainer>("pid-flow-pt-corr/FlowContainerPi");
  auto *hCent = file.Get<TH1D>("pid-flow-pt-corr/hCent");
  require(charged && unit && pion && hCent, "Existing output paths missing");
  require(!file.Get("pid-flow-pt-corr/nchCentralityCalibration"),
          "Imp mode must not write active Nch calibration metadata");
  near(hCent->GetEntries(), events.size(), "QA must count every input event");
  require(charged->GetProfile()->GetNbinsY() == 19 &&
              unit->GetProfile()->GetNbinsY() == 7,
          "No profiles may be removed or duplicated");

  for (std::size_t event = 0; event < events.size(); ++event) {
    const auto &e = events[event];
    const int a = e[0], b = e[1] + e[2] + e[3], c = e[4];
    const int n = a + b + c, m = a + e[1], p = c + e[3];
    const bool valid = n >= 3 && a > 0 && b > 1 && c > 0 && !(a == 1 && c == 1);
    const bool four = valid && m > 1 && p > 1;
    const double w3 = static_cast<double>(a) * b * c;
    const double d2 = static_cast<double>(n) * (n - 1);
    const double d4 = d2 * (n - 2) * (n - 3);
    const int x = charged->GetProfile()->GetXaxis()->FindBin(event + 1.5);
    const std::array<const char *, 9> labels{
        "covV2Pt", "c22TrackWeight", "hMeanPt", "hMeanPtFull", "ptAveFull",
        "ptSquareAveFull", "c22Full", "c24", "c22Gap0"};
    const double weights[]{valid ? w3 : 0., valid ? w3 : 0.,
                           valid ? static_cast<double>(b) : 0.,
                           valid ? static_cast<double>(n) : 0.,
                           valid ? d2 : 0., valid ? d2 : 0., valid ? d2 : 0.,
                           four ? d4 : 0., valid ? static_cast<double>(m) * p : 0.};
    for (std::size_t i = 0; i < labels.size(); ++i) {
      auto &profile = *charged->GetProfile();
      const int y = labelBin(profile, labels[i]);
      near(profile.GetBinEntries(profile.GetBin(x, y)), weights[i], labels[i]);
      double subsampleWeight = 0.;
      for (auto *object : *charged->GetSubProfiles()) {
        auto &sub = *static_cast<TProfile2D *>(object);
        subsampleWeight += sub.GetBinEntries(sub.GetBin(x, y));
      }
      near(subsampleWeight, weights[i], "Bootstrap must apply the same event gate");
    }
    // Quantities outside the nine comparison inputs must still accept sparse
    // events under their original conditions.
    auto &old = *charged->GetProfile();
    near(old.GetBinEntries(old.GetBin(x, labelBin(old, "ptAve"))),
         b > 1 ? b * (b - 1.) : 0., "Central pair profile changed");
    near(old.GetBinEntries(old.GetBin(x, labelBin(old, "c22"))),
         b > 0 ? a * c : 0., "Original charged gap sample changed");
    auto &u = *unit->GetProfile();
    near(u.GetBinEntries(u.GetBin(x, labelBin(u, "hMeanPt"))),
         b > 0 ? 1. : 0., "Unit-event branch acquired the comparison gate");
    near(u.GetBinEntries(u.GetBin(x, labelBin(u, "c24"))),
         b > 0 && n >= 4 ? 1. : 0., "Unit-event full four acquired the half cut");
    auto &pi = *pion->GetProfile();
    near(pi.GetBinEntries(pi.GetBin(x, labelBin(pi, "hMeanPt"))), b,
         "PID branch acquired the comparison gate");
  }
  std::cout << "PASS: comparison event gates, half-eta four-particle cut, "
               "bootstrap weights, unchanged unit/PID/other charged profiles, "
               "legacy imp centrality and output compatibility.\n";
}

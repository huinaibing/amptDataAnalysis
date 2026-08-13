#include "../corrConfigManager.h"
#include "../dataFrame/fileName.h"
#include "../dataLoader.h"
#include "../eventManager.h"
#include "../selection.h"
#include "../utils.h"
#include "PWGCF/GenericFramework/Core/FlowContainer.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "TCanvas.h"
#include "TPDGCode.h"
#include "TProfile.h"
#include "TProfile3D.h"
#include "TRandom3.h"
#include <iostream>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define GENERATE_AXIS(start, end, step) { \
    []() { \
        std::vector<double> v; \
        for (double x = start; x <= end + 1e-9; x += step) { \
            v.push_back(x); \
        } \
        return v; }()}

/**
 * @brief 记得改路径
 * @details 格式 (path/to/file, centrality, file number)
 *
 */
std::vector<CentralityConfig> cent_configs = CentConfigReader::load("./cent_cfg.json");

#define CF(cfg_type, out, name) CalculateAndFill(gfw, mgr, cfg_type, out, name, evtCent, rndm)

/**
 * @details things need to init
 * 1. magic number
 * 2. flow container
 * 3. gfw
 * 4. CorrConfigManager
 * 5. AMPTEventReader
 */
void calculate_v2_charged()
{
    /// @note configure
    AxisSpec axisMultiplicity{{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90}, "Centrality (%)"};
    AxisSpec axisBootstrap{{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30}, "bootstrap"};
    AxisSpec axisMeanPt{{GENERATE_AXIS(0, 3, 0.003)},
                        "meanPt for c22pt"};
    int cfgFlowNbootstrap = 30;
    // end configure

    /// @note hist for new method
    auto hPion = new TProfile3D("hPion", "hPion", axisMeanPt.getNbins(), axisMeanPt.binEdges.data(),
                                axisMultiplicity.getNbins(), axisMultiplicity.binEdges.data(),
                                axisBootstrap.getNbins(), axisBootstrap.binEdges.data());
    auto hChargedPionFull = new TProfile3D("hChargedPionFull", "hChargedPionFull", axisMeanPt.getNbins(), axisMeanPt.binEdges.data(),
                                           axisMultiplicity.getNbins(), axisMultiplicity.binEdges.data(),
                                           axisBootstrap.getNbins(), axisBootstrap.binEdges.data());
    auto hPionPion = new TProfile3D("hPionPion", "hPionPion", axisMeanPt.getNbins(), axisMeanPt.binEdges.data(),
                                    axisMultiplicity.getNbins(), axisMultiplicity.binEdges.data(),
                                    axisBootstrap.getNbins(), axisBootstrap.binEdges.data());
    auto hPionMeanpt = new TProfile3D("hPionMeanpt", "hPionMeanpt", axisMeanPt.getNbins(), axisMeanPt.binEdges.data(),
                                      axisMultiplicity.getNbins(), axisMultiplicity.binEdges.data(),
                                      axisBootstrap.getNbins(), axisBootstrap.binEdges.data());

    auto hKaon = new TProfile3D("hKaon", "hKaon", axisMeanPt.getNbins(), axisMeanPt.binEdges.data(),
                                axisMultiplicity.getNbins(), axisMultiplicity.binEdges.data(),
                                axisBootstrap.getNbins(), axisBootstrap.binEdges.data());
    auto hChargedKaonFull = new TProfile3D("hChargedKaonFull", "hChargedKaonFull", axisMeanPt.getNbins(), axisMeanPt.binEdges.data(),
                                           axisMultiplicity.getNbins(), axisMultiplicity.binEdges.data(),
                                           axisBootstrap.getNbins(), axisBootstrap.binEdges.data());
    auto hKaonKaon = new TProfile3D("hKaonKaon", "hKaonKaon", axisMeanPt.getNbins(), axisMeanPt.binEdges.data(),
                                    axisMultiplicity.getNbins(), axisMultiplicity.binEdges.data(),
                                    axisBootstrap.getNbins(), axisBootstrap.binEdges.data());
    auto hKaonMeanpt = new TProfile3D("hKaonMeanpt", "hKaonMeanpt", axisMeanPt.getNbins(), axisMeanPt.binEdges.data(),
                                      axisMultiplicity.getNbins(), axisMultiplicity.binEdges.data(),
                                      axisBootstrap.getNbins(), axisBootstrap.binEdges.data());

    auto hProton = new TProfile3D("hProton", "hProton", axisMeanPt.getNbins(), axisMeanPt.binEdges.data(),
                                  axisMultiplicity.getNbins(), axisMultiplicity.binEdges.data(),
                                  axisBootstrap.getNbins(), axisBootstrap.binEdges.data());
    auto hChargedProtonFull = new TProfile3D("hChargedProtonFull", "hChargedProtonFull", axisMeanPt.getNbins(), axisMeanPt.binEdges.data(),
                                             axisMultiplicity.getNbins(), axisMultiplicity.binEdges.data(),
                                             axisBootstrap.getNbins(), axisBootstrap.binEdges.data());
    auto hProtonProton = new TProfile3D("hProtonProton", "hProtonProton", axisMeanPt.getNbins(), axisMeanPt.binEdges.data(),
                                        axisMultiplicity.getNbins(), axisMultiplicity.binEdges.data(),
                                        axisBootstrap.getNbins(), axisBootstrap.binEdges.data());
    auto hProtonMeanpt = new TProfile3D("hProtonMeanpt", "hProtonMeanpt", axisMeanPt.getNbins(), axisMeanPt.binEdges.data(),
                                        axisMultiplicity.getNbins(), axisMultiplicity.binEdges.data(),
                                        axisBootstrap.getNbins(), axisBootstrap.binEdges.data());
    // end hist for new method

#pragma region // flow container init
    TRandom3 *fRndm = new TRandom3(0);
    FlowContainer *fFCCh = new FlowContainer("FlowContainerCharged");
    FlowContainer *fFCPi = new FlowContainer("FlowContainerPi");
    FlowContainer *fFCKa = new FlowContainer("FlowContainerKa");
    FlowContainer *fFCPr = new FlowContainer("FlowContainerPr");

    TObjArray *oba4Ch = new TObjArray();
    oba4Ch->Add(new TNamed("c22", "c22"));
    oba4Ch->Add(new TNamed("c32", "c32"));
    oba4Ch->Add(new TNamed("c24", "c24"));
    oba4Ch->Add(new TNamed("c34", "c34"));
    oba4Ch->Add(new TNamed("c22Full", "c22Full"));
    oba4Ch->Add(new TNamed("c22TrackWeight", "c22TrackWeight"));
    oba4Ch->Add(new TNamed("c32TrackWeight", "c32TrackWeight"));
    oba4Ch->Add(new TNamed("c24TrackWeight", "c24TrackWeight"));
    oba4Ch->Add(new TNamed("c34TrackWeight", "c34TrackWeight"));
    oba4Ch->Add(new TNamed("c22FullTrackWeight", "c22FullTrackWeight"));
    oba4Ch->Add(new TNamed("covV2Pt", "covV2Pt"));
    oba4Ch->Add(new TNamed("covV3Pt", "covV3Pt"));
    oba4Ch->Add(new TNamed("ptSquareAve", "ptSquareAve"));
    oba4Ch->Add(new TNamed("ptAve", "ptAve"));
    oba4Ch->Add(new TNamed("hMeanPt", "hMeanPt"));
    // end fill TObjArray for charged

    // init fFCCh
    fFCCh->SetName("FlowContainerCharged");
    fFCCh->Initialize(oba4Ch, axisMultiplicity, cfgFlowNbootstrap);
    // end init fFCCh

    // init fFCPID
    // note that need to add c22pure and c32pure
    TObjArray *oba4PID = reinterpret_cast<TObjArray *>(oba4Ch->Clone());
    oba4PID->Add(new TNamed("c22pure", "c22pure"));
    oba4PID->Add(new TNamed("c32pure", "c32pure"));
    oba4PID->Add(new TNamed("covV2PtPID", "covV2PtPID"));
    oba4PID->Add(new TNamed("c22TrackWeightPID", "c22TrackWeightPID"));

    fFCPi->SetName("FlowContainerPi");
    fFCPi->Initialize(oba4PID, axisMultiplicity, cfgFlowNbootstrap);

    fFCKa->SetName("FlowContainerKa");
    fFCKa->Initialize(oba4PID, axisMultiplicity, cfgFlowNbootstrap);

    fFCPr->SetName("FlowContainerPr");
    fFCPr->Initialize(oba4PID, axisMultiplicity, cfgFlowNbootstrap);
#pragma endregion // end flow container init

    /// @note gfw init
    GFW *gfw = new GFW();
    CorrConfigManager mgr(gfw);
    // end gfw init

    /// @note diff cent loop
    for (const auto &cfg : cent_configs)
    {
        AMPTEventReader reader(cfg.path.c_str(), cfg.n_files); // 使用配置里的 n_files
        /// @note event loop
        for (const auto &evt : reader)
        {
            gfw->Clear();
            float rndm = fRndm->Rndm();
            double evtCent = TMath::Pi() * evt.imp * evt.imp / 7.84;
            /// @note particle loop
            for (const auto &trk : evt.particles)
            {
                if (particleSelected(trk))
                {
                    gfw->Fill(trk.GetEta(), 0, trk.GetPhi(), 1, Mask::kRef);

                    if (TMath::Abs(trk.pdgPid) == 211)
                        gfw->Fill(trk.GetEta(), 0, trk.GetPhi(), 1, Mask::kPion | Mask::kPionOverlap);
                    if (TMath::Abs(trk.pdgPid) == 321)
                        gfw->Fill(trk.GetEta(), 0, trk.GetPhi(), 1, Mask::kKaon | Mask::kKaonOverlap);
                    if (TMath::Abs(trk.pdgPid) == 2212)
                        gfw->Fill(trk.GetEta(), 0, trk.GetPhi(), 1, Mask::kProton | Mask::kProtonOverlap);
                }
            }
            // end particle loop

            /// @note charged flow calculation
            {
                CF(CorrType::Ref08Gap22, fFCCh, "c22");
                CF(CorrType::Ref0Gap24, fFCCh, "c24");
                CF(CorrType::Ref0Gap22, fFCCh, "c22Full");

                CalculateCovV2ChargedPt(gfw, mgr, CorrType::Ref08Gap22, fFCCh, evt, evtCent, cut4Pt, rndm);
                CalculateC22TrackWeight(gfw, mgr, CorrType::Ref08Gap22, fFCCh, evt, evtCent, cut4Pt, rndm);

                double nParticlesCh = evt.nParticlesAfterCut(cut4Pt);
                if (nParticlesCh > 1)
                {
                    fFCCh->FillProfile("hMeanPt", evtCent, evt.GetMeanPt(cut4Pt), nParticlesCh, rndm);
                    fFCCh->FillProfile("ptAve", evtCent, evt.GetMeanPt(cut4Pt), nParticlesCh * nParticlesCh - nParticlesCh, rndm);
                    fFCCh->FillProfile("ptSquareAve", evtCent, evt.GetPtSquareAve(cut4Pt), nParticlesCh * nParticlesCh - nParticlesCh, rndm);
                }
            } // end charged flow calculation

            /// @note pion
            {
                CF(CorrType::Pion08Gap22a, fFCPi, "c22");
                CF(CorrType::Pion08Gap22b, fFCPi, "c22");

                CF(CorrType::PiPi08Gap22, fFCPi, "c22pure");

                CF(CorrType::Pion0Gap22a_Full, fFCPi, "c22Full");
                CF(CorrType::Pion0Gap22b_Full, fFCPi, "c22Full");

                CF(CorrType::Pion0Gap24a, fFCPi, "c24");
                CF(CorrType::Pion0Gap24b, fFCPi, "c24");

                CalculateCovV2ChargedPt(gfw, mgr, CorrType::PiPi08Gap22, fFCPi, evt, evtCent, cut4Pt, rndm);
                CalculateC22TrackWeight(gfw, mgr, CorrType::PiPi08Gap22, fFCPi, evt, evtCent, cut4Pt, rndm);

                std::function<bool(const Track &)> cut4PtPi = [](const Track &trk)
                { if(cut4Pt(trk)) if(TMath::Abs(trk.pdgPid) == PDG_t::kPiPlus) return true; return false; };

                CalculateCovV2ChargedPt(gfw, mgr, CorrType::PiPi08Gap22, fFCPi, evt, evtCent, cut4PtPi, rndm, "covV2PtPID");
                CalculateC22TrackWeight(gfw, mgr, CorrType::PiPi08Gap22, fFCPi, evt, evtCent, cut4PtPi, rndm, "c22TrackWeightPID");

                double nParticlesPi = evt.nParticlesAfterCut(cut4PtPi);
                if (nParticlesPi > 1)
                {
                    FillMeanptCentBSProfile(evtCent, rndm,
                                            evt.GetMeanPt(cut4PtPi), nParticlesPi,
                                            gfw, mgr, CorrType::Pion08Gap22a, CorrType::Pion08Gap22b, CorrType::PiPi08Gap22,
                                            hPion, hChargedPionFull, hPionPion, hPionMeanpt);
                    fFCPi->FillProfile("hMeanPt", evtCent, evt.GetMeanPt(cut4PtPi), nParticlesPi, rndm);
                    fFCPi->FillProfile("ptAve", evtCent, evt.GetMeanPt(cut4PtPi), nParticlesPi * nParticlesPi - nParticlesPi, rndm);
                    fFCPi->FillProfile("ptSquareAve", evtCent, evt.GetPtSquareAve(cut4PtPi), nParticlesPi * nParticlesPi - nParticlesPi, rndm);
                }
            }
            // end pion

            /// @note kaon
            {
                CF(CorrType::Kaon08Gap22a, fFCKa, "c22");
                CF(CorrType::Kaon08Gap22b, fFCKa, "c22");

                CF(CorrType::KaKa08Gap22, fFCKa, "c22pure");

                CF(CorrType::Kaon0Gap22a_Full, fFCKa, "c22Full");
                CF(CorrType::Kaon0Gap22b_Full, fFCKa, "c22Full");

                CF(CorrType::Kaon0Gap24a, fFCKa, "c24");
                CF(CorrType::Kaon0Gap24b, fFCKa, "c24");

                CalculateCovV2ChargedPt(gfw, mgr, CorrType::KaKa08Gap22, fFCKa, evt, evtCent, cut4Pt, rndm);
                CalculateC22TrackWeight(gfw, mgr, CorrType::KaKa08Gap22, fFCKa, evt, evtCent, cut4Pt, rndm);

                std::function<bool(const Track &)> cut4PtKa = [](const Track &trk)
                { if(cut4Pt(trk)) if(TMath::Abs(trk.pdgPid) == PDG_t::kKPlus) return true; return false; };

                CalculateCovV2ChargedPt(gfw, mgr, CorrType::KaKa08Gap22, fFCKa, evt, evtCent, cut4PtKa, rndm, "covV2PtPID");
                CalculateC22TrackWeight(gfw, mgr, CorrType::KaKa08Gap22, fFCKa, evt, evtCent, cut4PtKa, rndm, "c22TrackWeightPID");

                double nParticlesKa = evt.nParticlesAfterCut(cut4PtKa);
                if (nParticlesKa > 1)
                {
                    FillMeanptCentBSProfile(evtCent, rndm,
                                            evt.GetMeanPt(cut4PtKa), nParticlesKa,
                                            gfw, mgr, CorrType::Kaon08Gap22a, CorrType::Kaon08Gap22b, CorrType::KaKa08Gap22,
                                            hKaon, hChargedKaonFull, hKaonKaon, hKaonMeanpt);
                    fFCKa->FillProfile("hMeanPt", evtCent, evt.GetMeanPt(cut4PtKa), nParticlesKa, rndm);
                    fFCKa->FillProfile("ptAve", evtCent, evt.GetMeanPt(cut4PtKa), nParticlesKa * nParticlesKa - nParticlesKa, rndm);
                    fFCKa->FillProfile("ptSquareAve", evtCent, evt.GetPtSquareAve(cut4PtKa), nParticlesKa * nParticlesKa - nParticlesKa, rndm);
                }
            }
            // end kaon

            /// @note proton
            {
                CF(CorrType::Prot08Gap22a, fFCPr, "c22");
                CF(CorrType::Prot08Gap22b, fFCPr, "c22");

                CF(CorrType::PrPr08Gap22, fFCPr, "c22pure");

                CF(CorrType::Prot0Gap22a_Full, fFCPr, "c22Full");
                CF(CorrType::Prot0Gap22b_Full, fFCPr, "c22Full");

                CF(CorrType::Prot0Gap24a, fFCPr, "c24");
                CF(CorrType::Prot0Gap24b, fFCPr, "c24");

                CalculateCovV2ChargedPt(gfw, mgr, CorrType::PrPr08Gap22, fFCPr, evt, evtCent, cut4Pt, rndm);
                CalculateC22TrackWeight(gfw, mgr, CorrType::PrPr08Gap22, fFCPr, evt, evtCent, cut4Pt, rndm);

                std::function<bool(const Track &)> cut4PtPr = [](const Track &trk)
                { if(cut4Pt(trk)) if(TMath::Abs(trk.pdgPid) == PDG_t::kProton) return true; return false; };

                CalculateCovV2ChargedPt(gfw, mgr, CorrType::PrPr08Gap22, fFCPr, evt, evtCent, cut4PtPr, rndm, "covV2PtPID");
                CalculateC22TrackWeight(gfw, mgr, CorrType::PrPr08Gap22, fFCPr, evt, evtCent, cut4PtPr, rndm, "c22TrackWeightPID");

                double nParticlesPr = evt.nParticlesAfterCut(cut4PtPr);
                if (nParticlesPr > 1)
                {
                    FillMeanptCentBSProfile(evtCent, rndm,
                                            evt.GetMeanPt(cut4PtPr), nParticlesPr,
                                            gfw, mgr, CorrType::Prot08Gap22a, CorrType::Prot08Gap22b, CorrType::PrPr08Gap22,
                                            hProton, hChargedProtonFull, hProtonProton, hProtonMeanpt);
                    fFCPr->FillProfile("hMeanPt", evtCent, evt.GetMeanPt(cut4PtPr), nParticlesPr, rndm);
                    fFCPr->FillProfile("ptAve", evtCent, evt.GetMeanPt(cut4PtPr), nParticlesPr * nParticlesPr - nParticlesPr, rndm);
                    fFCPr->FillProfile("ptSquareAve", evtCent, evt.GetPtSquareAve(cut4PtPr), nParticlesPr * nParticlesPr - nParticlesPr, rndm);
                }
            }
            // end proton
        }
        // end event loop
    }
    // end diff cent loop

#pragma region // save to file
    TFile *outputAnalysisResult = new TFile("myAnalysisResult.root", "RECREATE");
    outputAnalysisResult->mkdir("pid-flow-pt-corr");
    outputAnalysisResult->cd("pid-flow-pt-corr");
    fFCCh->Write();
    fFCPi->Write();
    fFCKa->Write();
    fFCPr->Write();

    outputAnalysisResult->mkdir("pid-flow-pt-corr/meanptCentNbs");
    outputAnalysisResult->cd("pid-flow-pt-corr/meanptCentNbs");

    hPion->Write();
    hChargedPionFull->Write();
    hPionPion->Write();

    hKaon->Write();
    hChargedKaonFull->Write();
    hKaonKaon->Write();

    hProton->Write();
    hChargedProtonFull->Write();
    hProtonProton->Write();

    hPionMeanpt->Write();
    hKaonMeanpt->Write();
    hProtonMeanpt->Write();

    outputAnalysisResult->Close();
#pragma endregion

    std::cout << "==================================finish===================================" << std::endl;
}

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
#include "TRandom3.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

/**
 * @brief 记得改路径
 * @details 格式 (path/to/file, centrality, file number)
 *
 */
std::vector<CentralityConfig> cent_configs = {
    {"/home/huinaibing/ampt_result/cent50-60/Result/ampt_19370820_", 55, 100},
    {"/home/huinaibing/ampt_result/cent40-50/Result/ampt_19369195_", 45, 100},
    {"/home/huinaibing/ampt_result/cent30-40/Result/ampt_19367451_", 35, 100},
    {"/home/huinaibing/ampt_result/cent20_30/Result/ampt_19343841_", 25, 100},
    {"/home/huinaibing/ampt_result/Result0-10_4/ampt_16824297_", 5, 100},
    {"/home/huinaibing/ampt_result/Result0-10_3/ampt_16752378_", 5, 100},
    {"/home/huinaibing/ampt_result/Result0-10_2/ampt_16741210_", 5, 100},
    {"/home/huinaibing/ampt_result/Result0-10_1/ampt_16724824_", 5, 100}};

#define CF(cfg_type, out, name) CalculateAndFill(gfw, mgr, cfg_type, out, name, cfg.bin_val, rndm)

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
    AxisSpec axisMultiplicity{{0, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "Centrality (%)"};
    int cfgFlowNbootstrap = 30;
    // end configure

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
        AMPTEventReader reader(cfg.path, cfg.n_files); // 使用配置里的 n_files
        /// @note event loop
        for (const auto &evt : reader)
        {
            gfw->Clear();
            float rndm = fRndm->Rndm();
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

                CalculateCovV2ChargedPt(gfw, mgr, CorrType::Ref08Gap22, fFCCh, evt, cfg.bin_val, rndm);
                CalculateC22TrackWeight(gfw, mgr, CorrType::Ref08Gap22, fFCCh, evt, cfg.bin_val, rndm);

                if (evt.nParticles)
                {
                    fFCCh->FillProfile("hMeanPt", cfg.bin_val, evt.GetMeanPt(), evt.nParticles, rndm);
                    fFCCh->FillProfile("ptAve", cfg.bin_val, evt.GetMeanPt(), evt.nParticles, rndm);
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

                CalculateCovV2ChargedPt(gfw, mgr, CorrType::PiPi08Gap22, fFCPi, evt, cfg.bin_val, rndm);
                CalculateC22TrackWeight(gfw, mgr, CorrType::PiPi08Gap22, fFCPi, evt, cfg.bin_val, rndm);

                if (evt.nPidParticles(PDG_t::kPiPlus))
                {
                    fFCPi->FillProfile("hMeanPt", cfg.bin_val, evt.GetMeanPt(PDG_t::kPiPlus), evt.nPidParticles(PDG_t::kPiPlus), rndm);
                    fFCPi->FillProfile("ptAve", cfg.bin_val, evt.GetMeanPt(PDG_t::kPiPlus), evt.nPidParticles(PDG_t::kPiPlus), rndm);
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

                CalculateCovV2ChargedPt(gfw, mgr, CorrType::KaKa08Gap22, fFCKa, evt, cfg.bin_val, rndm);
                CalculateC22TrackWeight(gfw, mgr, CorrType::KaKa08Gap22, fFCKa, evt, cfg.bin_val, rndm);

                if (evt.nPidParticles(PDG_t::kKPlus))
                {
                    fFCKa->FillProfile("hMeanPt", cfg.bin_val, evt.GetMeanPt(PDG_t::kKPlus), evt.nPidParticles(PDG_t::kKPlus), rndm);
                    fFCKa->FillProfile("ptAve", cfg.bin_val, evt.GetMeanPt(PDG_t::kKPlus), evt.nPidParticles(PDG_t::kKPlus), rndm);
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

                CalculateCovV2ChargedPt(gfw, mgr, CorrType::PrPr08Gap22, fFCPr, evt, cfg.bin_val, rndm);
                CalculateC22TrackWeight(gfw, mgr, CorrType::PrPr08Gap22, fFCPr, evt, cfg.bin_val, rndm);

                if (evt.nPidParticles(PDG_t::kProton))
                {
                    fFCPr->FillProfile("hMeanPt", cfg.bin_val, evt.GetMeanPt(PDG_t::kProton), evt.nPidParticles(PDG_t::kProton), rndm);
                    fFCPr->FillProfile("ptAve", cfg.bin_val, evt.GetMeanPt(PDG_t::kProton), evt.nPidParticles(PDG_t::kProton), rndm);
                }
            }
            // end proton
        }
        // end event loop
    }
    // end diff cent loop

#pragma region // save to file
    TFile *outputAnalysisResult = new TFile("myAnalysisResult.root", "RECREATE");
    fFCCh->Write();
    fFCPi->Write();
    fFCKa->Write();
    fFCPr->Write();

    outputAnalysisResult->Close();
#pragma endregion
}

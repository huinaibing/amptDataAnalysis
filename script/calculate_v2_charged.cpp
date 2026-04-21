#include "../corrConfigManager.h"
#include "../dataLoader.h"
#include "../eventManager.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "TCanvas.h"
#include "TProfile.h"

void calculate_v2_charged()
{

    GFW *gfw = new GFW();
    CorrConfigManager mgr(gfw);
    auto &cfg_ref = mgr.Get(CorrType::Ref08Gap22);

    TProfile *res = new TProfile("res", "", 6, 0, 60);

    // 50-60
    {
        AMPTEventReader cent_50_60("/home/huinaibing/ampt_result/cent50-60/Result/ampt_19370820_", 100);
        for (const auto &evt : cent_50_60)
        {
            gfw->Clear();
            for (const auto &trk : evt.particles)
            {
                gfw->Fill(trk.GetEta(), 0, trk.GetPhi(), 1, 1);
            }

            double cum = gfw->Calculate(cfg_ref, 0, false).real();
            double npair = gfw->Calculate(cfg_ref, 0, true).real();
            if (npair != 0)
                res->Fill(55, cum / npair, npair);
        }
    }

    // cent 40-50
    {
        AMPTEventReader cent_40_50("/home/huinaibing/ampt_result/cent40-50/Result/ampt_19369195_", 100);
        for (const auto &evt : cent_40_50)
        {
            gfw->Clear();
            for (const auto &trk : evt.particles)
            {
                gfw->Fill(trk.GetEta(), 0, trk.GetPhi(), 1, 1);
            }

            double cum = gfw->Calculate(cfg_ref, 0, false).real();
            double npair = gfw->Calculate(cfg_ref, 0, true).real();
            if (npair != 0)
                res->Fill(45, cum / npair, npair);
        }
    }

    // 30-40
    {
        AMPTEventReader cent_30_40("/home/huinaibing/ampt_result/cent30-40/Result/ampt_19367451_", 100);
        for (const auto &evt : cent_30_40)
        {
            gfw->Clear();
            for (const auto &trk : evt.particles)
            {
                gfw->Fill(trk.GetEta(), 0, trk.GetPhi(), 1, 1);
            }

            double cum = gfw->Calculate(cfg_ref, 0, false).real();
            double npair = gfw->Calculate(cfg_ref, 0, true).real();
            if (npair != 0)
                res->Fill(35, cum / npair, npair);
        }
    }

    // 20-30
    {
        AMPTEventReader cent_20_30("/home/huinaibing/ampt_result/cent20_30/Result/ampt_19343841_", 100);
        for (const auto &evt : cent_20_30)
        {
            gfw->Clear();
            for (const auto &trk : evt.particles)
            {
                gfw->Fill(trk.GetEta(), 0, trk.GetPhi(), 1, 1);
            }

            double cum = gfw->Calculate(cfg_ref, 0, false).real();
            double npair = gfw->Calculate(cfg_ref, 0, true).real();
            if (npair != 0)
                res->Fill(25, cum / npair, npair);
        }
    }

    TCanvas *c1 = new TCanvas("c1", "", 800, 600);
    res->Draw();
}

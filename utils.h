#ifndef UTILS
#define UTILS
#include "PWGCF/GenericFramework/Core/FlowContainer.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "corrConfigManager.h"
#include "dataFrame/fileName.h"
#include "dataLoader.h"
#include "eventManager.h"
#include "selection.h"

/**
 * @brief 用来填充c22 c24 这类图
 *
 * @param gfw
 * @param mgr
 * @param cfg
 * @param fOut
 * @param profileName
 * @param bin_val
 * @param rndm
 */
void CalculateAndFill(
    GFW *gfw,
    const CorrConfigManager &mgr,
    CorrType cfg,
    FlowContainer *fOut,
    const char *profileName,
    double bin_val,
    float rndm)
{
    double cum = gfw->Calculate(mgr.Get(cfg), 0, false).real();
    double npair = gfw->Calculate(mgr.Get(cfg), 0, true).real();
    if (npair != 0)
    {
        fOut->FillProfile(profileName, bin_val, cum / npair, npair, rndm);
    }
}

void CalculateCovV2ChargedPt(
    GFW *gfw,
    const CorrConfigManager &mgr,
    CorrType cfg,
    FlowContainer *fOut,
    const Event &evt,
    double bin_val,
    std::function<bool(const Track &)> customCut,
    float rndm)
{
    double cum = gfw->Calculate(mgr.Get(cfg), 0, false).real();
    double npair = gfw->Calculate(mgr.Get(cfg), 0, true).real();
    if (npair == 0)
        return;

    double val = cum / npair;

    fOut->FillProfile("covV2Pt", bin_val, val * evt.GetMeanPt(customCut), npair * evt.nParticlesAfterCut(customCut), rndm);
}

void CalculateC22TrackWeight(
    GFW *gfw,
    const CorrConfigManager &mgr,
    CorrType cfg,
    FlowContainer *fOut,
    const Event &evt,
    double bin_val,
    std::function<bool(const Track &)> customCut,
    float rndm)
{
    double cum = gfw->Calculate(mgr.Get(cfg), 0, false).real();
    double npair = gfw->Calculate(mgr.Get(cfg), 0, true).real();
    if (npair != 0)
    {
        fOut->FillProfile("c22TrackWeight", bin_val, cum / npair, npair * evt.nParticlesAfterCut(customCut), rndm);
    }
}

#endif

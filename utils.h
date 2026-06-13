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

const int cfgFlowNbootstrap = 30;

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
    float rndm,
    const char *profileName = "covV2Pt")
{
    double cum = gfw->Calculate(mgr.Get(cfg), 0, false).real();
    double npair = gfw->Calculate(mgr.Get(cfg), 0, true).real();
    if (npair == 0)
        return;

    double val = cum / npair;

    fOut->FillProfile(profileName, bin_val, val * evt.GetMeanPt(customCut), npair * evt.nParticlesAfterCut(customCut), rndm);
}

void CalculateC22TrackWeight(
    GFW *gfw,
    const CorrConfigManager &mgr,
    CorrType cfg,
    FlowContainer *fOut,
    const Event &evt,
    double bin_val,
    std::function<bool(const Track &)> customCut,
    float rndm,
    const char *profileName = "c22TrackWeight")
{
    double cum = gfw->Calculate(mgr.Get(cfg), 0, false).real();
    double npair = gfw->Calculate(mgr.Get(cfg), 0, true).real();
    if (npair != 0)
    {
        fOut->FillProfile(profileName, bin_val, cum / npair, npair * evt.nParticlesAfterCut(customCut), rndm);
    }
}

double getPidC22InOneEvent(GFW *fGFW, const GFW::CorrConfig &corrconfA, const GFW::CorrConfig &corrconfB)
{
    double NpairA = fGFW->Calculate(corrconfA, 0, true).real();
    double NpairB = fGFW->Calculate(corrconfB, 0, true).real();

    if (NpairA == 0 && NpairB == 0)
        return 0;

    double ChC22A = NpairA ? fGFW->Calculate(corrconfA, 0, false).real() / NpairA : 0.;
    double ChC22B = NpairB ? fGFW->Calculate(corrconfB, 0, false).real() / NpairB : 0.;

    double ChC22 = (ChC22A * NpairA + ChC22B * NpairB) / (NpairA + NpairB);

    return ChC22;
}

/**
 * @brief 自己看变量名字
 *
 * @param cent
 * @param rndm
 * @param pidMeanpt
 * @param nPid
 * @param fGFW
 * @param mgr
 * @param A correlation name
 * @param B correlation name
 * @param POIREF TProfile3D
 * @param REF TProfile3D
 */
void FillMeanptCentBSProfile(const double &cent,
                             const double &rndm,
                             const double &pidMeanpt,
                             const double &nPid,

                             GFW *fGFW,
                             const CorrConfigManager &mgr,
                             CorrType A,
                             CorrType B,
                             CorrType pure,
                             TProfile3D *POIREF,
                             TProfile3D *REF,
                             TProfile3D *POIPOI)
{
    /// @note calculate <2> for charged
    double dnx, val;
    dnx = fGFW->Calculate(mgr.Get(CorrType::Ref08Gap22), 0, true).real();
    if (dnx == 0)
        return;
    val = fGFW->Calculate(mgr.Get(CorrType::Ref08Gap22), 0, false).real() / dnx;
    if (std::fabs(val) >= 1)
        return;
    // end calculate <2> for charged

    /// @note get PID particle's c22(POI - REF)
    double pidc22 = 0;
    double npairPid = 0;

    pidc22 = getPidC22InOneEvent(fGFW, mgr.Get(A), mgr.Get(B));
    if (pidc22 == 0)
        return;

    npairPid = fGFW->Calculate(mgr.Get(A), 0, true).real() + fGFW->Calculate(mgr.Get(B), 0, true).real();
    if (npairPid == 0)
        return;

    POIREF->Fill(pidMeanpt, cent, rndm * cfgFlowNbootstrap, pidc22, nPid * npairPid);
    REF->Fill(pidMeanpt, cent, rndm * cfgFlowNbootstrap, val, dnx * nPid);
    // end get PID particle's c22(POI - REF)

    /// @note get POIPOI
    double npairPure, valPure;
    npairPure = fGFW->Calculate(mgr.Get(pure), 0, true).real();
    if (npairPure == 0)
        return;
    valPure = fGFW->Calculate(mgr.Get(pure), 0, false).real() / npairPure;
    if (std::fabs(valPure) >= 1)
        return;
    // end get POIPOI

    POIPOI->Fill(pidMeanpt, cent, rndm * cfgFlowNbootstrap, valPure, npairPure);
}

#endif

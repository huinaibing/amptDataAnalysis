#ifndef SELECTION
#define SELECTION
#include "TMath.h"
#include "dataFrame/track.h"

/**
 * @brief track selection
 * @details 1. eta 2. pt
 * @param trk
 * @return
 */
bool particleSelected(const Track &trk)
{
    if (TMath::Abs(trk.GetEta()) > 0.8)
        return false;
    if (trk.GetPt() > 3)
        return false;
    if (trk.GetPt() < 0.2)
        return false;

    return true;
}

#endif

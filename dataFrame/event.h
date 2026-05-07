#ifndef EVENT
#define EVENT
#include "track.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
struct Event
{
    int eventID;
    double imp;

    std::vector<Track> particles;

    double GetMeanPt(std::function<bool(const Track &)> customCut) const
    {
        double meanPt = 0;
        for (const auto &track : this->particles)
        {
            if (customCut(track))
                meanPt += track.GetPt();
        }

        double nParticles = this->nParticlesAfterCut(customCut);
        if (nParticles == 0)
            return 0;
        return meanPt / nParticles;
    }

    int nParticlesAfterCut(std::function<bool(const Track &)> customCut) const
    {
        int nPid = 0;
        for (const auto &track : this->particles)
        {
            if (customCut(track))
            {
                nPid++;
            }
        }
        return nPid;
    }

    double GetPtSquareAve(std::function<bool(const Track &)> customCut) const
    {
        double nParticles = this->nParticlesAfterCut(customCut);
        if (nParticles == 0 || nParticles == 1)
            return 0;
        double ptSum = this->GetMeanPt(customCut) * nParticles;

        double ptSquareSum = 0;
        for (const auto &track : this->particles)
        {
            if (customCut(track))
                ptSquareSum += track.GetPt() * track.GetPt();
        }

        return (ptSum * ptSum - ptSquareSum) / (nParticles * nParticles - nParticles);
    }
};

#endif

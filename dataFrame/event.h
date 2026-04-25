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

        return meanPt / this->nParticlesAfterCut(customCut);
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
};

#endif

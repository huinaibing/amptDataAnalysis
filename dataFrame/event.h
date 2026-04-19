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
    int nParticles;
    double imp;

    std::vector<Track> particles;

    double GetMeanPt() const
    {
        double meanPt = 0;
        for (const auto &track : this->particles)
        {
            meanPt += track.GetPt();
        }

        return meanPt / this->nParticles;
    }
};

#endif

#ifndef TRACK
#define TRACK

struct Track
{
    int eventID;
    int nParticles;
    double imp;

    int pdgPid;
    double p_x, p_y, p_z;

    double GetPt() const
    {
        return std::sqrt(p_x * p_x + p_y * p_y);
    }

    double GetEta() const
    {
        double p = std::sqrt(p_x * p_x + p_y * p_y + p_z * p_z);
        if (p == std::abs(p_z))
            return INFINITY;
        return 0.5 * std::log((p + p_z) / (p - p_z));
    }

    double GetPhi() const
    {
        return std::atan2(p_y, p_x);
    }
};

#endif

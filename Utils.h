#include "PWGCF/GenericFramework/Core/GFW.h"
#include <string>
#include <utility>

/// @todo 重构这个文件


#define GFW_REGION_LIST \
    X(reffull)          \
    X(refN08)           \
    X(refP08)           \
    X(refN)             \
    X(refP)             \
    X(poiPiN08)         \
    X(poiPiP08)         \
    X(poiPiN)           \
    X(poiPiP)           \
    X(olPiN)            \
    X(olPiP)            \
    X(poiKaN08)         \
    X(poiKaP08)         \
    X(poiKaN)           \
    X(poiKaP)           \
    X(olKaN)            \
    X(olKaP)            \
    X(poiPrN08)         \
    X(poiPrP08)         \
    X(poiPrN)           \
    X(poiPrP)           \
    X(olPrN)            \
    X(olPrP)

#define GFW_CONFIG_LIST                                                     \
    X(Ref08Gap22, "refP08 {2} refN08 {-2}", "Ref08Gap22")                   \
    X(Ref0Gap24, "refN {2 2} refP {-2 -2}", "Ref0Gap24")                    \
    X(Ref0Gap22, "refN {2} refP {-2}", "Ref0Gap22")                         \
    X(Ref08Gap32, "refP08 {3} refN08 {-3}", "Ref08Gap32")                   \
    X(Ref08Gap34, "refP08 {3 3} refN08 {-3 -3}", "Ref08Gap34")              \
    X(Pion08gap22a, "poiPiN08 {2} refP08 {-2}", "Pion08gap22a")             \
    X(Pion08gap22b, "poiPiP08 {2} refN08 {-2}", "Pion08gap22b")             \
    X(Kaon08gap22a, "poiKaN08 {2} refP08 {-2}", "Kaon08gap22a")             \
    X(Kaon08gap22b, "poiKaP08 {2} refN08 {-2}", "Kaon08gap22b")             \
    X(Prot08gap22a, "poiPrN08 {2} refP08 {-2}", "Prot08gap22a")             \
    X(Prot08gap22b, "poiPrP08 {2} refN08 {-2}", "Prot08gap22b")             \
    X(Pion0gap24a, "poiPiN refN | olPiN {2 2} refP {-2 -2}", "Pion0gap24a") \
    X(Pion0gap24b, "poiPiP refP | olPiP {2 2} refN {-2 -2}", "Pion0gap24b") \
    X(Kaon0gap24a, "poiKaN refN | olKaN {2 2} refP {-2 -2}", "Kaon0gap24a") \
    X(Kaon0gap24b, "poiKaP refP | olKaP {2 2} refN {-2 -2}", "Kaon0gap24b") \
    X(Prot0gap24a, "poiPrN refN | olPrN {2 2} refP {-2 -2}", "Prot0gap24a") \
    X(Prot0gap24b, "poiPrP refP | olPrP {2 2} refN {-2 -2}", "Prot0gap24b") \
    X(Pion08gap32a, "poiPiN08 {3} refP08 {-3}", "Pion08gap32a")             \
    X(Pion08gap32b, "poiPiP08 {3} refN08 {-3}", "Pion08gap32b")             \
    X(Kaon08gap32a, "poiKaN08 {3} refP08 {-3}", "Kaon08gap32a")             \
    X(Kaon08gap32b, "poiKaP08 {3} refN08 {-3}", "Kaon08gap32b")             \
    X(Prot08gap32a, "poiPrN08 {3} refP08 {-3}", "Prot08gap32a")             \
    X(Prot08gap32b, "poiPrP08 {3} refN08 {-3}", "Prot08gap32b")             \
    X(Pion0gap34a, "poiPiN refN | olPiN {3 3} refP {-3 -3}", "Pion0gap34a") \
    X(Pion0gap34b, "poiPiP refP | olPiP {3 3} refN {-3 -3}", "Pion0gap34b") \
    X(Kaon0gap34a, "poiKaN refN | olKaN {3 3} refP {-3 -3}", "Kaon0gap34a") \
    X(Kaon0gap34b, "poiKaP refP | olKaP {3 3} refN {-3 -3}", "Kaon0gap34b") \
    X(Prot0gap34a, "poiPrN refN | olPrN {3 3} refP {-3 -3}", "Prot0gap34a") \
    X(Prot0gap34b, "poiPrP refP | olPrP {3 3} refN {-3 -3}", "Prot0gap34b") \
    X(PiPi08gap22, "poiPiN08 {2} poiPiP08 {-2}", "PiPi08gap22")             \
    X(KaKa08gap22, "poiKaN08 {2} poiKaP08 {-2}", "KaKa08gap22")             \
    X(PrPr08gap22, "poiPrN08 {2} poiPrP08 {-2}", "PrPr08gap22")             \
    X(PiPi08gap32, "poiPiN08 {3} poiPiP08 {-3}", "PiPi08gap32")             \
    X(KaKa08gap32, "poiKaN08 {3} poiKaP08 {-3}", "KaKa08gap32")             \
    X(PrPr08gap32, "poiPrN08 {3} poiPrP08 {-3}", "PrPr08gap32")             \
    X(Pion0gap22a, "poiPiN {2} refP {-2}", "Pion0gap22a")                   \
    X(Pion0gap22b, "poiPiP {2} refN {-2}", "Pion0gap22b")                   \
    X(Kaon0gap22a, "poiKaN {2} refP {-2}", "Kaon0gap22a")                   \
    X(Kaon0gap22b, "poiKaP {2} refN {-2}", "Kaon0gap22b")                   \
    X(Prot0gap22a, "poiPrN {2} refP {-2}", "Prot0gap22a")                   \
    X(Prot0gap22b, "poiPrP {2} refN {-2}", "Prot0gap22b")

namespace GFWUtils
{
    enum class GFWRegion
    {
#define X(name) name,
        GFW_REGION_LIST
#undef X
    };

    inline std::string GetGFWRegion(GFWRegion region)
    {
        switch (region)
        {
#define X(name)           \
    case GFWRegion::name: \
        return #name;
            GFW_REGION_LIST
#undef X
        default:
            return "";
        }
    }

    enum class GFWConfig
    {
#define X(enum_name, cuts, name) enum_name,
        GFW_CONFIG_LIST
#undef X
    };

    inline std::pair<std::string, std::string> GetGFWConfig(GFWConfig cfg)
    {
        switch (cfg)
        {
#define X(enum_name, cuts, name) \
    case GFWConfig::enum_name:   \
        return {cuts, name};
            GFW_CONFIG_LIST
#undef X
        default:
            return {"", ""};
        }
    }

    double GetCumulant(GFW *gfw, GFWConfig config)
    {
        std::cout << GFWUtils::GetGFWConfig(config).first << GFWUtils::GetGFWConfig(config).second;
        return gfw->Calculate(gfw->GetCorrelatorConfig(
                                  GFWUtils::GetGFWConfig(config).first, GFWUtils::GetGFWConfig(config).second),
                              0, false)
            .real();
    }

    double GetNpair(GFW *gfw, GFWConfig config)
    {
        return gfw->Calculate(gfw->GetCorrelatorConfig(
                                  GFWUtils::GetGFWConfig(config).first, GFWUtils::GetGFWConfig(config).second),
                              0, true)
            .real();
    }
}

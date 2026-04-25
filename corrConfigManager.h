/**
 * @file corrConfigManager.h
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2026-04-21
 * @example
 *  GFW *gfw = new GFW();
    CorrConfigManager mgr(gfw);
    auto &cfg_ref = mgr.Get(CorrType::Ref08Gap22);
    double cum = gfw->Calculate(cfg_ref, 0, false).real();
    double npair = gfw->Calculate(cfg_ref, 0, true).real();
 *
 * @copyright Copyright (c) 2026
 *
 */

#ifndef CORRCONFIGMANAGER_H
#define CORRCONFIGMANAGER_H

#include "PWGCF/GenericFramework/Core/GFW.h"
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// ==================== 1. 定义枚举：彻底替代魔法数字 (共41个) ====================
enum class CorrType
{
    // --- Reference 流 (0-4) ---
    Ref08Gap22,
    Ref0Gap24,
    Ref0Gap22,
    Ref08Gap32,
    Ref08Gap34,

    // --- 粒子鉴别：0.8 快度区，2阶谐波 (5-10) ---
    Pion08Gap22a,
    Pion08Gap22b,
    Kaon08Gap22a,
    Kaon08Gap22b,
    Prot08Gap22a,
    Prot08Gap22b,

    // --- 粒子鉴别：全快度区，4阶谐波 (11-16) ---
    Pion0Gap24a,
    Pion0Gap24b,
    Kaon0Gap24a,
    Kaon0Gap24b,
    Prot0Gap24a,
    Prot0Gap24b,

    // --- 粒子鉴别：0.8 快度区，3阶谐波 (17-22) ---
    Pion08Gap32a,
    Pion08Gap32b,
    Kaon08Gap32a,
    Kaon08Gap32b,
    Prot08Gap32a,
    Prot08Gap32b,

    // --- 粒子鉴别：全快度区，6阶谐波 (23-28) ---
    Pion0Gap34a,
    Pion0Gap34b,
    Kaon0Gap34a,
    Kaon0Gap34b,
    Prot0Gap34a,
    Prot0Gap34b,

    // --- 相同粒子关联 (29-34) ---
    PiPi08Gap22,
    KaKa08Gap22,
    PrPr08Gap22,
    PiPi08Gap32,
    KaKa08Gap32,
    PrPr08Gap32,

    // --- 粒子鉴别：全快度区，2阶谐波 (35-40) ---
    Pion0Gap22a_Full,
    Pion0Gap22b_Full, // 注意：加了_Full避免与上面重名
    Kaon0Gap22a_Full,
    Kaon0Gap22b_Full,
    Prot0Gap22a_Full,
    Prot0Gap22b_Full,

    kTotal // 标记总数，不要删
};

// ==================== 2. 定义元数据结构体 ====================
struct CorrMeta
{
    CorrType type;
    const char *config_str; // GFW 配置字符串
    const char *name;       // 配置名字
    // 分类信息，方便后续批量处理
    enum class Particle
    {
        kRef,
        kPi,
        kKa,
        kPr
    };
    Particle particle;
    int harmonic; // 阶数 (2, 3)
};

struct Mask
{
    static constexpr int kRef = 1 << 0;           // 1
    static constexpr int kPion = 1 << 1;          // 2
    static constexpr int kKaon = 1 << 2;          // 4
    static constexpr int kProton = 1 << 3;        // 8
    static constexpr int kPionOverlap = 1 << 4;   // 16
    static constexpr int kKaonOverlap = 1 << 5;   // 32
    static constexpr int kProtonOverlap = 1 << 6; // 64
};

// ==================== 3. 管理器类 (核心) ====================
class CorrConfigManager
{
public:
    /**
     * @brief 构造函数：一键完成 Region 初始化 和 CorrConfig 生成
     * @param fGFW 指针
     * @param etaGap 默认 0.4
     * @param etaMax 默认 0.8
     */
    CorrConfigManager(GFW *fGFW, double etaGap = 0.4, double etaMax = 0.8)
        : m_gfw(fGFW)
    {
        // 第一步：初始化所有 Region
        InitRegions(etaGap, etaMax);
        // 第二步：初始化所有 CorrConfig
        InitConfigs();
        fGFW->CreateRegions();
    }

    // --- 核心访问接口 1：按枚举获取配置 (最推荐，编译期检查) ---
    const GFW::CorrConfig &Get(CorrType type) const
    {
        return m_corrconfigs[static_cast<size_t>(type)];
    }

    // --- 核心访问接口 2：按名字字符串获取配置 ---
    const GFW::CorrConfig &Get(const std::string &name) const
    {
        auto it = m_name_to_index.find(name);
        if (it == m_name_to_index.end())
        {
            throw std::runtime_error("CorrConfig name not found: " + name);
        }
        return m_corrconfigs[it->second];
    }

    // --- 批量处理接口：获取某一类粒子的所有配置 ---
    std::vector<CorrType> GetByParticle(CorrMeta::Particle p) const
    {
        std::vector<CorrType> res;
        for (const auto &meta : m_metas)
        {
            if (meta.particle == p)
                res.push_back(meta.type);
        }
        return res;
    }

    size_t Size() const { return m_corrconfigs.size(); }

private:
    GFW *m_gfw;
    std::vector<GFW::CorrConfig> m_corrconfigs;
    std::vector<CorrMeta> m_metas;
    std::unordered_map<std::string, size_t> m_name_to_index;

    void InitRegions(double etaGap, double etaMax)
    {
        // --- Reference ---
        m_gfw->AddRegion("reffull", -etaMax, etaMax, 1, Mask::kRef);
        m_gfw->AddRegion("refN08", -etaMax, -etaGap, 1, Mask::kRef);
        m_gfw->AddRegion("refP08", etaGap, etaMax, 1, Mask::kRef);
        m_gfw->AddRegion("refN", -etaMax, 0, 1, Mask::kRef);
        m_gfw->AddRegion("refP", 0, etaMax, 1, Mask::kRef);

        // --- Pions ---
        m_gfw->AddRegion("poiPiN08", -etaMax, -etaGap, 1, Mask::kPion);
        m_gfw->AddRegion("poiPiP08", etaGap, etaMax, 1, Mask::kPion);
        m_gfw->AddRegion("poiPiN", -etaMax, 0, 1, Mask::kPion);
        m_gfw->AddRegion("poiPiP", 0, etaMax, 1, Mask::kPion);

        // --- Overlap Pions ---
        m_gfw->AddRegion("olPiN", -etaMax, 0, 1, Mask::kPionOverlap);
        m_gfw->AddRegion("olPiP", 0, etaMax, 1, Mask::kPionOverlap);

        // --- Kaons ---
        m_gfw->AddRegion("poiKaN08", -etaMax, -etaGap, 1, Mask::kKaon);
        m_gfw->AddRegion("poiKaP08", etaGap, etaMax, 1, Mask::kKaon);
        m_gfw->AddRegion("poiKaN", -etaMax, 0, 1, Mask::kKaon);
        m_gfw->AddRegion("poiKaP", 0, etaMax, 1, Mask::kKaon);

        // --- Overlap Kaons ---
        m_gfw->AddRegion("olKaN", -etaMax, 0, 1, Mask::kKaonOverlap);
        m_gfw->AddRegion("olKaP", 0, etaMax, 1, Mask::kKaonOverlap);

        // --- Protons ---
        m_gfw->AddRegion("poiPrN08", -etaMax, -etaGap, 1, Mask::kProton);
        m_gfw->AddRegion("poiPrP08", etaGap, etaMax, 1, Mask::kProton);
        m_gfw->AddRegion("poiPrN", -etaMax, 0, 1, Mask::kProton);
        m_gfw->AddRegion("poiPrP", 0, etaMax, 1, Mask::kProton);

        // --- Overlap Protons ---
        m_gfw->AddRegion("olPrN", -etaMax, 0, 1, Mask::kProtonOverlap);
        m_gfw->AddRegion("olPrP", 0, etaMax, 1, Mask::kProtonOverlap);
    }

    void InitConfigs()
    {
        // 在这里集中定义所有 41 个配置！
        // 格式：{ 枚举类型, 配置字符串, 名字, 粒子类型, 阶数 }
        using P = CorrMeta::Particle;
        m_metas = {
            // --- Reference (0-4) ---
            {CorrType::Ref08Gap22, "refP08 {2} refN08 {-2}", "Ref08Gap22", P::kRef, 2},
            {CorrType::Ref0Gap24, "refN {2 2} refP {-2 -2}", "Ref0Gap24", P::kRef, 2},
            {CorrType::Ref0Gap22, "refN {2} refP {-2}", "Ref0Gap22", P::kRef, 2},
            {CorrType::Ref08Gap32, "refP08 {3} refN08 {-3}", "Ref08Gap32", P::kRef, 3},
            {CorrType::Ref08Gap34, "refP08 {3 3} refN08 {-3 -3}", "Ref08Gap34", P::kRef, 3},

            // --- 0.8 eta, 2nd harmonic (5-10) ---
            {CorrType::Pion08Gap22a, "poiPiN08 {2} refP08 {-2}", "Pion08gap22a", P::kPi, 2},
            {CorrType::Pion08Gap22b, "poiPiP08 {2} refN08 {-2}", "Pion08gap22b", P::kPi, 2},
            {CorrType::Kaon08Gap22a, "poiKaN08 {2} refP08 {-2}", "Kaon08gap22a", P::kKa, 2},
            {CorrType::Kaon08Gap22b, "poiKaP08 {2} refN08 {-2}", "Kaon08gap22b", P::kKa, 2},
            {CorrType::Prot08Gap22a, "poiPrN08 {2} refP08 {-2}", "Prot08gap22a", P::kPr, 2},
            {CorrType::Prot08Gap22b, "poiPrP08 {2} refN08 {-2}", "Prot08gap22b", P::kPr, 2},

            // --- Full eta, 4th harmonic (11-16) ---
            {CorrType::Pion0Gap24a, "poiPiN refN | olPiN {2 2} refP {-2 -2}", "Pion0gap24a", P::kPi, 2},
            {CorrType::Pion0Gap24b, "poiPiP refP | olPiP {2 2} refN {-2 -2}", "Pion0gap24b", P::kPi, 2},
            {CorrType::Kaon0Gap24a, "poiKaN refN | olKaN {2 2} refP {-2 -2}", "Kaon0gap24a", P::kKa, 2},
            {CorrType::Kaon0Gap24b, "poiKaP refP | olKaP {2 2} refN {-2 -2}", "Kaon0gap24b", P::kKa, 2},
            {CorrType::Prot0Gap24a, "poiPrN refN | olPrN {2 2} refP {-2 -2}", "Prot0gap24a", P::kPr, 2},
            {CorrType::Prot0Gap24b, "poiPrP refP | olPrP {2 2} refN {-2 -2}", "Prot0gap24b", P::kPr, 2},

            // --- 0.8 eta, 3rd harmonic (17-22) ---
            {CorrType::Pion08Gap32a, "poiPiN08 {3} refP08 {-3}", "Pion08gap32a", P::kPi, 3},
            {CorrType::Pion08Gap32b, "poiPiP08 {3} refN08 {-3}", "Pion08gap32b", P::kPi, 3},
            {CorrType::Kaon08Gap32a, "poiKaN08 {3} refP08 {-3}", "Kaon08gap32a", P::kKa, 3},
            {CorrType::Kaon08Gap32b, "poiKaP08 {3} refN08 {-3}", "Kaon08gap32b", P::kKa, 3},
            {CorrType::Prot08Gap32a, "poiPrN08 {3} refP08 {-3}", "Prot08gap32a", P::kPr, 3},
            {CorrType::Prot08Gap32b, "poiPrP08 {3} refN08 {-3}", "Prot08gap32b", P::kPr, 3},

            // --- Full eta, 6th harmonic (23-28) ---
            {CorrType::Pion0Gap34a, "poiPiN refN | olPiN {3 3} refP {-3 -3}", "Pion0gap34a", P::kPi, 3},
            {CorrType::Pion0Gap34b, "poiPiP refP | olPiP {3 3} refN {-3 -3}", "Pion0gap34b", P::kPi, 3},
            {CorrType::Kaon0Gap34a, "poiKaN refN | olKaN {3 3} refP {-3 -3}", "Kaon0gap34a", P::kKa, 3},
            {CorrType::Kaon0Gap34b, "poiKaP refP | olKaP {3 3} refN {-3 -3}", "Kaon0gap34b", P::kKa, 3},
            {CorrType::Prot0Gap34a, "poiPrN refN | olPrN {3 3} refP {-3 -3}", "Prot0gap34a", P::kPr, 3},
            {CorrType::Prot0Gap34b, "poiPrP refP | olPrP {3 3} refN {-3 -3}", "Prot0gap34b", P::kPr, 3},

            // --- Same particle correlations (29-34) ---
            {CorrType::PiPi08Gap22, "poiPiN08 {2} poiPiP08 {-2}", "PiPi08gap22", P::kPi, 2},
            {CorrType::KaKa08Gap22, "poiKaN08 {2} poiKaP08 {-2}", "KaKa08gap22", P::kKa, 2},
            {CorrType::PrPr08Gap22, "poiPrN08 {2} poiPrP08 {-2}", "PrPr08gap22", P::kPr, 2},
            {CorrType::PiPi08Gap32, "poiPiN08 {3} poiPiP08 {-3}", "PiPi08gap32", P::kPi, 3},
            {CorrType::KaKa08Gap32, "poiKaN08 {3} poiKaP08 {-3}", "KaKa08gap32", P::kKa, 3},
            {CorrType::PrPr08Gap32, "poiPrN08 {3} poiPrP08 {-3}", "PrPr08gap32", P::kPr, 3},

            // --- Full eta, 2nd harmonic (35-40) ---
            {CorrType::Pion0Gap22a_Full, "poiPiN {2} refP {-2}", "Pion0gap22a", P::kPi, 2},
            {CorrType::Pion0Gap22b_Full, "poiPiP {2} refN {-2}", "Pion0gap22b", P::kPi, 2},
            {CorrType::Kaon0Gap22a_Full, "poiKaN {2} refP {-2}", "Kaon0gap22a", P::kKa, 2},
            {CorrType::Kaon0Gap22b_Full, "poiKaP {2} refN {-2}", "Kaon0gap22b", P::kKa, 2},
            {CorrType::Prot0Gap22a_Full, "poiPrN {2} refP {-2}", "Prot0gap22a", P::kPr, 2},
            {CorrType::Prot0Gap22b_Full, "poiPrP {2} refN {-2}", "Prot0gap22b", P::kPr, 2},
        };

        // 自动生成 corrconfigs 和 name map
        m_corrconfigs.reserve(m_metas.size());
        for (size_t i = 0; i < m_metas.size(); ++i)
        {
            const auto &meta = m_metas[i];
            m_corrconfigs.push_back(m_gfw->GetCorrelatorConfig(meta.config_str, meta.name, false));
            m_name_to_index[meta.name] = i;
        }
    }
};

#endif // CORRCONFIGMANAGER_H

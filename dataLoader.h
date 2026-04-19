#ifndef FILELOADER
#define FILELOADER
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TString.h"
#include "dataFrame/track.h"
#include <iostream>
#include <memory> // 用于智能指针
#include <stdexcept>
#include <vector>

/**
 * @brief 该类用于读取所有的粒子数据，注意到对其iter会逐个返回一条track，不区分event
 *
 */
class AMPTDataReader
{
protected:
    std::unique_ptr<TChain> chain; // 使用智能指针管理内存

    /**
     * @brief iterator for data reader
     *
     */
    class Iterator
    {
    private:
        AMPTDataReader *m_tree;
        Long64_t m_index;

    public:
        Iterator(AMPTDataReader *tree, Long64_t index) : m_tree(tree), m_index(index) {}

        // 返回const引用避免拷贝（如果Track允许）
        const Track &operator*() const
        {
            // 这里需要将当前Track缓存为成员变量，或者确保Track可以被引用返回
            // 为了保持原有逻辑，这里仍返回值，但建议优化Track的设计
            return *new Track{
                m_tree->eventID,
                m_tree->nParticles,
                m_tree->imp,

                m_tree->pdgPid,
                m_tree->p_x, m_tree->p_y, m_tree->p_z};
        }

        // 前置++：读取下一个事例
        Iterator &operator++()
        {
            m_index++;
            if (m_index < m_tree->chain->GetEntries())
            {
                m_tree->chain->GetEntry(m_index);
            }
            return *this;
        }

        // 比较迭代器
        bool operator!=(const Iterator &other) const
        {
            return m_index != other.m_index;
        }
    };

public:
    int eventID, nParticles, pdgPid;
    double p_x, p_y, p_z, imp;

    /**
     * @brief Construct a new AMPTDataReader object
     * @example
     * AMPTDataReader tree("/home/huinaibing/ampt_result/cent50-60/Result/ampt_19370820_", 100);
     * for (auto &evt : tree)
     * {
     *     std::cout << "Event " << evt.pdgPid << std::endl;
     * }
     *
     * @param filePrefix
     * @param nFiles
     * @param fileSuffix
     */
    AMPTDataReader(const char *filePrefix, int nFiles)
        : chain(std::make_unique<TChain>("particles")) // 初始化智能指针
    {
        for (int i = 0; i < nFiles; i++)
        {
            TString fname = TString::Format("%s%d.root", filePrefix, i);
            if (chain->Add(fname, 0))
            { // 0=不立即加载所有数据，速度快
                std::cout << "✓ 加入文件: " << fname << std::endl;
            }
            else
            {
                std::cerr << "✗ 跳过不存在的文件: " << fname << std::endl;
            }
        }
        std::cout << "\n总共加载 " << chain->GetEntries() << " 个事例" << std::endl;

        // 设置分支地址
        chain->SetBranchAddress("eventID", &eventID);
        chain->SetBranchAddress("nParticles", &nParticles);
        chain->SetBranchAddress("imp", &imp);
        chain->SetBranchAddress("pdgPid", &pdgPid);
        chain->SetBranchAddress("p_x", &p_x);
        chain->SetBranchAddress("p_y", &p_y);
        chain->SetBranchAddress("p_z", &p_z);

        // 只有在有数据时才读取第一个entry
        if (chain->GetEntries() > 0)
        {
            chain->GetEntry(0);
        }
    }

    // 禁用拷贝构造和拷贝赋值（避免double free）
    AMPTDataReader(const AMPTDataReader &) = delete;
    AMPTDataReader &operator=(const AMPTDataReader &) = delete;

    // 启用移动构造和移动赋值（可选，提高性能）
    AMPTDataReader(AMPTDataReader &&) = default;
    AMPTDataReader &operator=(AMPTDataReader &&) = default;

    Iterator begin() { return Iterator(this, 0); }
    Iterator end() { return Iterator(this, chain->GetEntries()); }
    Long64_t size() const { return chain->GetEntries(); }

    ~AMPTDataReader() = default; // 智能指针会自动释放chain
};

#endif

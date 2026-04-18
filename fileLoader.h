#ifndef FILELOADER
#define FILELOADER
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TString.h"
#include <iostream>
#include <stdexcept>
#include <vector>

class AMPTDataReader
{
protected:
    TChain *chain;

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

        // 解引用：返回当前树对象（可直接访问变量）
        AMPTDataReader &operator*() { return *m_tree; }
        AMPTDataReader *operator->() { return m_tree; }

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
    AMPTDataReader(const char *filePrefix, int nFiles, const char *fileSuffix = ".root")
    {
        chain = new TChain("particles"); // 你的Tree名字是particles
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

        chain->SetBranchAddress("eventID", &eventID);
        chain->SetBranchAddress("nParticles", &nParticles);
        chain->SetBranchAddress("imp", &imp);
        chain->SetBranchAddress("pdgPid", &pdgPid);
        chain->SetBranchAddress("p_x", &p_x);
        chain->SetBranchAddress("p_y", &p_y);
        chain->SetBranchAddress("p_z", &p_z);

        chain->GetEntry(0);
    }

    Iterator begin() { return Iterator(this, 0); }
    Iterator end() { return Iterator(this, chain->GetEntries()); }
    Long64_t size() const { return chain->GetEntries(); }

    AMPTDataReader() { delete chain; }
};

#endif

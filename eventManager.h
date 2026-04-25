/**
 * @file eventManager.h
 * @author huinaibing
 * @brief 用于将track按event分类
 * @version 0.1
 * @date 2026-04-21
 * @example
 *    AMPTEventReader cent_30_40("/home/huinaibing/ampt_result/cent30-40/Result/ampt_19367451_", 100);
 *    for (const auto &evt : cent_30_40)
      {
        gfw->Clear();
        for (const auto &trk : evt.particles)
        {
            gfw->Fill(trk.GetEta(), 0, trk.GetPhi(), 1, 1);
        }
    }
 *
 * @copyright Copyright (c) 2026
 *
 */

#ifndef EVENTMANAGER
#define EVENTMANAGER
#include "dataFrame/event.h"
#include "dataLoader.h"
#include <iostream>
#include <memory> // 引入智能指针
#include <vector>

class AMPTEventReader
{
private:
    std::unique_ptr<AMPTDataReader> fReader;
    std::vector<Event> fEvents;

public:
    class Iterator
    {
    private:
        AMPTEventReader *m_parent;
        size_t m_index;

    public:
        Iterator(AMPTEventReader *parent, size_t index) : m_parent(parent), m_index(index) {}

        // 返回const引用避免拷贝（如果Event允许）
        const Event &operator*() const { return m_parent->fEvents[m_index]; }

        Iterator &operator++()
        {
            m_index++;
            return *this;
        }
        bool operator!=(const Iterator &other) const { return m_index != other.m_index; }
    };

    AMPTEventReader(const char *filePrefix, int nFiles) : fReader(std::make_unique<AMPTDataReader>(filePrefix, nFiles)) // 初始化智能指针
    {
        int eventTmpId = -1;
        for (const auto &track : *fReader)
        {
            if (track.eventID != eventTmpId)
            {
                // 新事件：直接构造Event并添加track
                eventTmpId = track.eventID;
                fEvents.emplace_back(); // 直接在vector末尾构造Event（比push_back更高效）
                fEvents.back().particles.emplace_back(track);
            }
            else
            {
                // 同一事件：添加到当前Event的particles中
                fEvents.back().particles.emplace_back(track);
            }
        }

        eventTmpId = 0;
        for (auto &evt : this->fEvents)
        {
            evt.imp = evt.particles.back().imp;
            evt.eventID = eventTmpId++;
        }
    }

    AMPTEventReader(const AMPTEventReader &) = delete;
    AMPTEventReader &operator=(const AMPTEventReader &) = delete;
    AMPTEventReader(AMPTEventReader &&) = default;
    AMPTEventReader &operator=(AMPTEventReader &&) = default;

    const Event &GetEvent(size_t m_index) const
    {
        if (m_index >= fEvents.size())
        {
            throw std::out_of_range("Event index out of range");
        }
        return fEvents[m_index];
    }

    // 迭代器接口
    Iterator begin() { return Iterator(this, 0); }
    Iterator end() { return Iterator(this, fEvents.size()); }
    size_t size() const { return fEvents.size(); }
};
#endif

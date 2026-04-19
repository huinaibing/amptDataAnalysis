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
    AMPTDataReader *fReader;    // 用智能指针管理AMPTDataReader
    std::vector<Event> fEvents; // 直接用vector对象，避免指针管理

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

    AMPTEventReader(const char *filePrefix, int nFiles)
    {
        this->fReader = new AMPTDataReader(filePrefix, nFiles);
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
            evt.nParticles = evt.particles.size();
            evt.imp = evt.particles.back().imp;
            evt.eventID = eventTmpId++;
        }
    }

    // 补充GetEvent实现（返回const引用更高效）
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

    // 不需要显式析构函数，智能指针和vector会自动管理内存
};
#endif

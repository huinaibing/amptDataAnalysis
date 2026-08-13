#ifndef AMPT_DATA_ANALYSIS_EVENT_MANAGER_H
#define AMPT_DATA_ANALYSIS_EVENT_MANAGER_H

#include "dataFrame/event.h"
#include "dataLoader.h"

#include <cstddef>

/**
 * Stream AMPT particles one event at a time. Memory use is proportional to
 * the largest event rather than to all events in the selected ROOT files.
 */
class AMPTEventReader {
public:
  class Iterator {
  public:
    Iterator() = default;

    Iterator(AMPTDataReader *reader, bool atEnd)
        : mCurrent(reader ? (atEnd ? reader->end() : reader->begin())
                          : AMPTDataReader::Iterator{}),
          mEnd(reader ? reader->end() : AMPTDataReader::Iterator{}),
          mAtEnd(atEnd || !reader || mCurrent == mEnd) {
      if (!mAtEnd) {
        readNextEvent();
      }
    }

    const Event &operator*() const { return mEvent; }
    const Event *operator->() const { return &mEvent; }

    Iterator &operator++() {
      readNextEvent();
      return *this;
    }

    bool operator==(const Iterator &other) const {
      if (mAtEnd || other.mAtEnd) {
        return mAtEnd == other.mAtEnd;
      }
      return mCurrent == other.mCurrent;
    }

    bool operator!=(const Iterator &other) const { return !(*this == other); }

  private:
    void readNextEvent() {
      if (mCurrent == mEnd) {
        mAtEnd = true;
        mEvent = {};
        return;
      }

      const Track first = *mCurrent;
      const int sourceFile = first.sourceFile;
      const int sourceEventID = first.eventID;

      mEvent = {};
      mEvent.eventID = sourceEventID;
      mEvent.imp = first.imp;
      if (first.nParticles > 0) {
        mEvent.particles.reserve(static_cast<std::size_t>(first.nParticles));
      }

      while (mCurrent != mEnd) {
        const Track track = *mCurrent;
        if (track.sourceFile != sourceFile || track.eventID != sourceEventID) {
          break;
        }
        mEvent.imp = track.imp;
        mEvent.particles.emplace_back(track);
        ++mCurrent;
      }

      mAtEnd = false;
    }

    AMPTDataReader::Iterator mCurrent;
    AMPTDataReader::Iterator mEnd;
    Event mEvent;
    bool mAtEnd = true;
  };

  AMPTEventReader(const std::string &filePrefix, int nFiles)
      : mReader(filePrefix, nFiles) {}

  AMPTEventReader(const AMPTEventReader &) = delete;
  AMPTEventReader &operator=(const AMPTEventReader &) = delete;

  Iterator begin() { return Iterator(&mReader, false); }
  Iterator end() { return Iterator(&mReader, true); }
  Long64_t particleRows() const { return mReader.entries(); }

private:
  AMPTDataReader mReader;
};

#endif // AMPT_DATA_ANALYSIS_EVENT_MANAGER_H

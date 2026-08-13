#ifndef AMPT_DATA_ANALYSIS_DATA_LOADER_H
#define AMPT_DATA_ANALYSIS_DATA_LOADER_H

#include "dataFrame/track.h"

#include "TChain.h"
#include "TString.h"

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

/**
 * Read particle rows from an AMPT TChain without copying them into an
 * intermediate container. The iterator returns Track values because ROOT
 * reuses the branch buffers for every entry.
 */
class AMPTDataReader {
public:
  class Iterator {
  public:
    Iterator() = default;
    Iterator(AMPTDataReader *reader, Long64_t index)
        : mReader(reader), mIndex(index) {}

    Track operator*() const {
      return {mReader->mEventID,
              mReader->mNParticles,
              mReader->mImpactParameter,
              mReader->mPdgPid,
              mReader->mPx,
              mReader->mPy,
              mReader->mPz,
              mReader->mChain->GetTreeNumber()};
    }

    Iterator &operator++() {
      ++mIndex;
      if (mIndex < mReader->entries()) {
        mReader->loadEntry(mIndex);
      }
      return *this;
    }

    bool operator==(const Iterator &other) const {
      return mReader == other.mReader && mIndex == other.mIndex;
    }

    bool operator!=(const Iterator &other) const { return !(*this == other); }

  private:
    AMPTDataReader *mReader = nullptr;
    Long64_t mIndex = 0;
  };

  AMPTDataReader(const std::string &filePrefix, int nFiles)
      : mChain(std::make_unique<TChain>("particles")) {
    if (nFiles < 0) {
      throw std::invalid_argument("nFiles must not be negative");
    }

    int addedFiles = 0;
    for (int i = 0; i < nFiles; ++i) {
      const TString fileName =
          TString::Format("%s%d.root", filePrefix.c_str(), i);
      if (mChain->Add(fileName, 0) > 0) {
        ++addedFiles;
        std::cout << "Added input file: " << fileName << '\n';
      } else {
        std::cerr << "Skipped missing input file: " << fileName << '\n';
      }
    }

    bindRequiredBranch("eventID", &mEventID);
    bindRequiredBranch("nParticles", &mNParticles);
    bindRequiredBranch("imp", &mImpactParameter);
    bindRequiredBranch("pdgPid", &mPdgPid);
    bindRequiredBranch("p_x", &mPx);
    bindRequiredBranch("p_y", &mPy);
    bindRequiredBranch("p_z", &mPz);

    std::cout << "Loaded " << addedFiles << " files with " << entries()
              << " particle rows\n";
    if (entries() > 0) {
      loadEntry(0);
    }
  }

  AMPTDataReader(const AMPTDataReader &) = delete;
  AMPTDataReader &operator=(const AMPTDataReader &) = delete;
  AMPTDataReader(AMPTDataReader &&) = delete;
  AMPTDataReader &operator=(AMPTDataReader &&) = delete;

  Iterator begin() {
    if (entries() > 0) {
      loadEntry(0);
    }
    return Iterator(this, 0);
  }

  Iterator end() { return Iterator(this, entries()); }
  Long64_t entries() const { return mChain->GetEntries(); }

private:
  template <typename T> void bindRequiredBranch(const char *name, T *address) {
    if (!mChain->GetBranch(name)) {
      throw std::runtime_error(std::string("Missing branch '") + name +
                               "' in particles tree");
    }
    if (mChain->SetBranchAddress(name, address) < 0) {
      throw std::runtime_error(std::string("Cannot bind branch '") + name +
                               "'");
    }
  }

  void loadEntry(Long64_t index) {
    if (mChain->GetEntry(index) <= 0) {
      throw std::runtime_error("Cannot read particle entry " +
                               std::to_string(index));
    }
  }

  std::unique_ptr<TChain> mChain;
  int mEventID = 0;
  int mNParticles = 0;
  int mPdgPid = 0;
  double mPx = 0.;
  double mPy = 0.;
  double mPz = 0.;
  double mImpactParameter = 0.;
};

#endif // AMPT_DATA_ANALYSIS_DATA_LOADER_H

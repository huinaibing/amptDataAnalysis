#include "../analysisUtils.h"

#include <iostream>
#include <stdexcept>

void test_fileLoader(
    const char *filePrefix =
        "/home/huinaibing/ampt_result/cent50-60/Result1/ampt_19370820_",
    int nFiles = 1) {
  const auto &config = ampt_analysis::defaultConfig();
  AMPTEventReader events(filePrefix, nFiles);

  Long64_t eventCount = 0;
  Long64_t particleCount = 0;
  for (const auto &event : events) {
    const auto moments = ampt_analysis::getMeanPtMoments(event, 0, config);
    if (moments.count > event.particles.size()) {
      throw std::runtime_error("Selected particle count exceeds event size");
    }
    ++eventCount;
    particleCount += static_cast<Long64_t>(event.particles.size());
  }

  if (eventCount == 0 || particleCount != events.particleRows()) {
    throw std::runtime_error("AMPT event streaming lost particle rows");
  }

  std::cout << "Validated " << eventCount << " events and " << particleCount
            << " particle rows." << std::endl;
}

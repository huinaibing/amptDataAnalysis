#include "../external/o2physics/FlowContainerRoot.h"

#include "TFileMerger.h"
#include "TClass.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

// Load FlowContainer before merging so its nominal and bootstrap profiles
// are combined by FlowContainer::Merge, including with a plain ROOT install.
void merge_results(const char *outputFile, const char *inputListFile) {
  // ClassDefInline supplies I/O but does not register a custom Merge callback.
  FlowContainer::Class()->SetMerge(
      [](void *object, TCollection *inputs, TFileMergeInfo *) -> Long64_t {
        return static_cast<FlowContainer *>(object)->Merge(inputs);
      });

  std::ifstream input(inputListFile);
  if (!input) {
    throw std::runtime_error(std::string("Cannot open merge input list: ") +
                             inputListFile);
  }

  TFileMerger merger(false, false);
  if (!merger.OutputFile(outputFile, "RECREATE")) {
    throw std::runtime_error(std::string("Cannot create merged output: ") +
                             outputFile);
  }
  std::string file;
  int count = 0;
  while (std::getline(input, file)) {
    if (file.empty()) {
      continue;
    }
    if (!merger.AddFile(file.c_str())) {
      throw std::runtime_error("Cannot add merge input: " + file);
    }
    ++count;
  }
  if (count == 0 || !merger.Merge()) {
    throw std::runtime_error("Failed to merge analysis results");
  }
  std::cout << "Merged " << count << " files into " << outputFile << '\n';
}

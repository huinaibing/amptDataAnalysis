#ifndef AMPT_DATA_ANALYSIS_CENTRALITY_CONFIG_H
#define AMPT_DATA_ANALYSIS_CENTRALITY_CONFIG_H

#include "json.hpp"

#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct CentralityConfig {
  std::string path;
  double bin_val;
  int n_files;
};

class CentConfigReader {
public:
  static std::vector<CentralityConfig> load(const std::string &json_path) {
    std::ifstream f(json_path);
    if (!f.is_open()) {
      throw std::runtime_error("无法打开配置文件: " + json_path);
    }

    nlohmann::json document;
    f >> document;
    if (!document.is_array()) {
      throw std::runtime_error("配置文件顶层必须是数组: " + json_path);
    }

    std::vector<CentralityConfig> result;
    result.reserve(document.size());
    for (const auto &item : document) {
      CentralityConfig cfg{item.at("path").get<std::string>(),
                           item.at("bin_val").get<double>(),
                           item.at("n_files").get<int>()};
      if (cfg.path.empty() || cfg.n_files < 0) {
        throw std::runtime_error("配置包含空路径或负数 n_files: " + json_path);
      }
      result.emplace_back(std::move(cfg));
    }
    return result;
  }
};

#endif // AMPT_DATA_ANALYSIS_CENTRALITY_CONFIG_H

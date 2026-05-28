#ifndef FILENAME
#define FILENAME

#include "json.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

struct CentralityConfig
{
    std::string path;
    double bin_val;
    int n_files; // 新增：文件数量
};

using json = nlohmann::json;

// JSON 配置读取类
class CentConfigReader
{
public:
    // 核心：读取JSON → 返回 vector<CentralityConfig>
    static std::vector<CentralityConfig> load(const std::string &json_path)
    {
        std::vector<CentralityConfig> result;
        std::ifstream f(json_path);

        if (!f.is_open())
            throw std::runtime_error("无法打开配置文件: " + json_path);

        json j;
        f >> j;

        for (const auto &item : j)
        {
            CentralityConfig cfg;

            // 重点：string 转 const char*
            cfg.path = item["path"];
            cfg.bin_val = item["bin_val"].get<double>();
            cfg.n_files = item["n_files"].get<int>();

            result.push_back(cfg);
        }

        return result;
    }
};

#endif

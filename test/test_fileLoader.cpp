#include "../fileLoader.h"

void test_fileLoader()
{
    AMPTDataReader tree("/home/huinaibing/ampt_result/cent50-60/Result/ampt_19370820_", 100);

    for (auto &evt : tree)
    {
        //std::cout << "Event " << evt.pdgPid << std::endl;
    }
}

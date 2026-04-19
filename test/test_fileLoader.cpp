#include "../dataLoader.h"
#include "../eventManager.h"

void test_fileLoader()
{
    // AMPTDataReader tree("/home/huinaibing/ampt_result/cent50-60/Result/ampt_19370820_", 100);

    // for (const auto &evt : tree)
    // {
    //     std::cout << "Event " << evt.eventID << " " << evt.p_x << std::endl;
    // }

    AMPTEventReader events("/home/huinaibing/ampt_result/cent50-60/Result/ampt_19370820_", 100);
    for (const auto &evt : events)
    {
        std::cout << "meanpt " << evt.GetMeanPt() << std::endl;
    }
}

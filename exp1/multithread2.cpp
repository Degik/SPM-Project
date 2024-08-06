#include <vector>
#include <future>
#include <iostream>
#include <thread>
#include "hpc_helpers.hpp"

uint64_t fibo(uint64_t number){
    uint64_t a_0 = 0;
    uint64_t a_1 = 1;
    for(uint64_t i = 0; i < number; i++){
        const uint64_t tmp = a_0;
        a_0 = a_1;
        a_1 += tmp;
    }
    return a_0;
}


int main(int argc, char * argv[]){

    TIMERSTART(init)
    const uint64_t nmumberThreads = 32;
    // create
    std::vector<std::future<uint64_t>> results;

    for (uint64_t idx = 0; idx < nmumberThreads; idx++){
        results.emplace_back(std::async(std::launch::async, fibo, idx));
    }

    for (auto& result: results){
        std::cout << "Risultato: " << result.get() << std::endl;
    }
    TIMERSTOP(init)
}
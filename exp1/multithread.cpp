#include <cstdint>
#include <iostream>
#include <thread>
#include <vector>
#include <future>
#include "hpc_helpers.hpp"

template <typename id_name>

void say_hello(id_name id){
    std::cout << "Hello from thread: " << id << std::endl;
}

template <typename value_t,
          typename idx_t>

void fibo(value_t number, std::promise<value_t> && result){
    value_t a_0 = 0;
    value_t a_1 = 1;
    for(idx_t idx = 0; idx < number; idx++){
        const value_t tmp = a_0;
        a_0 = a_1;
        a_1 += tmp;
    }

    result.set_value(a_0);
}


int main(int argc, char *argv[]){
    if (argc <= 1){
        std::cout << "give me the number of threads" << std::endl;
    }

    const uint64_t numThreads = std::stol(argv[1]);

    std::cout << "Number: " << numThreads << std::endl;

    // Promise and future
    std::vector<std::future<uint64_t>> results;
    std::vector<std::thread> threads;

    TIMERSTART(init)

    for(uint64_t id = 0; id < numThreads; id++){
        // Create future and promise

        std::promise<uint64_t> promise;
        results.emplace_back(promise.get_future());

        threads.emplace_back(
            fibo<uint64_t, uint64_t>, id, std::move(promise)
        );
    }

    for(auto& thread: threads){
        thread.join();
    }

    static uint64_t counter = 0;

    for(auto& res: results){
        counter++;
        std::cout << "Result for (" << counter << "): " << res.get() << std::endl;
    }

    TIMERSTOP(init)

}
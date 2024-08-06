#include <iostream>
#include <thread>
#include <future>
#include <vector>
#include <cstdint>


// Define generic data type
template <typename value_t,
          typename index_t>
// fibo function take n and promise as result
void fibo(value_t n, std::promise<value_t> && result){
    value_t a_0 = 0;
    value_t a_1 = 1;

    for (index_t index = 0; index < n; index++){
        const value_t tmp = a_0;
        a_0 = a_1;
        a_1 += tmp;
    }
    result.set_value(a_0);
}

int main(int argc, char* argv[]){
    // Number of threads
    const uint64_t number_threads = 32;
    // Vector of threads
    std::vector<std::thread> threads;
    // Vector with futures, each future have datatype uint64_t
    std::vector<std::future<uint64_t>> results;
    
    for (uint64_t id; id < number_threads; id++){
        // Create promise for the thread_i
        std::promise<uint64_t> promise;
        // Store associated future
        results.emplace_back(promise.get_future());
        // Create thread and call on it the fibo function
        // fibo take (value_t n, value_t promise)
        threads.emplace_back(
            fibo<uint64_t, uint64_t>, id, std::move(promise)
        );
    }

    for (auto& result: results){
        std::cout << result.get() << std::endl;
    }

    for (auto& thread: threads){
        thread.join();
    }
}
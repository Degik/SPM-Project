#include <vector>
#include <mutex>
#include <chrono>
#include <thread>
#include <future>
#include <iostream>
#include <condition_variable>
#include "hpc_helpers.hpp"

int main(int argc, char * argv[]){
    const uint64_t numberThreads = 32;

    std::mutex mutex;
    std::condition_variable cv;

    std::promise<void> promise;
    auto shared_future = promise.get_future().share();

    auto students = [&] (uint64_t myid) -> void {

        shared_future.get();
        std::cout << "Thread(" << myid << "): make break" << std::endl;
        
    };

    std::vector<std::thread> threadsStudents;
    
    for (uint64_t idx = 0; idx < numberThreads; idx++){
        threadsStudents.emplace_back(
            students, idx
        );
    }
    // Sleep for 2 seconds
    std::this_thread::sleep_for(std::chrono::seconds(2));

    promise.set_value();

    for (auto& thread: threadsStudents){
        thread.join();
    }

}
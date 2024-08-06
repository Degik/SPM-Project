#include <vector>
#include <thread>
#include <mutex>
#include <future>
#include <string>
#include <fstream>
#include "hpc_helpers.hpp"


template <typename value_t, typename index_t>
void dump_binary(const value_t * data,
                 const index_t length,
                 std::string filename){

                    std::ofstream ofile(filename.c_str(), std::ios::binary);
                    ofile.write((char*) data, sizeof(value_t)*length);
                    ofile.close();
                 }

template <typename value_t, typename index_t>
void load_binary(const value_t * data,
                 const index_t length,
                 std::string filename){

                    std::ifstream ifile(filename.c_str(), std::ios::binary);
                    ifile.read((char*) data, sizeof(value_t)*length);
                    ifile.close();
                 }

std::mutex mutex;
template <typename value_t, typename index_t>
void dynamic_all_pairs(std::vector<value_t>& mnist,
                       std::vector<value_t>& all_pair,
                       index_t rows,
                       index_t cols,
                       index_t nthreads,
                       index_t chunk_size=64/sizeof(value_t)){

                        index_t global_lower = 0;
                        auto dynamic_block_cyclic = [&] (const index_t& id) -> void {
                           index_t lower = 0;
                           while(lower < rows){

                              {
                                 std::lock_guard<std::mutex> guard(mutex);
                                 lower = global_lower;
                                 global_lower += chunk_size;
                              }

                              const index_t upper = std::min(lower+chunk_size, rows);

                              for (index_t i = lower; i < upper; i++){
                                 for (index_t I = 0; I <= i; I++){
                                    value_t accum = value_t(0);
                                    for (index_t j = 0; j < cols; j++){
                                       value_t residue = mnist[i*cols+j] - mnist[I*cols+j];
                                       accum += residue * residue;
                                    }

                                    all_pair[i*rows+I] = all_pair[I*rows+i] = accum;
                                 }
                              }
                           }
                        };

                        std::vector<std::thread> threads;

                        for (uint64_t idx = 0; idx < nthreads; idx++){
                           threads.emplace_back(
                              dynamic_block_cyclic, idx
                           );
                        }

                        for(auto& thread: threads){
                           thread.join();
                        }

                     }

int main(int argc, char * argv[]){
    typedef no_init_t<float> value_t;
    typedef uint64_t index_t;

    const index_t rows = 70000;
    const index_t cols = 28*28;
    const index_t length = rows*cols;

    TIMERSTART(load_data_from_disk);
    std::vector<value_t> mnist(length); // 70000 * (28*28)
    load_binary(mnist.data(), length, "mnist_dataset.bin");
    TIMERSTOP(load_data_from_disk);

    const index_t nthreads = 4;
    const index_t chunk_size = 1;

    TIMERSTART(compute_distance);
    std::vector<value_t> all_pair(length);
    dynamic_all_pairs(mnist, all_pair, rows, cols, nthreads, chunk_size);
    TIMERSTOP(compute_distance);
}
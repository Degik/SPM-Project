#include <cmath>
#include <mutex>
#include <thread>
#include <vector>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <condition_variable>

// FastFlow
#include <ff/ff.hpp>
#include <ff/farm.hpp>

//#define DEBUG

using vector_d = std::vector<double>;

const int MIN_CHUNK_SIZE = 512;

/*!
    \name DiagonalTask
    \brief DiagonalTask
    \note DiagonalTask - Calculate the current element in the diagonal
*/
struct DiagonalTask {
    int k;
    int m_start;
    int m_end;
    vector_d &M;
    uint16_t N;
    std::atomic<int> &tasks;
};

/*!
    \name FillMatrix
    \brief Fill the matrix M with the values
    \note Fill the matrix M with the values
*/
vector_d* FillMatrix(vector_d *M, uint16_t N) {
    for(int m = 0; m < N; m++) {
        (*M)[m*N+m] = static_cast<double>(m+1)/N;
    }
    return M;
}

/*!
    \name SaveMatrixPtrToFile
    \param M vector_d M
    \param N uint16_t N
    \param filename string filename
    \brief Save the matrix M to a file
    \note Save the matrix M to a file with the name filename
*/
void SaveMatrixToFile(vector_d *M, uint16_t N, std::string filename) {
    std::ofstream file;
    file.open(filename);
    file << std::fixed << std::setprecision(6);
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            file << (*M)[i*N+j] << " "; // M[i][j]
        }
        file << std::endl;
    }
    file.close();
}

/*!
    \name DiagonalWorker
    \brief DiagonalWorker
    \note DiagonalWorker - Calculate the current element in the diagonal
*/
struct DiagonalWorker : ff::ff_node_t<DiagonalTask> {
    DiagonalTask* svc(DiagonalTask *task) {
        for (int m = task->m_start; m < task->m_end; m++) { // Calculate the current element in the diagonal (Chunck block)
            int row = m * task->N;
            int col_t = (m + task->k) * task->N;
            double element = 0.0;
            for (int i = 0; i < task->k; i++) {
                element += task->M[row + i + m] * task->M[col_t + i + m + 1];
            }
            double new_element = std::cbrt(element);
            task->M[row + m + task->k] = new_element;
            task->M[col_t + m] = new_element;
        }
        task->tasks--; // Decrease the number of tasks
        delete task;
        return GO_ON;
    }
};

/*!
    \name DiagonalEmitter
    \brief DiagonalEmitter
    \note DiagonalEmitter - Emit the tasks for the workers to calculate the diagonal elements
*/
struct DiagonalEmitter : ff::ff_monode_t<int, DiagonalTask> {
    vector_d &M;
    uint16_t N;
    int W;

    DiagonalEmitter(vector_d &M, uint16_t N, int num_workers) : M(M), N(N), W(W) {}

    DiagonalTask* svc(int*) {
        for (int k = 1; k < N; k++) {
            int total_m = N - k;
            int num_tasks = std::min(W, total_m);                                           // Number of tasks
            int chunk_size = (total_m + num_tasks - 1) / num_tasks;                         // Chunk size for each task
            std::atomic<int> tasks(num_tasks);                                              // Atomic counter for the tasks           
            int m_start = 0, m_end = 0;
            for (int task_id = 0; task_id < num_tasks; task_id++) {                         // Create the tasks and send them to the workers
                m_start = task_id * chunk_size;                                             // Start index  
                m_end = std::min(m_start + chunk_size, total_m);                            // End index
                ff_send_out(new DiagonalTask{k, m_start, m_end, M, N, std::ref(tasks)});    // Send the task to the worker
            }

            // Wait for the tasks to finish
            while (tasks.load() > 0) {
                std::this_thread::yield();
            }
        }
        broadcast_task(EOS);
        return EOS;
    }
};

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " N (Size N*N) W (Workers)" << std::endl;
        return -1;
    }

    unsigned int max_threads = std::thread::hardware_concurrency();
    #ifdef DEBUG
        std::cout << "Max threads: " << max_threads << std::endl;
    #endif

    const uint16_t N = atoi(argv[1]);
    const uint16_t W = atoi(argv[2]);

    if(W > max_threads-1){
        std::cout << "The number of workers is higher than the number of threads available" << std::endl;
        return -1;
    }

    // Create and fill the matrix
    // Process to create the matrix
    auto start = std::chrono::high_resolution_clock::now();

    vector_d M = vector_d(N*N, 0.0);
    // Fill the matrix
    FillMatrix(&M, N);
    #ifdef DEBUG
        SaveMatrixToFile(&M, N, "matrix_farm_chunk_normal.txt");
    #endif
    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> passed_time = stop - start;
    std::cout << "Matrix created and filled in: " << passed_time.count() << " seconds" << std::endl;


    ff::ffTime(ff::START_TIME);

    std::vector<std::unique_ptr<ff::ff_node>> workers;
    for(int i = 0; i < W; i++) {
        workers.push_back(std::make_unique<DiagonalWorker>());
    }

    ff::ff_Farm<DiagonalTask> farm(std::move(workers));
    farm.remove_collector();
    farm.wrap_around();
    farm.set_scheduling_ondemand();

    DiagonalEmitter emitter(M, N, W);
    farm.add_emitter(emitter);

    if(farm.run_and_wait_end() < 0) {
        ff::error("Running farm");
        return -1;
    }

    #ifdef DEBUG
        SaveMatrixToFile(&M, N, "matrix_farm_chunk_results.txt");
    #endif
    ff::ffTime(ff::STOP_TIME);
    std::cout << "Time passed to calculate the wavefront: " << ff::ffTime(ff::GET_TIME)/1000.0 << " seconds" << std::endl;

    return 0;
}
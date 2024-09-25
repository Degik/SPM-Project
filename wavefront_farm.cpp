#include <cmath>
#include <vector>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <algorithm>

#include <ff/ff.hpp>
#include <ff/farm.hpp>
#include <ff/parallel_for.hpp>

//#define DEBUG

using vector_d = std::vector<double>;

typedef struct DiagonalTask{
    int k;
    int m;
    int row;   // m+N;
    int col_t; // (m+k)*N
    vector_d &M;
    std::atomic<int> &tasks;
};

// FillMatrix
/*!
    \name FillMatrix
    \brief Fill the matrix M with the values
    \note Fill the matrix M with the values
*/
vector_d* FillMatrix(vector_d *M, uint16_t N){
    // Fill the diagonal elements (i,j) (where i == j) with (m+1)/N
    for(int m = 0; m < N; m++){
        (*M)[m*N+m] = static_cast<double>(m+1)/N; // M[m][m] = (m+1)/N
    }
    return M;
}

/*!
    \name SaveMatrixPtrToFile
    \param M shared_ptr<vector_d> M
    \param N uint16_t N
    \param filename string filename
    \brief Save the matrix M to a file
    \note Save the matrix M to a file with the name filename
*/
void SaveMatrixToFile(vector_d *M, uint16_t N, std::string filename){
    std::ofstream file;
    file.open(filename);
    file << std::fixed << std::setprecision(6);
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            file << (*M)[i*N+j] << " ";  // M[i][j]
        }
        file << std::endl;
    }
    file.close();
}

struct DiagonalWorker: ff::ff_node_t<DiagonalTask, void>{
    void *svc(DiagonalTask *task){
        double element = 0.0;
        for(int i = 0; i < task->k; i++){
            element += (task->M)[task->row+i+task->m] * (task->M)[task->col_t+i+task->m+1];
        }
        // Calculate the cubic root
        double new_element = std::cbrt(element);
        //Update the matrix
        (task->M)[task->row+task->m+task->k] = new_element;
        (task->M)[task->col_t+task->m] = new_element;

        // Decrease the number of tasks
        task->tasks--;

        return GO_ON;
    }
};


/*!
    \name DiagonalEmitter
    \brief DiagonalEmitter
    \note DiagonalEmitter
*/
struct DiagonalEmitter: ff::ff_monode_t<int, DiagonalTask>{
    vector_d &M;
    uint16_t N;

    DiagonalEmitter(vector_d &M, uint16_t N) : M(M), N(N) {}

    DiagonalTask *svc(int*){
        // Send to the worker the index for the dot product zone
        for(int k = 1; k < N; k++){
            std::atomic<int> tasks = {N-k};
            for(int m = 0; m < N-k; m++){
                ff_send_out(new DiagonalTask{k, m, m*N, (m+k)*N, M, std::ref(tasks)});
            }
            // Wait for the tasks to finish
            while(tasks > 0){
                continue;
            }
        }
        broadcast_task(EOS);
        return EOS;
    }
};


int main(int argc, char* argv[]){
    // N, W
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << "N (Size N*N) W (Workers)" << std::endl;
        return -1;
    }

    const uint16_t N = atoi(argv[1]);
    const uint16_t W = atoi(argv[2]);

    // Create and fill the matrix
    // Process to create the matrix
    auto start = std::chrono::high_resolution_clock::now();

    vector_d M = vector_d(N*N, 0.0);
    // Fill the matrix
    FillMatrix(&M, N);
    #ifdef DEBUG
        SaveMatrixToFile(&M, N, "matrix_farm_normal.txt");
    #endif
    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> passed_time = stop - start;
    std::cout << "Matrix created and filled in: " << passed_time.count() << " seconds" << std::endl;

    // Start the timer
    ff::ffTime(ff::START_TIME);
    // Create workers
    std::vector<std::unique_ptr<ff::ff_node>> workers;
    for(int i = 0; i < W; i++){
        workers.push_back(std::make_unique<DiagonalWorker>());
    }
    // Create the farm
    ff::ff_Farm<DiagonalTask> farm(std::move(workers));
    farm.remove_collector();
    farm.wrap_around();
    farm.set_scheduling_ondemand();
    // Create the emitter
    DiagonalEmitter emitter(M, N);
    //
    farm.add_emitter(emitter);
    // Run the farm
    if(farm.run_and_wait_end() < 0){
        ff::error("Running farm");
        return -1;
    }
    
    #ifdef DEBUG
        SaveMatrixToFile(&M, N, "matrix_farm_results.txt");
    #endif
    ff::ffTime(ff::STOP_TIME);
    std::cout << "Time passed to calculate the wavefront: " << ff::ffTime(ff::GET_TIME)/1000.0 << " seconds" << std::endl;
    return 0;

}

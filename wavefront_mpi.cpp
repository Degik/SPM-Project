#include <mpi.h>
#include <cmath>
#include <chrono>
#include <vector>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

//#define DEBUG
#define TAG_TERMINATE 1
#define TAG_TASK 0

using vector_d = std::vector<double>;

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

/*!
    \name DotProduct
    \param M vector_d M
    \param k int k
    \param m int m
    \param N int N
    \brief Compute the dot product of the m-element of the k-th diagonal
    \note Compute the dot product of the m-element of the k-th diagonal
*/
double DotProductWithCbrt (const vector_d &M, int k, int m, int N){
    double sum = 0.0;
    for(int i = 0; i < k; i++){
        sum += M[m*N+i+m] * M[(m+i+1)*N+m+k];
    }
    return cbrt(sum);
}


/*!
    \name Task
    \brief Used for send the task to compute the m-element
*/
struct Task{
    int m;
};

/*!
    \name Task_Result
    \brief Used for store the result of the task
*/
struct Task_Result{
    int m;
    double value;
};


int main(int argc, char* argv[]){

    if(argc != 2){
        printf("Usage: %s <N>\n", argv[0]);
        return -1;
    }

    int N = atoi(argv[1]);
    if(N <= 1){
        printf("N must be greater than 1\n");
        return -1;
    }
    
    //Timer to measure the creation and filling of the matrix
    auto start = std::chrono::high_resolution_clock::now();
    // Create the matrix M
    vector_d M = vector_d(N*N, 0.0);
    // Fill the matrix M with the values
    FillMatrix(&M, N);
    // Save the matrix to a file if BENCHMARK is not defined
    #ifdef DEBUG
    SaveMatrixToFile(&M, N, "matrix_mpi.txt");
    #endif
    //
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds_matrix = end - start;
    printf("Matrix created and filled in: %f seconds\n", elapsed_seconds_matrix.count());

    //
    MPI_Init(&argc, &argv);                                    // Initialize the MPI environment
    //Timer to measure the wavefront algorithm
    double start_mpi_timer = MPI_Wtime();


    int rank, number_of_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);                      // Get the rank of the process
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);       // Get the number of processes

    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);              // Broadcast the size of the matrix
    MPI_Bcast(M.data(), N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);   // Broadcast the matrix

    // Iterate over the k (diagonal distance)
    for (int k = 1; k < N; k++){
        // The master process separates the work and sends it to the other processes
        if (rank == 0) {
                std::vector<Task> task_list;
                // Iterate over the m diagonal element of the k-th diagonal
                for(int m = 0; m < N-k; m++){
                    task_list.push_back({m});
                }

                // Send the work to the other processes
                MPI_Status status;
                int workers = number_of_processes - 1;
                int active_workers = 0;

                // Send task to other processes
                int next_task = 0;
                for (int i = 1; i < number_of_processes && next_task < task_list.size(); i++){
                    MPI_Send(&task_list[next_task], sizeof(Task), MPI_BYTE, i, TAG_TASK, MPI_COMM_WORLD);
                    next_task++;
                    active_workers++;
                }

                // Take the rusults from the staks
                while (active_workers > 0){
                    Task_Result task_result;
                    MPI_Recv(&task_result, sizeof(Task_Result), MPI_BYTE, MPI_ANY_TAG, MPI_ANY_SOURCE, MPI_COMM_WORLD, &status);

                    //Update the value
                    int result_idx = task_result.m * N + task_result.m + k;
                    M[result_idx] = task_result.value;

                    //Take the worker rank
                    int woker_rank = status.MPI_SOURCE;
                    // Check if there is another task to assign
                    if (next_task < task_list.size()){
                        MPI_Send(&task_list[next_task], sizeof(Task), MPI_BYTE, woker_rank, TAG_TASK, MPI_COMM_WORLD);
                        next_task++;
                    }else{
                        MPI_Send(NULL, 0, MPI_BYTE, woker_rank, TAG_TERMINATE, MPI_COMM_WORLD);
                        active_workers--;
                    }
                }

                //Reults vector
                std::vector<double> k_diagonal(N-k, 0.0);
                //Take all the new results
                for (int i = 0; i < N - k; i++){
                    k_diagonal[i] = M[i * N + i + k];
                }

                #ifdef DEBUG
                std::cout << "Master Worker - " << "Sending the new k-" << k << "-diagonal" << std::endl;
                std::cout << "{";
                std::cout << std::fixed << std::showpoint;
                std::cout << std::setprecision(6);
                for (int i = 0; i < k_diagonal.size(); i++){
                    std::cout << " ";
                    std::cout << k_diagonal[i]; 
                }
                std::cout << " }" << std::endl;
                #endif

                MPI_Bcast(k_diagonal.data(), N-k, MPI_DOUBLE, 0, MPI_COMM_WORLD);   // Sends to all the computed k_diagonal

        } else {
            while (true){
                // The other processes receive the work from the master process
                Task task;
                MPI_Status status;
                MPI_Recv(&task, sizeof(Task), MPI_BYTE, MPI_ANY_TAG, MPI_ANY_SOURCE, MPI_COMM_WORLD, &status);

                if (status.MPI_TAG == TAG_TERMINATE){ break; }          // In case we interrupt the worker

                // Compute the dot product and the cubic root on the sum
                double sum = DotProductWithCbrt(M, k, task.m, N);       

                #ifdef DEBUG
                std::cout << "I'm the worker " << rank << " and i've calculated the sum like: " << sum << std::endl;
                #endif

                Task_Result result = {task.m, sum};                     // Data to send

                MPI_Send(&result, sizeof(Task_Result), MPI_BYTE, 0, TAG_TASK, MPI_COMM_WORLD);
            }
            std::vector<double> k_diagonal(N-k, 0.0);
            MPI_Bcast(k_diagonal.data(), N*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);    // Receive the new k_diagonal
            // Update the matrix with the current k_diagonal
            for (int i = 0; i < N - k; ++i) {
                M[i * N + i + k] = k_diagonal[i];
            }
        }
        // Wait for all workers
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (rank == 0){
        #ifdef DEBUG
            SaveMatrixToFile(&M, N, "matrix_mpi_result.txt");
        #endif
        double end_mpi_timer = MPI_Wtime();
        double passed_time = end_mpi_timer - start_mpi_timer;
        std::cout << "Time to compute the matrix: " << passed_time << std::endl;
    }
    MPI_Finalize();                                                     // Finalize the MPI environment
    return 0;
}

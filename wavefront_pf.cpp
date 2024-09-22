#include <vector>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <iostream>


#include <ff/ff.hpp>
#include <ff/farm.hpp>
#include <ff/parallel_for.hpp>

//#define DEBUG

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


int main(int argc, char* argv[]){
    // N, W
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << "N (Size N*N) W (Workers)" << std::endl;
        return -1;
    }

    const uint16_t N = atoi(argv[1]);
    const uint16_t W = atoi(argv[2]);
  
    // Process to create the matrix
    auto start = std::chrono::high_resolution_clock::now();

    vector_d M = vector_d(N*N, 0.0);
    // Fill the matrix
    FillMatrix(&M, N);
    #ifdef DEBUG-
        SaveMatrixToFile(&M, N, "matrix_pf_normal.txt");
    #endif
    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> passed_time = stop - start;
    std::cout << "Matrix created and filled in: " << passed_time.count() << " seconds" << std::endl;

    // Start the timer
    ff::ffTime(ff::START_TIME);

    ff::ParallelFor pf(W);
    for (int k = 1; k < N; k++){
        pf.parallel_for(0, N-k, [&](const long m){
            double element = 0.0;
            for(int i = 0; i < k; i++){
                element += M[m*N+i+m] * M[(m+i+1)*N+m+k];
            }
            M[m*N+m+k] = cbrt(element);
        });
    }

    #ifdef DEBUG
        SaveMatrixToFile(&M, N, "matrix_pf_results.txt");
    #endif
    ff::ffTime(ff::STOP_TIME);
    std::cout << "Time passed to calculate the wavefront: " << ff::ffTime(ff::GET_TIME)/1000.0 << " seconds" << std::endl;
    return 0;

}

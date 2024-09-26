#include <cmath>
#include <vector>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <iostream>

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
    // N
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << "N (Size N*N)" << std::endl;
        return -1;
    }

    const uint16_t N = atoi(argv[1]);

    // Process to create the matrix
    auto start = std::chrono::high_resolution_clock::now();

    vector_d M = vector_d(N*N, 0.0);
    // Fill the matrix
    FillMatrix(&M, N);
    #ifdef DEBUG
        SaveMatrixToFile(&M, N, "matrix_seq_cache_normal.txt");
    #endif
    auto stop = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> passed_time = stop - start;
    std::cout << "Matrix created and filled in: " << passed_time.count() << " seconds" << std::endl;

    //Wavefront sequential
    auto start_wavefront = std::chrono::high_resolution_clock::now();

    for (int k = 1; k < N; k++){
        for(int m = 0; m < N-k; m++){
            int row = m*N;
            int col_t = (m+k)*N;
            double element = 0.0;
            for (int i = 0; i < k; i++){
                element += M[row+i+m] * M[col_t+i+m+1]; //M[m][i+m] * M[m+k][i+m+1]
            }
            double new_element = std::cbrt(element);
            M[row+m+k] = new_element;
            M[col_t+m] = new_element; // Update the element for the transpose matrix
        }
    }

    #ifdef DEBUG
        SaveMatrixToFile(&M, N, "matrix_seq_cache_results.txt");
    #endif
    auto stop_wavefront = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = stop_wavefront - start_wavefront;
    std::cout << "Time passed to calculate the wavefront: " << elapsed_time.count() << " seconds" << std::endl;
    return 0;
}

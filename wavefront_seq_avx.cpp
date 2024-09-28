#include <cmath>
#include <vector>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <immintrin.h>

//#define DEBUG

using vector_d = std::vector<double>;

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

void CreateMatrix(vector_d &M, uint16_t N){
    // Crete the value
    __m256d zero = _mm256_set1_pd(0.0);
    //
    int i = 0;
    int total = N * N;
    // https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#ig_expand=10,0,5812,6534&text=_mm256_storeu_pd
    for(; i < total - 4; i += 4){
        _mm256_storeu_pd(&M[i], zero);
    }
}

void FillMatrix(vector_d *M, uint16_t N){

    for(int m = 0; m < N; m++){
        (*M)[m*N+m] = static_cast<double>(m+1)/N; // M[m][m] = (m+1)/N
    }
}

void ComputeWavefrontAVX(vector_d *M, uint16_t N){
    for (int k = 1; k < N; k++){
        for(int m = 0; m < N-k; m++){
            int row = m*N;
            int col_t = (m+k)*N;
            double element = 0.0;

            int i = 0;
            // Initialize sum_vec to zero
            __m256d sum_vec = _mm256_setzero_pd();
            // Process 4 elements at a time
            for (; i <= k - 4; i += 4){
                __m256d vec1 = _mm256_loadu_pd(&(*M)[row + i + m]);         // Load 4 elements from the row
                __m256d vec2 = _mm256_loadu_pd(&(*M)[col_t + i + m + 1]);   // Load 4 elements from the column_transpose
                __m256d prod = _mm256_mul_pd(vec1, vec2);                   // Multiply the two vectors -> prod = [a*b, c*d, e*f, g*h]
                sum_vec = _mm256_add_pd(sum_vec, prod);                     // Move the results to the sum_vector
            }
            // Extract the sum from the vector
            __m128d sum_high = _mm256_extractf128_pd(sum_vec, 1);           // Extract the last 128 bits - Latency 3 cycles - Throughput 1 cycle
            __m128d sum_low = _mm256_castpd256_pd128(sum_vec);              // Take the first 128 bits without cost - Latency 1 cycle
            __m128d sum_128 = _mm_add_pd(sum_low, sum_high);                // Sum the two 128 bits vectos -> sum_low = [a, b] sum_high = [c, d] -> sum_128 = [a+c, b+d]
            sum_128 = _mm_hadd_pd(sum_128, sum_128);                        // Horizontal add -> sum_128_before = [a+c, b+d] -> sum_128 = [a+c+b+d, a+c+b+d]
            double sum;
            _mm_store_sd(&sum, sum_128);                                    // Store the result in a double
            element += sum;

            // Process the elements out of the block of 4
            for (; i < k; i++){
                element += (*M)[row + i + m] * (*M)[col_t + i + m + 1];
            }

            double new_element = std::cbrt(element);
            (*M)[row + m + k] = new_element;        // Update the element
            (*M)[col_t + m] = new_element;          // Update the element for the transpose matrix
        }
    }
}


int main(int argc, char *argv[]){
    // N
    if (argc != 2){
        std::cout << "Usage: " << argv[0] << "N (Size N*N)" << std::endl;
        return -1;
    }

    uint16_t N = atoi(argv[1]);
    
    auto start_timer = std::chrono::high_resolution_clock::now();

    vector_d M(N*N);

    CreateMatrix(M, N);
    FillMatrix(&M, N);
    
    #ifdef DEBUG
        SaveMatrixToFile(&M, N, "wavefront_seq_avx_normal.txt");
    #endif

    auto stop_timer = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = stop_timer - start_timer;
    std::cout << "Matrix created and filled in: " << time.count() << " seconds" << std::endl;

    start_timer = std::chrono::high_resolution_clock::now();
    //
    ComputeWavefrontAVX(&M, N);
    //
    #ifdef DEBUG
        SaveMatrixToFile(&M, N, "wavefront_seq_avx_results.txt");
    #endif
    //
    stop_timer = std::chrono::high_resolution_clock::now();
    time = stop_timer - start_timer;
    std::cout << "Time passed to calculate the wavefront: " << time.count() << " seconds" << std::endl;
    
    return 0;
}
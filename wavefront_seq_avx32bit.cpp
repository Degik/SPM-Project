#include <cmath>
#include <vector>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <immintrin.h>

//#define DEBUG

using vector_d = std::vector<float>;

/*!
    \name SaveMatrixPtrToFile
    \param M vector_d M
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
    \name CreateMatrix
    \param M vector_d M
    \param N uint16_t N
    \brief Create the matrix M with the size N*N
    \note Create the matrix M with the size N*N
*/
void CreateMatrix(vector_d &M, uint16_t N){
    // Crete the value
    __m256 zero = _mm256_set1_ps(0.0f);
    //
    int i = 0;
    int total = N * N;
    // https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#ig_expand=10,0,5812,6534&text=_mm256_storeu_pd
    for(; i < total - 8; i += 8){
        _mm256_storeu_ps(&M[i], zero);
    }
}

/*!
    \name FillMatrix
    \param M vector_d M
    \param N uint16_t N
    \brief Fill the matrix M with the values
    \note Fill the matrix M with the values
*/
void FillMatrix(vector_d *M, uint16_t N){

    for(int m = 0; m < N; m++){
        (*M)[m*N+m] = static_cast<float>(m+1)/N; // M[m][m] = (m+1)/N
    }
}

/*!
    \name ComputeWavefrontAVX
    \param M vector_d M
    \param N uint16_t N
    \brief Compute the wavefront using AVX
    \note Compute the wavefront using AVX - AVX 256 bits - 8 elements at a time
*/
void ComputeWavefrontAVX(vector_d *M, uint16_t N){
    for (int k = 1; k < N; k++){
        for(int m = 0; m < N-k; m++){
            int row = m*N;
            int col_t = (m+k)*N;
            float element = 0.0f;

            int i = 0;
            // Initialize sum_vec to zero
            __m256 sum_vec = _mm256_setzero_ps();
            // Process 8 elements at a time
            for (; i <= k - 8; i += 8){
                __m256 vec1 = _mm256_loadu_ps(&(*M)[row + i + m]);         // Load 8 elements from the row
                __m256 vec2 = _mm256_loadu_ps(&(*M)[col_t + i + m + 1]);   // Load 8 elements from the column_transpose
                __m256 prod = _mm256_mul_ps(vec1, vec2);                   // Multiply the two vectors -> prod = [a*b, c*d, e*f, g*h, i*j, k*l, m*n, o*p]
                sum_vec = _mm256_add_ps(sum_vec, prod);                    // Move the results to the sum_vector
            }
            // Extract the sum from the vector
            __m128 sum_high = _mm256_extractf128_ps(sum_vec, 1);           // Extract the last 128 bits - Latency 3 cycles - Throughput 1 cycle
            __m128 sum_low = _mm256_castps256_ps128(sum_vec);              // Take the first 128 bits without cost - Latency 1 cycle
            __m128 sum_128 = _mm_add_ps(sum_low, sum_high);                // Sum the two 128 bits vectos -> sum_low = [a, b, c, d] sum_high = [e, f, g, h] -> sum_128 = [a+e, b+f, c+g, d+h]
            sum_128 = _mm_hadd_ps(sum_128, sum_128);                       // Horizontal add -> sum_128_before = [a+e, b+f, c+g, d+h] -> sum_128 = [a+e+b+f, c+g+d+h]
            // sum the elements of the vector
            sum_128 = _mm_hadd_ps(sum_128, sum_128);                       // Horizontal add -> sum_128_before = [a+e+b+f, c+g+d+h] -> sum_128 = [a+e+b+f+c+g+d+h]
            float sum;
            _mm_store_ss(&sum, sum_128);                                   // Store the result in a double
            element += sum;

            // Process the elements out of the block of 8
            for (; i < k; i++){
                element += (*M)[row + i + m] * (*M)[col_t + i + m + 1];
            }

            float new_element = std::cbrtf(element);
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
        SaveMatrixToFile(&M, N, "wavefront_seq_avx32bit_normal.txt");
    #endif

    auto stop_timer = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = stop_timer - start_timer;
    std::cout << "Matrix created and filled in: " << time.count() << " seconds" << std::endl;

    start_timer = std::chrono::high_resolution_clock::now();
    //
    ComputeWavefrontAVX(&M, N);
    //
    #ifdef DEBUG
        SaveMatrixToFile(&M, N, "wavefront_seq_avx32bit_results.txt");
    #endif
    //
    stop_timer = std::chrono::high_resolution_clock::now();
    time = stop_timer - start_timer;
    std::cout << "Time passed to calculate the wavefront: " << time.count() << " seconds" << std::endl;
    
    return 0;
}
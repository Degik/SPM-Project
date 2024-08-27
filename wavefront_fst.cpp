#include <sys/time.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <cstdlib>
#include <numeric>
#include <vector>
#include <thread>
#include <memory>
#include <mutex>
#include <list>

#include <ff/ff.hpp>
#include <ff/parallel_for.hpp>
#include <ff/pipeline.hpp>

using namespace std;
using namespace ff;

using matrix_d = vector<vector<double>>;
using vector_d = vector<double>;

typedef struct Resources{
    uint16_t Z; // Z = Workers for the M-Diagonal-Stage
    uint16_t D; // D = Workers for the Dot-Product-Stage
} Resources;

Resources CalculateResources(uint16_t W, uint16_t K, uint16_t N){
    Resources resources = {0, 0};
    if (K < W){
        resources.D = W / 2;
        resources.Z = W / 2;
    } else {
        resources.D = (float)(W / 2) + ((float)(W / 2) * ((float)(K - 1) / (float)K));
        resources.Z = W - resources.D;
    }
    return resources;
}



// FillMatrix
vector<vector<double>>& FillMatrix(vector<vector<double>>& M, uint16_t N, uint16_t W){
    ParallelFor pf(W); // Create the parallel for object with W workers
    // Fill the diagonal elements (i,j) (where i == j) with (m+1)/N
    pf.parallel_for(0, N, [&](const long m){
        M[m][m] = static_cast<double>(m+1)/N;
    });
    return M;
}

// Print matrix M
void PrintMatrix(vector<vector<double>>& M, uint16_t N){
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            cout << M[i][j] << " ";
        }
        cout << endl;
    }
}

void SaveMatrixToFile(vector<vector<double>>& M, uint16_t N){
    ofstream file;
    file.open("matrix.txt");
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            file << M[i][j] << " ";
        }
        file << endl;
    }
    file.close();
}

struct CreateMatrix: ff_node_t<int, vector<vector<double>>> {
    uint16_t N, W;
    CreateMatrix(uint16_t size, uint16_t workers): N(size), W(workers) {} // Constructor
    
    vector<vector<double>> *svc(int *task){ // Service
        vector<vector<double>> M (N, vector<double>(N, 0.0)); // Create the matrix
        M = FillMatrix(M, N, W);     // Fill the matrix
        //PrintMatrix(M, N);           // Print the matrix
        //SaveMatrixToFile(M, N);      // Save the matrix to a file
        return new vector<vector<double>>(M);
    }
    
};


/*!
   \brief Calculate the partial dot product of v1,v2
*/
double PartialDotProduct(const vector_d &v1, const vector_d &v2, uint16_t size){
    double partial_sum = 0.0;
    
    for(uint16_t idx = 0; idx < size; idx++){
        partial_sum += v1[idx] * v2[idx];
    }
    return partial_sum;
}


/*!
    \name SplitVector
    \note This function takes one vector and returns a list with the sub-vectors
       \n The last element can have a different size

*/
list<vector_d> SplitVector(vector_d v, uint16_t D, uint16_t K){
    size_t base_size = D / K;
    size_t remainder = D % K;
    size_t start = 0;

    // All sub vectors are stored in this list
    list<vector_d> splitted_vectors;

    for (size_t i = 0; i < K; ++i) {
        size_t current_size = base_size + (i < remainder ? 1 : 0);
        splitted_vectors.push_back(vector_d(v.begin() + start, v.begin() + start + current_size));
        start += current_size;
    }

    return splitted_vectors;
}

// FastFlow Functions
//
// DotProduct_Worker
// Take the tuple (v1, v2, start, ending)
/*!
    \brief Calculate the dot product fof the sub vector
    \return double value (dot product)
*/
struct DotProduct_Worker: ff_node_t<tuple<vector_d, vector_d, int>, double>{
    double* svc(tuple<vector_d, vector_d, int> *task){
        const vector_d& v1 = std::get<0>(*task);
        const vector_d& v2 = std::get<1>(*task);
        const int size = std::get<2>(*task);

        double partial_result = PartialDotProduct(v1, v2, size);
        ff_send_out(new double(partial_result)); // send the result

        return GO_ON;
    }

};

/*!
    \name DotProductStage
    \brief Take the v1 and v2 and split them for calculate the dot product
    \note From the v1 and v2 vectors split each one vector into v1-sub-vectors and v2-sub-vectors for N-vectors
       \n It depends from the D (workers)  (N) number of split = K / D or N = K (Optimal case)
       \n Call DotProduct_Worker for calculate each sub dot product on this sub vectors
       \n Take all the results and calculate cbrt(sum) and then return the value
*/
struct DotProduct_Stage: ff_node_t<tuple<vector_d, vector_d>, double> {
    uint16_t N, K, W; // Size, K, Number of workers
    uint16_t D;       // Number of workers for M_Diagonal

    DotProduct_Stage(uint16_t N, uint16_t K, uint16_t W, uint16_t Z): N(N), K(K), W(W), D(D) {}

    // Take the vectors (v1, v2) created on the previous stage
    double* svc(tuple<vector_d, vector_d> *vectors){
        const vector_d v1 = std::get<0>(*vectors);
        const vector_d v2 = std::get<1>(*vectors);
        vector<unique_ptr<ff_node>> workers;
        //
        list<vector_d> sub_v1_list(D);
        list<vector_d> sub_v2_list(D);
        //
        for(uint16_t w = 0; w < D; w++){
                workers.push_back(make_unique<DotProduct_Worker>());
        }

        ff_Farm<double> farm(move(workers));

        // if K is less or equal with D we can have 1:(1,1) 
        // [one workers:for one elememt] for calculate the dot product
        sub_v1_list = SplitVector(v1, D, K);
        sub_v2_list = SplitVector(v2, D, K);

        farm.wrap_around();
        farm.set_scheduling_ondemand();
        farm.run_then_freeze();

        for(uint16_t w = 0; w < D; w++){
            //Take the begin of the lists each time
            const auto begin_v1 = sub_v1_list.begin();
            const auto begin_v2 = sub_v2_list.begin();
            //
            advance(begin_v1, w); // Move the begin iterator to the w-th element
            advance(begin_v2, w); // Move the begin iterator to the w-th element
            //
            const uint16_t size = begin_v1->size(); // v1 and v2 have the same size
            farm.ff_send_out(new tuple<vector_d, vector_d, int>(*begin_v1, *begin_v2, size));
        }
        
        // Take the results from the workers
        double sum = 0.0;
        

        // Update the matrix
    }
};

struct M_Diagonal_Stage: ff_node_t<matrix_d, matrix_d> { // Take matrix M and return M'
    uint16_t N, K, W;
    uint16_t Z, D;

    M_Diagonal_Stage(uint16_t N, uint16_t K, uint16_t W, uint16_t Z, uint16_t D): N(N), K(K), W(W), Z(Z), D(D) {}

    vector<vector<double>> *svc(vector<vector<double>> *M){
        
        vector<unique_ptr<ff_node>> workers;
        // Cycle with m with m = [0, n-k[
        for (uint16_t w = 0; w < Z; w++){
            workers.push_back(make_unique<DotProduct_Stage>(N, K, W, D));
        }

        ff_Farm<double> farm(move(workers));



    }
};

int main(int argc, char* arg[]){
    // N, K, W
    if (argc != 4) {
        cout << "Usage: " << arg[0] << "N (Size N*N) K (K value) W (Workers)" << endl;
        return -1;
    }

    const uint16_t N = atoi(arg[1]);
    const uint16_t K = atoi(arg[2]);
    const uint16_t W = atoi(arg[3]);

    cout << "N: " << N << " K: " << K << " W: " << W << endl;

    // Calculate how to split the workers W
    Resources resources = CalculateResources(W, K, N);
    
    cout << "DiagonalStage-Workers Z: " << resources.Z << " DotProductStage-Workers D: " << resources.D << endl;

    ffTime(START_TIME);

    // Stage 1
    // Create the matrix M
    // Fill the matrix M with the values (m+1)/N
    CreateMatrix s1(N, W);
    // Stage 2
    // Create the dot product of the matrix M
    
    // Stage 3
    // Merge the results of the dot product

    ff_Pipe<> pipe(s1);
    /*if (pipe.run_and_wait_end() < 0) {
        error("Running pipe\n");
        return -1;
    }
    cout << "End" << endl;*/
    int task = 0;
    s1.svc(&task);
    ffTime(STOP_TIME);
    cout << "Time: " << ffTime(GET_TIME)/1000.0 << endl;

    // CREATE MATRIX M (Stage 1)
    // STAGE 1 create all subelements of the matrix for the diagoanl
    // STAGE 2 create the dot product of the matrix

    // SPLIT THE PROBLEME WITH THE PIPE()
    // CREATE THE FARM FOR GENERETING THE DOT PRODUCT
}

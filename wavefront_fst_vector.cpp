#include <sys/time.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <iomanip>
#include <cstdlib>
#include <numeric>
#include <vector>
#include <thread>
#include <memory>
#include <mutex>
#include <list>

#include <ff/ff.hpp>
#include <ff/farm.hpp>
#include <ff/parallel_for.hpp>
#include <ff/pipeline.hpp>

using namespace std;
using namespace ff;

using vector_d = vector<double>;
using tuple_for_worker = tuple<double*, double*, int, int>;
using tuple_dot_product = tuple<shared_ptr<vector_d>, int, uint16_t, uint16_t, uint16_t, uint16_t>;

typedef struct Resources{
    uint16_t Z; // Z = Workers for the M-Diagonal-Stage
    uint16_t D; // D = Workers for the Dot-Product-Stage
} Resources;


// CalculateResources
/*!
    \name CalculateResources
    \brief Calculate the resources for the workers
    \note Calculate the resources for the workers based on the values of W, K, N and return the resources for the workers (Z, D)
*/
Resources CalculateResources(uint16_t W, uint16_t K, uint16_t N) {
    Resources resources;
    resources.Z = max(1, static_cast<int>(W * static_cast<float>(N - K) / N));
    resources.D = W - resources.Z;
    return resources;
}



// FillMatrix
/*!
    \name FillMatrix
    \brief Fill the matrix M with the values
    \note Fill the matrix M with the values
*/
shared_ptr<vector_d> FillMatrix(shared_ptr<vector_d> M, uint16_t N, uint16_t W){
    ParallelFor pf(W); // Create the parallel for object with W workers
    // Fill the diagonal elements (i,j) (where i == j) with (m+1)/N
    pf.parallel_for(0, N, [&](const long m){
        (*M)[m*N+m] = static_cast<double>(m+1)/N; // M[m][m] = (m+1)/N
    });
    return M;
}

// Print matrix M
/*!
    \name PrintMatrix
    \brief Print the matrix M
    \note Print the matrix M
*/
void PrintMatrix(shared_ptr<vector_d> M, uint16_t N){
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            printf("%.6f ", (*M)[i*N+j]); // M[i][j]
        }
        printf("\n");
    }
}

/*!
    \name SaveMatrixPtrToFile
    \param M shared_ptr<vector_d> M
    \param N uint16_t N
    \param filename string filename
    \brief Save the matrix M to a file
    \note Save the matrix M to a file with the name filename
*/
void SaveMatrixToFile(shared_ptr<vector_d> M, uint16_t N, string filename){
    ofstream file;
    file.open(filename);
    vector_d& matrix = *M;
    file << fixed << setprecision(6);
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            file << (*M)[i*N+j] << " ";  // M[i][j]
        }
        file << endl;
    }
    file.close();
}

// FastFlow Functions
//
/*!
    \name CreateMatrix
    \brief Create the matrix M
    \note Create the matrix M with the size N and fill it with the values
*/
struct CreateMatrix: ff_node_t<int, shared_ptr<vector_d>> {
    uint16_t N, W;
    CreateMatrix(uint16_t size, uint16_t workers): N(size), W(workers) {} // Constructor
    
    shared_ptr<vector_d> *svc(int *task) override { // Service
        static auto M = make_shared<vector_d>(N*N, 0.0);
        //vector_d M (N, vector_d(N, 0.0)); // Create the matrix
        M = FillMatrix(M, N, W);                   // Fill the matrix
        vector_d& matrix = *M;
        //PrintMatrix(M, N);                         // Print the matrix
        //SaveMatrixToFile(M, N, "matrix.txt");      // Save the matrix to a file
        ff_send_out(&M);                           // Send the matrix to the next stage
        return EOS;
    }
    
};

/*!
    \name Sink
    \brief Calculate the cbrt(sum) and update the matrix
*/
struct Sink: ff_node_t<double, double>{
    int i, j;    // Indices of the matrix to update
    uint16_t N; // Size of the matrix
    size_t limit;
    int count = 0; // Ex v1.size() = 2, count stop at 2
    double sum = 0.0;
    shared_ptr<vector_d> M; // Reference to the matrix
    Sink(shared_ptr<vector_d> M, uint16_t N, int i, int j, size_t limit) : M(M), N(N), i(i), j(j), limit(limit) {
        // cout << "Indirizzo di M nel Sink: " << M.get() << endl;
        // cout << "Limit: " << limit << endl;
    }

    double* svc(double* task){
        if (count != limit){
            sum += *task;
            delete task;
            count++;
            return GO_ON;
        } else {
            return EOS;
        }
    }

    void svc_end(){ // vector_d, int(i), int(j)
        // Update the matrix
        const double element = cbrt(sum);
        (*M)[i*N+j] = element;   // M[i][j]
        //(*M)[j*N+i] = element;   // M[j][i]
        // cout << "Ho aggiornato la matrice con il valore: " << element << " nella posizione M[" << i << "][" << j << "]" << endl;
    }
};

/*!
    \name DotProduct_Emitter
    \brief Split the vectors v1 and v2 for the workers
    \note Split the vectors v1 and v2 for the workers
*/// DotProduct_Emitter
struct DotProduct_Emitter: ff_node_t<tuple_for_worker> {
    const int m, N;
    uint16_t K, W, D;
    shared_ptr<vector_d> M;

    DotProduct_Emitter(uint16_t K, uint16_t W, uint16_t D, shared_ptr<vector_d> M, int m, int N)
        : K(K), W(W), D(D), M(M), m(m), N(N) {}

    tuple_for_worker* svc(tuple_for_worker*) {
        const int limit = std::min(static_cast<int>(K), N - m);
        // Calculate the number of workers to use
        // D = min(D, limit)
        uint16_t actual_workers = std::min(D, static_cast<uint16_t>(limit));
        
        uint16_t base_size = limit / actual_workers;   // Base size
        uint16_t remainder = limit % actual_workers;   // Remainder

        uint16_t start = 0;
        
        for (uint16_t w = 0; w < actual_workers; ++w) {
            uint16_t current_size = base_size + (w < remainder ? 1 : 0);
            uint16_t end = start + current_size;

            double* v1 = &(*M)[m * N + (m + start)];
            double* v2 = &(*M)[(m + start + 1) * N + (m + K)];

            ff_send_out(new tuple_for_worker(v1, v2, end - start, N));

            start = end;
        }

        return EOS;
    }
};

/*!
    \name Diagonal_Emitter
    \brief Take the matrix M and split it for the workers
    \note Take the matrix M and split it for the workers
*/
struct Diagonal_Emitter: ff_node_t<tuple_dot_product>{
    shared_ptr<vector_d> M;
    uint16_t N, K, W, D;

    Diagonal_Emitter(shared_ptr<vector_d> M, uint16_t N, uint16_t K, uint16_t W, uint16_t D)
        : M(M), N(N), K(K), W(W), D(D) {}

    tuple_dot_product* svc(tuple_dot_product*){
        // cout << "m = [0, " << N - K << "[" << endl;
        //vector_d& matrix = *M;
        for (int m = 0; m < N - K; ++m){

            // cout << "Sending the tuple to the farm" << endl;
            ff_send_out(new tuple_dot_product(M, m, K, N, W, D));
            // cout << "Tuple sent" << endl;
        }
        return EOS;
    }
};

// DotProduct_Worker
// Take the tuple (v1, v2, start, ending)
/*!
    \brief Calculate the dot product fof the sub vector
    \return double value (dot product)
*/
struct DotProduct_Worker: ff_node_t<tuple_for_worker, double>{
    double* svc(tuple_for_worker *task){
        double* v1 = std::get<0>(*task);
        double* v2 = std::get<1>(*task);
        const int length = std::get<2>(*task); // length of the sub-vector
        const int N = std::get<3>(*task); // N
        double partial_result = 0.0; // Partial result

        for (int i = 0, j = 0; i < length; i++) {
            partial_result += v1[i] * v2[j];
            j += N;
        }

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
struct DotProduct_Stage: ff_node_t<tuple_dot_product, double> {
    // Take the vectors (v1, v2) created on the previous stage
    double* svc(tuple_dot_product *task){
        shared_ptr<vector_d> M = std::get<0>(*task);
        const int m = std::get<1>(*task); // m
        const uint16_t K = std::get<2>(*task); // K
        const uint16_t N = std::get<3>(*task); // N
        const uint16_t W = std::get<4>(*task); // Number of workers
        const uint16_t D = std::get<5>(*task); // Number of workers for the DotProduct


        vector<unique_ptr<ff_node>> workers;

        for(uint16_t w = 0; w < D; w++){
            workers.push_back(make_unique<DotProduct_Worker>());
        }
        // cout << "v1_size or Limit: " << v1_size << endl;
        Sink sink(M, N, m, m+K, K);                             // Create the sink
        int int_N = static_cast<int>(N);
        DotProduct_Emitter emitter(K, W, D, M, m, int_N); // Create the emitter
        ff_Farm<tuple<vector_d, vector_d, int>> farm(move(workers));
        // Add the emitter and the sink
        farm.add_emitter(emitter);
        farm.add_collector(sink);
        // cout << "Farm created (DotProduct)" << endl;

        //farm.wrap_around();
        farm.set_scheduling_ondemand();
        if (farm.run_and_wait_end() < 0){
            error("Running farm (DotProduct)\n");
            return EOS;
        }
        return GO_ON;
    }
};

struct M_Diagonal_Stage: ff_node_t<shared_ptr<vector_d>, shared_ptr<vector_d>> { // Take matrix M and return M'

private:
    uint16_t N, K, W;
    uint16_t Z, D;
public:
    M_Diagonal_Stage(uint16_t N, uint16_t K, uint16_t W, uint16_t Z, uint16_t D): N(N), K(K), W(W), Z(Z), D(D) {}

    shared_ptr<vector_d> *svc(shared_ptr<vector_d> *M) {
        
        vector<unique_ptr<ff_node>> workers;
        // Cycle with m with m = [0, n-k[
        for (uint16_t w = 0; w < Z; w++){
            // cout << "Worker - M: " << w << endl;
            workers.push_back(make_unique<DotProduct_Stage>());
        }
        Diagonal_Emitter emitter(*M, N, K, W, D);
        ff_Farm<tuple_dot_product, shared_ptr<vector_d>> farm(move(workers), emitter);
        // cout << "Farm created (M-Diagonal)" << endl;
        farm.wrap_around();
        //farm.remove_collector();
        farm.set_scheduling_ondemand();
        
        if (farm.run_and_wait_end() < 0){
            error("Running farm (M-Diagonal)\n");
            return nullptr;
        }
        ff_send_out(M);
        return EOS;
    }
};

struct SaveMatrix_Stage: ff_node_t<shared_ptr<vector_d>, void>{
    void *svc(shared_ptr<vector_d> *task) {
        //vector_d& matrix = **task;
        uint16_t size = (*task)->size();
        //SaveMatrixToFile(*task, sqrt(size), "matrix_prime.txt");
        return EOS;
    }
};

int main(int argc, char* arg[]){
    // N, W
    if (argc != 3) {
        cout << "Usage: " << arg[0] << "N (Size N*N) W (Workers)" << endl;
        return -1;
    }

    const uint16_t N = atoi(arg[1]);
    const uint16_t W = atoi(arg[2]);

    if(N <= 1 || W <= 1){
        cout << "N and W must be greater than 1" << endl;
        return -1;
    }

    cout << "N: " << N << " W: " << W << endl;


    ffTime(START_TIME);

    // Stage 1
    CreateMatrix s1(N, W);
    ff_Pipe<> pipe(s1);
    //pipe.add_stage(s1);
    //Stage 2 -> N - K
    vector<unique_ptr<ff_node>> workers;
    for (uint16_t k = 1; k < N; k++){
        Resources resources = CalculateResources(W, k, N);
        // cout << "K = " << k << endl;
        // cout << "DiagonalStage-Workers Z: " << resources.Z << " DotProductStage-Workers D: " << resources.D << endl;
        // cout << endl;
        ff_node* worker = new M_Diagonal_Stage(N, k, W, resources.Z, resources.D);
        pipe.add_stage(*worker);
    }
    // Last Stage
    SaveMatrix_Stage s3;
    pipe.add_stage(s3);
    if (pipe.run_and_wait_end() < 0){
        error("Running pipe\n");
        return -1;
    }
    ffTime(STOP_TIME);
    cout << "Time: " << ffTime(GET_TIME)/1000.0 << endl;
}

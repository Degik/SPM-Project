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
using matrix_d = vector<vector<double>>;
using tuple_dot_product = tuple<vector_d, vector_d, shared_ptr<matrix_d>, uint16_t, uint16_t, uint16_t, int, int>;

mutex mtx; // Mutex for the print

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
/*!
    \name FillMatrix
    \brief Fill the matrix M with the values
    \note Fill the matrix M with the values
*/
shared_ptr<matrix_d> FillMatrix(shared_ptr<matrix_d> M, uint16_t N, uint16_t W){
    ParallelFor pf(W); // Create the parallel for object with W workers
    // Fill the diagonal elements (i,j) (where i == j) with (m+1)/N
    pf.parallel_for(0, N, [&](const long m){
        (*M)[m][m] = static_cast<double>(m+1)/N;
    });
    return M;
}

// Print matrix M
/*!
    \name PrintMatrix
    \brief Print the matrix M
    \note Print the matrix M
*/
void PrintMatrix(shared_ptr<matrix_d> M, uint16_t N){
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            printf("%.6f ", (*M)[i][j]);
        }
        printf("\n");
    }
}

/*!
    \name SaveMatrixPtrToFile
    \param M shared_ptr<matrix_d> M
    \param N uint16_t N
    \param filename string filename
    \brief Save the matrix M to a file
    \note Save the matrix M to a file with the name filename
*/
void SaveMatrixToFile(shared_ptr<matrix_d> M, uint16_t N, string filename){
    ofstream file;
    file.open(filename);
    //matrix_d& matrix = *M;
    file << fixed << setprecision(6);
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            file << (*M)[i][j] << " ";
        }
        file << endl;
    }
    file.close();
}

/*!
    \name PartialDotProduct
    \brief Calculate the partial dot product of the vectors
    \note Calculate the partial dot product of the vectors v1 and v2 with the size
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
    uint16_t base_size = (K < D) ? D / K : K / D;
    uint16_t remainder = (K < D) ? D % K : K % D;
    uint16_t start = 0;

    // All sub vectors are stored in this list
    list<vector_d> splitted_vectors;

    // If K is less or equal with D we can have 1:(1,1)
    if (K <= 1){ splitted_vectors.push_back(v); return splitted_vectors; }

    for (size_t i = 0; i < K; ++i) {
        size_t current_size = base_size + (i < remainder ? 1 : 0);
        splitted_vectors.push_back(vector_d(v.begin() + start, v.begin() + start + current_size));
        start += current_size;
    }

    return splitted_vectors;
}

// FastFlow Functions
//
/*!
    \name CreateMatrix
    \brief Create the matrix M
    \note Create the matrix M with the size N and fill it with the values
*/
struct CreateMatrix: ff_node_t<int, shared_ptr<matrix_d>> {
    uint16_t N, W;
    CreateMatrix(uint16_t size, uint16_t workers): N(size), W(workers) {} // Constructor
    
    shared_ptr<matrix_d> *svc(int *task) override { // Service
        static auto M = make_shared<matrix_d>(N, vector_d(N, 0.0));
        //matrix_d M (N, vector_d(N, 0.0)); // Create the matrix
        M = FillMatrix(M, N, W);                   // Fill the matrix
        matrix_d& matrix = *M;
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
    shared_ptr<matrix_d> M; // Reference to the matrix
    int i, j;    // Indices of the matrix to update
    double sum = 0.0;
    int count = 0; // Ex v1.size() = 2, count stop at 2
    size_t limit;
    Sink(shared_ptr<matrix_d> M, int i, int j, size_t) : M(M), i(i), j(j), limit(limit) {
        //cout << "Indirizzo di M nel Sink: " << M.get() << endl;
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

    void svc_end(){ // matrix_d, int(i), int(j)
        lock_guard<mutex> lock(mtx);
        // Update the matrix
        const double element = cbrt(sum);
        (*M)[i][j] = element;
        (*M)[j][i] = element;
        //cout << "Ho aggiornato la matrice con il valore: " << element << " nella posizione M[" << i << "][" << j << "]" << endl;
    }
};

/*!
    \name DotProduct_Emitter
    \brief Split the vectors v1 and v2 for the workers
    \note Split the vectors v1 and v2 for the workers
*/
struct DotProduct_Emitter: ff_node_t<tuple<vector_d, vector_d, int>>{
    uint16_t K, W, D;
    const vector_d& v1, v2;
    
    DotProduct_Emitter(uint16_t K, uint16_t W, uint16_t D, const vector_d& v1, const vector_d& v2)
        : K(K), W(W), D(D), v1(v1), v2(v2) {}
    
    tuple<vector_d, vector_d, int>* svc(tuple<vector_d, vector_d, int>*){
        //
        list<vector_d> sub_v1_list(D);
        list<vector_d> sub_v2_list(D);
        //
        // if K is less or equal with D we can have 1:(1,1) 
        // [one workers:for one elememt] for calculate the dot product
        sub_v1_list = SplitVector(v1, D, K);
        sub_v2_list = SplitVector(v2, D, K);
        //
        auto begin_v1 = sub_v1_list.begin();
        auto begin_v2 = sub_v2_list.begin();

        uint16_t limit = min((size_t)D, sub_v1_list.size());

        for(uint16_t w = 0; w < limit; w++, begin_v1++, begin_v2++){
            // Take the sub-vectors and send them to the workers
            const uint16_t size = begin_v1->size(); // v1 and v2 have the same size
            ff_send_out(new tuple<vector_d, vector_d, int>(*begin_v1, *begin_v2, size));
            //cout << "Tuple sent (DotProduct)" << endl;
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
    shared_ptr<matrix_d> M;
    uint16_t N, K, W, D;

    Diagonal_Emitter(shared_ptr<matrix_d> M, uint16_t N, uint16_t K, uint16_t W, uint16_t D)
        : M(M), N(N), K(K), W(W), D(D) {}

    tuple_dot_product* svc(tuple_dot_product*){
        //cout << "m = [0, " << N - K << "[" << endl;
        //matrix_d& matrix = *M;
        for (int m = 0; m < N - K; m++){
            //cout << "Taking v1 and v2 vectors for m: " << m << endl;
            vector_d v1;
            vector_d v2;
            cout << fixed << showpoint;
            cout << setprecision(6);
            for (int i = 0; i < K; i++){
                v1.push_back((*M)[m][m+i]);
                //cout << "M[" << m << "][" << m+i << "]: " << (*M)[m][m+i] << endl;
                v2.push_back((*M)[m+i+1][m+K]);
                //cout << "M[" << m+i+1 << "][" << m+K << "]: " << (*M)[m+i+1][m+K] << endl;
            }
            // cout << "v1: { ";
            // for(auto i : v1){
            //     cout << i << " ";
            // }
            // cout << "}" << endl;

            // cout << "v2: { ";
            // for(auto i : v2){
            //     cout << i << " ";
            // }
            // cout << "}" << endl;
            // cout << "Sending the tuple to the farm" << endl;
            ff_send_out(new tuple_dot_product(v1, v2, M, K, W, D, m, m+K));
            //cout << "Tuple sent" << endl;
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
struct DotProduct_Worker: ff_node_t<tuple<vector_d, vector_d, int>, double>{
    double* svc(tuple<vector_d, vector_d, int> *task){
        const vector_d& v1 = std::get<0>(*task);
        const vector_d& v2 = std::get<1>(*task);
        const int size = std::get<2>(*task);

        double partial_result = PartialDotProduct(v1, v2, size);
        //cout << "Partial result: " << partial_result << endl;
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
        const vector_d& v1 = std::get<0>(*task);
        const vector_d& v2 = std::get<1>(*task);
        shared_ptr<matrix_d> M = std::get<2>(*task);
        const uint16_t K = std::get<3>(*task); // K
        const uint16_t W = std::get<4>(*task); // Number of workers
        const uint16_t D = std::get<5>(*task); // Number of workers for the DotProduct
        const int i = std::get<6>(*task); // i
        const int j = std::get<7>(*task); // j


        vector<unique_ptr<ff_node>> workers;

        for(uint16_t w = 0; w < D; w++){
            workers.push_back(make_unique<DotProduct_Worker>());
        }

        // Matrix M is shared between the workers
        const uint16_t v1_size = v1.size();
        Sink sink(M, i, j, v1_size);                             // Create the sink
        DotProduct_Emitter emitter(K, W, D, v1, v2); // Create the emitter
        ff_Farm<tuple<vector_d, vector_d, int>> farm(move(workers));
        // Add the emitter and the sink
        farm.add_emitter(emitter);
        farm.add_collector(sink);
        //cout << "Farm created (DotProduct)" << endl;

        //farm.wrap_around();
        farm.set_scheduling_ondemand();
        if (farm.run_and_wait_end() < 0){
            error("Running farm (DotProduct)\n");
            return EOS;
        }
        return GO_ON;
    }
};

struct M_Diagonal_Stage: ff_node_t<shared_ptr<matrix_d>, shared_ptr<matrix_d>> { // Take matrix M and return M'

private:
    uint16_t N, K, W;
    uint16_t Z, D;
public:
    M_Diagonal_Stage(uint16_t N, uint16_t K, uint16_t W, uint16_t Z, uint16_t D): N(N), K(K), W(W), Z(Z), D(D) {}

    shared_ptr<matrix_d> *svc(shared_ptr<matrix_d> *M) {
        
        vector<unique_ptr<ff_node>> workers;
        // Cycle with m with m = [0, n-k[
        for (uint16_t w = 0; w < Z; w++){
            //cout << "Worker - M: " << w << endl;
            workers.push_back(make_unique<DotProduct_Stage>());
        }
        Diagonal_Emitter emitter(*M, N, K, W, D);
        ff_Farm<tuple_dot_product, shared_ptr<matrix_d>> farm(move(workers), emitter);
        //cout << "Farm created (M-Diagonal)" << endl;
        farm.wrap_around();
        //farm.remove_collector();
        farm.set_scheduling_ondemand();
        
        if (farm.run_and_wait_end() < 0){
            error("Running farm (M-Diagonal)\n");
            return nullptr;
        }

        // cout << "Farm runned (M-Diagonal)" << endl;

        // for (uint16_t m = 0 ; m < N - K; m++){
        //     cout << "Taking v1 and v2 vectors for m: " << m << endl;
        //     vector_d v1;
        //     vector_d v2;
        //     // Take the vectors v1 and v2 starting from position m and m+K to K
        //     // Fill the v1 vector
        //     for (uint16_t i = 0; i < K; i++){
        //         v1.push_back((*M)[m][i]);
        //         cout << "M[" << m << "][" << i << "]: " << (*M)[m][i] << endl;
        //     }
        //     for (uint16_t i = 0; i < K; i++){
        //         v2.push_back((*M)[m+K][i]);
        //         cout << "M[" << i << "][" << m+K << "]: " << (*M)[i][m+K] << endl;
        //     }
        //     cout << "Sending the tuple to the farm" << endl;
        //     auto task = new tuple_dot_product(v1, v2, *M, K, W, D, m, m+K);
        //     farm.ff_send_out(task);
        //     cout << "Tuple sent" << endl;
        // }
        // farm.wait();
        // delete M;
        //matrix_d& matrix = **M;
        ff_send_out(M);
        return EOS;
    }
};

struct SaveMatrix_Stage: ff_node_t<shared_ptr<matrix_d>, void>{
    void *svc(shared_ptr<matrix_d> *task) {
        //matrix_d& matrix = **task;
        uint16_t size = (*task)->size();
        //SaveMatrixToFile(*task, size, "matrix_prime.txt");
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
        //cout << "K = " << k << endl;
        //cout << "DiagonalStage-Workers Z: " << resources.Z << " DotProductStage-Workers D: " << resources.D << endl;
        //cout << endl;
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


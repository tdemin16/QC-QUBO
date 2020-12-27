#include "../lib/generators.h"
#include "../lib/lib.h"

using namespace std;
using namespace Eigen;

/*
Use ./deploy.sh to genereate an annealer executable
Use ./simulate.sh to generate a classical eecutable
These two will compile in 64 bit architecture using g++ -std=c++2a

BINARY computes using vectors of {0, 1}^n
SPIN computes using vecotrs of {-1, 1}^n

Use solve() to compute a QUBO problem
- Input:
    MatrixXf Q that contains a QUBO problem
    int number of iteration
    bool mode (BINARY or SPIN) -> SPIN is not required for NPP and QAP
    int number of annealers run (effective only on annealer's executable code)
    bool log
- Output
    VectorXf containing the solution

NPP class contains number partitioning problem generator and solution mapper
QAP class contains quadratic assignment problem generator and solution 

NPP::number_partitioning_problem will generate a NPP problem
- Input:
    Empty MatrixXf Q
    Empty vector<int> nums (which dimension will be used to gen the sequence)
    int range (will be generated numbers in the range [1, range])
- Output:
    float c which is the sum over nums
NPP::diff
- Input:
    Q
    Solution vector 
    c
- Output:
    Difference btw the two most similar subsets

QAP::quadratic_assignment_problem will generate a QAP problem
- Input:
    Empty MatrixXf Q
    strig containing input file ("../test/qapk.txt")
- Output:
    float penalty
QAP::y
- Input:
    Q
    Solution vector
    penalty
*/

#define RANGE 10  // Will generate numbers in [1, RANGE]
#define N 10      // Problem dimension
#define IT 1000   // Algorithm iteration
#define K 5       // Annealer's run
#define LOG false  // Log true/false

int main() {
    // Start timer
    auto start = chrono::steady_clock::now();

    MatrixXf Q;                                                  // This will contain the QUBO problem
    vector<int> nums(N);                                         // vector used to generate npp
    float c = NPP::number_partitioning_problem(Q, nums, RANGE);  // npp generation
    VectorXf x = solve(Q, IT, BINARY, K, LOG);                   // Compute solution
    float diff = NPP::diff(Q, x, c);                             // map back the solution

    auto end = chrono::steady_clock::now();  // end timer
    chrono::duration<double> difference = end - start;

    cout << difference.count() << endl;

    NPP::to_file(difference, N, RANGE, diff, x, nums);

    return 0;
}
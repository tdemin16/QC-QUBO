#include "../lib/lib.h"
#include "../lib/npp.h"
#include "../lib/qap.h"
#include "../lib/tsp.h"

using namespace std;
using namespace Eigen;

/*
Use ./deploy.sh to genereate an annealer executable
Use ./simulate.sh to generate a classical excutable
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
    int max coefficient -> will store the max coefficient, not necessary used for statistics
    float lambda factor -> used to compute penalty
    string containing input file ("../test/qapk.txt")
- Output:
    float penalty
QAP::y
- Input:
    Q
    Solution vector
    penalty
*/

#define IT 4  // Algorithm iteration
#define N 10
#define K 5       // Number of measurements per problem
#define RANGE 1000000
#define LOG true  // Log true/false
#define LOG_FILE true

int main() {
    string filename = to_string(time(0));
    // Start timer
    auto start = chrono::steady_clock::now();

    MatrixXd Q;  // This will contain the QUBO problem
    vector<ll> nums(N);

    ll c = NPP::number_partitioning_problem(Q, nums, RANGE);

    VectorXd x(N);
    x = solve(Q, IT, BINARY, K, LOG, LOG_FILE, filename);

    auto end = chrono::steady_clock::now();  // end timer
    chrono::duration<double> difference = end - start;

    cout << endl
         << difference.count() << "s" << endl;


    ll diff = NPP::diff(Q, x, c);

    NPP::to_file(difference, IT, fQ(Q, x), N, RANGE, diff, x, nums, filename);

    return 0;
}
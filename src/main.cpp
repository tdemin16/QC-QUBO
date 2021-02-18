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

#define IT 2000   // Algorithm iteration
#define K 5       // Annealer's run
#define LOG true  // Log true/false

int main() {
    string filename = to_string(time(0));
    // Start timer
    auto start = chrono::steady_clock::now();

    MatrixXd Q; // This will contain the QUBO problem
    double max_coeff;
    float lambda = 2.25;
    //C.E. Nugent, T.E. Vollmann and J. Ruml
    //Y. Li and P.M. Pardalos
    string file = "../test/nug16a.txt";
    double pen = QAP::quadratic_assignment_problem(Q, max_coeff, lambda, file);
    VectorXd x = solve(Q, IT, BINARY, K, LOG, filename); // Compute solution
    double y = QAP::y(Q, x, pen);

    auto end = chrono::steady_clock::now(); // end timer
    chrono::duration<double> difference = end - start;

    cout << difference.count() << endl;

    QAP::to_file(difference, IT, K, lambda, pen, fQ(Q, x), file, y, x, filename);

    return 0;
}
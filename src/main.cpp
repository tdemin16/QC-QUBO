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

#define IT 50  // Algorithm iteration
#define N 4
#define K 5  // Number of measurements per problem
#define RANGE 5

int main() {
    const string filename = to_string(time(0));
    // Start timer
    const auto start = chrono::steady_clock::now();

    MatrixXd Q;  // This will contain the QUBO problem

    vector<Point> points = {Point(6.72128556, 3.78364346),
                            Point(5.79083958, 3.90574701),
                            Point(5.8849624 , 7.81206684),
                            Point(5.24617768, 5.25306609),
                            Point(9.27320268, 8.78431339),
                            Point(4.20709519, 8.69949136),
                            Point(2.46245669, 8.24656887),
                            Point(4.4052503 , 3.75829681),
                            Point(7.30407378, 6.13043625),
                            Point(8.35143241, 2.43962078)};

    MatrixXd D = TSP::build_tsp(points);

    TSP::travelling_salesman_problem(Q, D, points.size());

    VectorXd x;
    x = solve(Q, IT, K, LOG_CONSOLE | LOG_FILE, filename);

    const auto end = chrono::steady_clock::now();  // end timer
    const chrono::duration<double> difference = end - start;

    cout << endl
         << difference.count() << "s" << endl;

    vector<ll> solution;
    solution = TSP::decode_solution(x, true);
    cout << "QALS solution" << endl << "QALS: [";
    for(auto it:solution) {
        cout << it << " ";
    }
    cout << "] " << TSP::cost_route(D, solution) << endl;
    cout << "Calculation time: "<< difference.count() << endl;

#ifdef SIMULATION
    cout << "BRUTE:" << TSP::tsp_brute(D) << endl;
#endif

    return 0;
}
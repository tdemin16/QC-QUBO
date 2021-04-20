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

#define IT 300  // Algorithm iteration
#define K 5  // Number of measurements per problem

int main() {
    const string filename = to_string(time(0));
    // Start timer
    auto start = chrono::steady_clock::now();

    MatrixXd Q;  // This will contain the QUBO problem

    vector<Point> points = {Point(8.8580261 , 9.19800721),
                            Point(3.58473975, 1.65934421),
                            Point(4.6050223 , 0.72658109),
                            Point(0.01825975, 3.78999345),
                            Point(0.16127204, 5.97264417),
                            Point(7.03489676, 7.40024074),
                            Point(0.86655297, 0.49004868),
                            Point(0.35013547, 7.08491828),
                            Point(7.02722521, 5.99131315),
                            Point(4.32180413, 1.06361234),
                            Point(3.81010642, 7.63817388),
                            Point(4.45345072, 4.45003767)};

    MatrixXd D = TSP::build_tsp(points);

    TSP::travelling_salesman_problem(Q, D, points.size());

    VectorXd x;
    x = solve(Q, IT, K, filename, LOG_CONSOLE);

    auto end = chrono::steady_clock::now();  // end timer
    chrono::duration<double> difference = end - start;

    cout << endl
         << difference.count() << "s" << endl;

    vector<ll> solution;
    solution = TSP::decode_solution(x, true);
    cout << "QALS solution" << endl
         << "QALS: [";
    for (lu it = 0; it < solution.size(); it++) {
        cout << solution[it];
        if (it < solution.size() - 1) cout << ", ";
    }
    cout << "] " << TSP::cost_route(D, solution) << endl;
    cout << "Calculation time: " << difference.count() << endl;

    if (solution.size() <= 16) {
        start = chrono::steady_clock::now();
        solution = TSP::tsp_brute(D);
        end = chrono::steady_clock::now();
        difference = end - start;

        cout << "Brute Force solution" << endl
             << "Brute Force: [";
        for (lu it = 0; it < solution.size(); it++) {
            cout << solution[it];
            if (it < solution.size() - 1) cout << ", ";
        }
        cout << "] " << TSP::cost_route(D, solution) << endl;
        cout << "Calculation time: " << difference.count() << endl;
    }

    return 0;
}
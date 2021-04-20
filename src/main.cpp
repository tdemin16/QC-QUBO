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

#define IT 150  // Algorithm iteration
#define K 5  // Number of measurements per problem

int main() {
    const string filename = to_string(time(0));
    // Start timer
    auto start = chrono::steady_clock::now();

    MatrixXd Q;  // This will contain the QUBO problem

    vector<Point> points = {Point(7.21362447, 4.72776071),
                            Point(7.48093569, 8.15635289),
                            Point(9.72339245, 7.86815924),
                            Point(2.68741026, 2.43794644),
                            Point(6.86692814, 2.86802978),
                            Point(8.11962058, 4.72617735),
                            Point(9.66730978, 1.6601765 ),
                            Point(0.48265261, 9.97703727),
                            Point(4.85291697, 3.57640542),
                            Point(4.64010762, 0.77011666)};

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
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
#define K 5     // Number of measurements per problem
#define N_TESTS 10

int main() {
    chrono::_V2::steady_clock::time_point start;
    chrono::_V2::steady_clock::time_point end;
    chrono::duration<double> difference;

    const string filename = to_string(time(0));

    MatrixXd Q;  // This will contain the QUBO problem

    vector<Point> points = {Point(0.66083966, 9.36939755),
                            Point(1.98485729, 7.62491491),
                            Point(1.54206421, 1.18410071),
                            Point(2.01555644, 3.15713817),
                            Point(7.83888128, 8.77009394),
                            Point(1.4779611, 4.16581664),
                            Point(0.6508892, 6.31063212),
                            Point(6.6267559, 5.45120931),
                            Point(9.73821452, 2.20299234),
                            Point(3.50140032, 5.36660266)};

    MatrixXd D = TSP::build_tsp(points);

    TSP::travelling_salesman_problem(Q, D, points.size());

    VectorXd x;
    vector<ll> solution;
    for (int i = 0; i < N_TESTS; i++) {
        //start = chrono::steady_clock::now();

        //x = solve(Q, IT, K, filename);

        //end = chrono::steady_clock::now();  // end timer
        //difference = end - start;
        //cout << endl
        //     << difference.count() << "s" << endl;

        solution = TSP::decode_solution(x, true);
        cout << "QALS with validation" << endl
             << "QALS: [";
        for (lu it = 0; it < solution.size(); it++) {
            cout << solution[it];
            if (it < solution.size() - 1) cout << ", ";
        }
        cout << "] " << TSP::cost_route(D, solution) << endl;

        solution = TSP::decode_solution(x, false);
        cout << "QALS solution" << endl
             << "QALS: [";
        for (lu it = 0; it < solution.size(); it++) {
            cout << solution[it];
            if (it < solution.size() - 1) cout << ", ";
        }
        cout << "] " << TSP::cost_route(D, solution) << endl;
        cout << "Calculation time: " << difference.count() << endl
             << endl;
    }

    //if (solution.size() <= 16) {
    //    start = chrono::steady_clock::now();
    //    solution = TSP::tsp_brute(D);
    //
    //    end = chrono::steady_clock::now();
    //    difference = end - start;
    //
    //    cout << endl
    //         << "Brute Force solution" << endl
    //         << "Brute Force: [";
    //    for (lu it = 0; it < solution.size(); it++) {
    //        cout << solution[it];
    //        if (it < solution.size() - 1) cout << ", ";
    //    }
    //    cout << "] " << TSP::cost_route(D, solution) << endl;
    //    cout << "Calculation time: " << difference.count() << endl;
    //}

    return 0;
}
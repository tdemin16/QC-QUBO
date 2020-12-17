#include "../lib/generators.h"
#include "../lib/lib.h"

using namespace std;
using namespace Eigen;

/*
BINARY computes using vectors of {0, 1}^n
SPIN computes using vecotrs of {-1, 1}^n

Use solve() to compute a QUBO problem
*/

int main() {
    MatrixXf Q;
    float c = QAP::quadratic_assignment_problem(Q, "../test/qap3.txt");
    VectorXf x = solve(Q, 10, BINARY, 5, true);
    cout << QAP::y(Q, x, c) << endl;

    return 0;
}
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
    float c = NPP::number_partitioning_problem(Q, "../test/npp16.txt");
    VectorXf x = solve(Q, 1000, BINARY, 5, true);
    cout << NPP::diff(Q, x, c) << endl;

    return 0;
}
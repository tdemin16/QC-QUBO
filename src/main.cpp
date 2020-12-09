#include "../lib/lib.h"
#include "../lib/generators.h"

using namespace std;
using namespace Eigen;

/*
BINARY computes using vectors of {0, 1}^n
SPIN computes using vecotrs of {-1, 1}^n
*/

int main() {
    MatrixXf Q;  // QUBO Problem Matrix
    
    int c = number_partitioning_problem(Q, "../test/npp8.txt");

    VectorXf sol = solve(Q, BINARY);

    cout << sol.transpose() << endl;
    
    return 0;
}
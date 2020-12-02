#include "../lib/lib.h"
#include "../lib/generators.h"

using namespace std;
using namespace Eigen;

int main() {
    MatrixXf Q;  // QUBO Problem Matrix
    
    int c = number_partitioning_problem(Q);

    VectorXf sol = solve(Q);

    cout << "z*: " << sol.transpose() << endl;

    return 0;
}
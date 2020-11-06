#include <time.h>

#include <iostream>
#include <vector>

#include "Eigen/Core"
#include "Eigen/SparseCore"
#include "lib.h"
#include "randutils.hpp"
#include "unsupported/Eigen/KroneckerProduct"

using namespace std;

#define n 16

mt19937 e_uniform_g;
mt19937 e_uniform_h;
mt19937 e_uniform_shuffle;
mt19937 e_uniform_ann;
mt19937 e_uniform_pert;
mt19937 e_uniform_vector;
unsigned int seed_g;
unsigned int seed_h;
unsigned int seed_shuffle;
unsigned int seed_ann;
unsigned int seed_pert;
unsigned int seed_vector;
random_device rd;
uniform_real_distribution<float> d_real_uniform(0.0, 1.0);
uniform_int_distribution<int> d_int_uniform(0, 2048);  // 2048 max number of nodes

int main() {
    srand(time(0));
    init_seeds();

    MatrixXf Q(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            if (j == i)
                Q(i, j) = (float)i / 2;  // i/2 on the diagonal (after will be added up by itself)
            else
                Q(i, j) = (float)((rand() % 201) - 100) / 10;  // random value in [1, 10] on half the matrix
        }
    }
    MatrixXf temp = Q;
    temp.transposeInPlace();
    Q += temp;  // Add itself but transposed -> symmetric matrix

    SparseMatrix<float> A = init_A(n);
    
    SparseMatrix<float> theta2;
    theta2 = g_strong(Q, A, 0.75);

    return 0;
}
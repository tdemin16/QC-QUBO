#include <cstdio>
#include <iostream>
#include <vector>

#include "Eigen/Core"
#include "Eigen/SparseCore"
#include "lib.h"
#include "randutils.hpp"
#include "unsupported/Eigen/KroneckerProduct"

using namespace std;
using namespace Eigen;

VectorXf solve(MatrixXf);

mt19937 e_uniform_g;
mt19937 e_uniform_h;
mt19937 e_uniform_ann;
mt19937 e_uniform_pert;
unsigned int seed_g;
unsigned int seed_h;
unsigned int seed_ann;
unsigned int seed_pert;
random_device rd;
uniform_real_distribution<float> d_real_uniform(0.0, 1.0);

int main() {
    srand(time(0));
    int n = 16;        // number of coefficients (has to be 8*x)
    MatrixXf Q(n, n);  //QUBO Problem Matrix
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            if (j == i)
                Q(i, j) = (float)i / 2;  // i/2 on the diagonal (after will be added up by itself)
            else
                Q(i, j) = (float)((rand() % (100)) + 10) / 10;  // random value in [1, 10] on half the matrix
        }
    }
    MatrixXf temp = Q;
    temp.transposeInPlace();
    Q += temp;  // Add itself but transposed -> symmetric matrix

    solve(Q);

    return 0;
}

VectorXf solve(MatrixXf Q) {
    //Init
    int n = Q.rows();

    SparseMatrix<float> A = init_A(n);  //Chimera topology

    //Input
    float pmin = 0.1f;     // minimum probability 0 < pδ < 0.5 of permutation modification
    float eta = 0.1f;      // probability decreasing rate η > 0
    float q = 0.1f;        // candidate perturbation probability q > 0
    float lambda0 = 1.0f;  // initial balancing factor λ0 > 0
    int k = 1;             // number of annealer runs k ≥ 1
    float p = 1.0f;        // probability of an element to be considered for shuffling
    int N = 10;            // Decreasing time

    //Termination Parameters
    int imax = 1000;
    int Nmax = 50;
    int dmin = 50;

    MatrixXf In(n, n);  //Identity matrix
    In.setIdentity();

    printf("pmin:%f, eta:%f, q:%f, lambda0: %f\n", pmin, eta, q, lambda0);
    printf("k:%d\n", k);
    printf("imax:%d, N:%d, dmin:%d\n", imax, N, dmin);
    cout << endl
         << "Q" << endl
         << Q << endl
         << endl;
    cout << "A" << endl
         << A << endl
         << endl;

    //Algorithm
    MatrixXf Q_first(n, n);
    MatrixXf P(n, n), P1(n, n), P2(n, n), P_star(n, n);
    SparseMatrix<float> theta1(n, n), theta2(n, n), theta_first;
    VectorXf z_star(n), z_first(n), z1(n), z2(n);
    MatrixXf z_diag(n, n);
    MatrixXf S(n, n);  //Tabu Matrix
    float f1, f2, f_star, f_first;
    int e = 0;
    int d = 0;
    float lambda = lambda0;
    bool perturbed = false;  // True when h perturbs the candidate

    P = In;

    cout << endl
         << "---Test g---" << endl
         << "P" << endl
         << P << endl
         << endl
         << g(P, p) << endl
         << endl
         << "---End Test g, ok if P is permuted---" << endl
         << endl;

    cout << endl
         << "---Test shuffle---" << endl;

    map<int, int> m;
    for (int i = 0; i < n; i++) {
        if (d_real_uniform(e_uniform_g) <= 0.5f) {  // checks if the extracted numeber is equal or less than pr
            m.insert(pair<int, int>(i, i));         //if it is, inserts the number
        }
    }
    for (auto it : m) cout << "{" << it.first << "}->{" << it.second << "}" << endl;
    cout << endl;
    shuffle(m);
    for (auto it : m) cout << "{" << it.first << "}->{" << it.second << "}" << endl;

    cout << endl
         << "---End Test shuffle, ok if m is permuted---" << endl
         << endl
         << endl;

    cout << "---Test fQ---" << endl
         << endl;
    z1 << -1, 1, -1, 1, -1, 1, -1, -1, 1, 1, -1, 1, -1, -1, 1, 1;
    cout << "z=" << z1.transpose() << endl
         << "fQ=" << fQ(Q, z1) << endl
         << endl
         << "---End Test fq, ok if z^T * Q * z" << endl
         << endl
         << endl;

    cout << "---Test h---" << endl
         << endl
         << z1.transpose() << endl;
    h(z1, 0.5f);
    cout << z1.transpose() << endl
         << endl
         << "---End Test h, ok if z1 is half permuted permuted---" << endl
         << endl
         << endl;
         
    z2 = z1;
    cout << "---Test vector comparison---" << endl
         << endl
         << z1.transpose() << endl
         << z1.transpose() << endl
         << comp_vectors(z1, z1) << endl
         << endl
         << z1.transpose() << endl;
    h(z2, 0.5f);
    cout << z2.transpose() << endl
         << comp_vectors(z1, z2)
         << endl
         << "---End Test comparison---" << endl
         << endl
         << endl;

    return VectorXf();
}
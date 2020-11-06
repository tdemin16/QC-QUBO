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

#define n 16 // number of coefficients (has to be 8*x)

VectorXf solve(MatrixXf);

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
    MatrixXf Q(n, n);  // QUBO Problem Matrix

#ifndef DEBUG
    srand(time(0));
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
#else
    Q << 0,  0.2,  9.7,  9.5, -4.6, -8.2,  4.7, -3.4, -3.5,    6,  0.2,  6.3,  8.6,  8.7,  7.8,  3.9,
       0.2,    1, -0.1, -7.1, -5.5,  9.3,  6.8, -7.7,  6.7,   -3, -3.5, -6.9,    9,  4.6, -6.2, -0.4,
       9.7, -0.1,    2, -1.9, -0.3,  9.8, -2.3, -0.9, -9.9, -5.5, -6.2, -3.2,    1,  4.8,    2, -7.7,
       9.5, -7.1, -1.9,    3,  3.3,  0.7,  5.1, -7.9, -9.5,    8, -8.3,  9.8,  4.7,  9.1,  6.4,  6.7,
      -4.6, -5.5, -0.3,  3.3,    4, -4.5,  4.6,  5.7,  -10,  3.4, -4.8, -1.9,   -7,   10,  0.8,  2.2,
      -8.2,  9.3,  9.8,  0.7, -4.5,    5, -9.9,  5.4,    6, -8.1, -8.7,  0.7,  3.9, -6.4,    9, -5.5,
       4.7,  6.8, -2.3,  5.1,  4.6, -9.9,    6,  3.7, -8.9,  -10,  1.6,  7.9,  4.8, -8.8,  6.9,  1.2,
      -3.4, -7.7, -0.9, -7.9,  5.7,  5.4,  3.7,    7,    8, -7.7, -9.3, -1.4,  7.4,  4.1,  3.8, -9.5,
      -3.5,  6.7, -9.9, -9.5,  -10,    6, -8.9,    8,    8,  2.1,  3.7,  1.3, -5.8, -1.2, -8.4,  5.2,
         6,   -3, -5.5,    8,  3.4, -8.1,  -10, -7.7,  2.1,    9,  0.7,  8.1, -4.2,  9.7,  6.7,  9.9,
       0.2, -3.5, -6.2, -8.3, -4.8, -8.7,  1.6, -9.3,  3.7,  0.7,   10,  9.2,  0.4,    6,  9.3,    7,
       6.3, -6.9, -3.2,  9.8, -1.9,  0.7,  7.9, -1.4,  1.3,  8.1,  9.2,   11,  3.8,    4,  8.3,  0.6,
       8.6,    9,    1,  4.7,   -7,  3.9,  4.8,  7.4, -5.8, -4.2,  0.4,  3.8,   12, -9.9,  1.2, -2.1,
       8.7,  4.6,  4.8,  9.1,   10, -6.4, -8.8,  4.1, -1.2,  9.7,    6,    4, -9.9,   13,  5.9,  9.8,
       7.8, -6.2,    2,  6.4,  0.8,    9,  6.9,  3.8, -8.4,  6.7,  9.3,  8.3,  1.2,  5.9,   14, -9.8,
       3.9, -0.4, -7.7,  6.7,  2.2, -5.5,  1.2, -9.5,  5.2,  9.9,    7,  0.6, -2.1,  9.8, -9.8,   15;
#endif

    VectorXf sol = solve(Q);

    cout << "z*: " << sol.transpose() << endl;

    return 0;
}

VectorXf solve(MatrixXf Q) {
    //Init

    if (n % 8 != 0 || n == 0) {
        cout << "[Warning] n must be multiple of 8" << endl;
        exit(1);
    }

    SparseMatrix<float> A = init_A(n);  //Chimera topology
    init_seeds();

    //Input
    float pmin = 0.1f;     // minimum probability 0 < pδ < 0.5 of permutation modification
    float eta = 0.01f;     // probability decreasing rate η > 0
    float q = 0.1f;        // candidate perturbation probability q > 0
    float lambda0 = 1.0f;  // initial balancing factor λ0 > 0
    int k = 1;             // number of annealer runs k ≥ 1
    int N = 20;            // Decreasing time

    //Termination Parameters
    int imax = 2000;  // Max number of iteration
    int Nmax = 50;    // Max number of: solution equal to the best one + solution worse than the best one
    int dmin = 30;    // Number of solution that are worse than the best beyond which the best solution is not valid anymore

    MatrixXf In(n, n);  //Identity matrix
    In.setIdentity();

    cout << "Q" << endl
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
    vector<int> perm(n);
    double f1, f2, f_star, f_first;
    double p;                 // probability of an element to be considered for shuffling
    int e = 0;                // Number of solutions that equal the best
    int d = 0;                // Number of solutions that are sub optimal
    double lambda = lambda0;  // Tabu weight
    bool perturbed;           // True when h perturbs the candidate
    bool simul_ann;

    P = In;
    p = 1.0f;

    P1 = g(P, p);
    P2 = g(P, p);

    theta1 = (P1.transpose() * Q * P1).cwiseProduct(A);
    theta2 = (P2.transpose() * Q * P2).cwiseProduct(A);

#ifdef SIMULATION
    z1 = P1.transpose() * min_energy(theta1);
    z2 = P2.transpose() * min_energy(theta2);
#else
    //Call annealer
#endif

    f1 = fQ(Q, z1);
    f2 = fQ(Q, z2);

    if (f1 < f2) {    // f1 is better
        z_star = z1;  // Best
        f_star = f1;
        P_star = P1;
        z_first = z2;  // Worst
    } else {           // f2 is better
        z_star = z2;
        f_star = f2;
        P_star = P2;
        z_first = z1;
    }

    // f1 and f2 are floats -> float comparison
    if (abs(f1 - f2) > __FLT_EPSILON__) {  // if(f1 != f2)
        z_diag = z_first.asDiagonal();     // Matrix where the diagonal is made by elemnt of z_first
        S = kroneckerProduct(z1, z1.transpose()) - In + z_diag;
    } else {
        S = MatrixXf::Zero(n, n);
    }

    int i = 0;
    do {
        perturbed = false;
        simul_ann = false;

        Q_first = Q + lambda * S;

        if (!(i % N)) p = p - (p - pmin) * eta;  // 0 mod N va considerato come 0?

        P = g(P_star, p);
        theta_first = (P.transpose() * Q * P).cwiseProduct(A);
#ifdef SIMULATION
        z_first = P.transpose() * min_energy(theta_first);
#else
        //Call annealer
#endif
        if (d_real_uniform(e_uniform_pert) <= q) {
            h(z_first, p);  // possibly perturb the candidate
            perturbed = true;
        }

        if (!comp_vectors(z_first, z_star)) {
            f_first = fQ(Q, z_first);
            if (f_first < f_star) {     // f_first is better
                swap(z_first, z_star);  // z_first is better
                f_star = f_first;
                P_star = P;
                e = 0;
                d = 0;
                z_diag = z_first.asDiagonal();
                S = S + kroneckerProduct(z1, z1.transpose()) - In + z_diag;

            } else {
                d++;
                if (d_real_uniform(e_uniform_ann) <= simulated_annealing(f_first, f_star, p) && f_first > f_star) {
                    swap(z_first, z_star);
                    f_star = f_first;
                    P_star = P;
                    e = 0;
                    simul_ann = true;
                }
            }

            lambda = min(lambda0, i, e);

        } else {
            e++;
        }

#ifdef SIMULATION
        log(z_star, f_star, lambda, p, e, d, perturbed, simul_ann, i + 1);
#endif
        i++;
    } while (i <= imax && (e + d < Nmax || d >= dmin));

    printf("pmin:%f\teta:%f\tq:%f\tlambda0:%f\tN:%d\n", pmin, eta, q, lambda0, N);
    printf("k:%d\n", k);
    printf("imax:%d, Nmax:%d, dmin:%d\n", imax, Nmax, dmin);
    printf("e:%d\td:%d\ti:%d", e, d, i);

    cout << endl
         << "f*: " << f_star << endl
         << endl;

    return z_star;
}
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
    int n = 8;  // number of coefficients (has to be 8*x)

    if (n % 8 == 0 && n != 0) {
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

        VectorXf sol = solve(Q);

        cout << "Solution: " << sol.transpose() << endl;
    } else {
        cout << "n must be multiple of 8" << endl;
    }

    return 0;
}

VectorXf solve(MatrixXf Q) {
    //Init
    int n = Q.rows();
    SparseMatrix<float> A = init_A(n);  //Chimera topology

    // Toy ------------------------------
    VectorXf z1(n);
    VectorXf z2(n);
    /*for (int i = 0; i < n; i++) {
        if (rand() % (2))
            z1(i) = 1;
        else
            z1(i) = -1;
    }*/
    for (int i = 0; i < n; i++) {
        if (rand() % (2))
            z2(i) = 1;
        else
            z2(i) = -1;
    }
    //-----------------------------------

    //Input
    float pmin = 0.1f;     // minimum probability 0 < pδ < 0.5 of permutation modification
    float eta = 0.1f;      // probability decreasing rate η > 0
    float q = 0.1f;        // candidate perturbation probability q > 0
    float lambda0 = 1.0f;  // initial balancing factor λ0 > 0
    int k = 1;             // number of annealer runs k ≥ 1
    float p = 1.0f;        // probability of an element to be considered for shuffling
    int N = 10;            // Decreasing time

    //Termination Parameters
    int imax = 10;
    int Nmax = 1000;
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
    VectorXf z_star(n), z_first(n);
    MatrixXf z_diag(n, n);
    MatrixXf S(n, n);  //Tabu Matrix
    float f1, f2, f_star, f_first;
    int e = 0;
    int d = 0;
    float lambda = lambda0;

    P = In;
    P1 = g(P, p);
    P2 = g(P, p);

    theta1 = (P1.transpose() * Q * P1).cwiseProduct(A);
    theta2 = (P2.transpose() * Q * P2).cwiseProduct(A);

    z1 = P1.transpose() * min_energy(theta1);
    z2 = P2.transpose() * min_energy(theta2);

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
        Q_first = Q + lambda * S;

        if (!(i % N)) p = p - (p - pmin) * eta;  // 0 mod N va considerato come 0?

        P = g(P_star, p);
        theta_first = (P.transpose() * Q * P).cwiseProduct(A);
        z_first = P.transpose() * min_energy(theta_first);

        if (d_real_uniform(e_uniform_pert) <= q) h(z_first, p);  // possibly perturb the candidate

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
                if (d_real_uniform(e_uniform_ann) <= simulated_annealing(f_first, f_star, p)) {
                    swap(z_first, z_star);
                    f_star = f_first;
                    P_star = P;
                    e = 0;
                }
            }

            lambda = min(lambda0, i, e);

        } else {
            e++;
        }

        i++;
    } while (i <= imax && (e + d < Nmax || d > dmin));
    cout << endl
         << f_star << endl
         << endl;
    return z_star;
}
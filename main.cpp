#include <cstdio>
#include <iostream>
#include <vector>

#include "Eigen/Core"
#include "lib.h"
#include "poly.h"
#include "unsupported/Eigen/KroneckerProduct"

using namespace std;
using namespace Eigen;

VectorXf solve(MatrixXf, MatrixXf, VectorXf, VectorXf, int);

int main() {
    int n = 4;  // number of coefficients

    MatrixXf Q(n, n);  //Toy Problem
    Q << 1.0f, 3.0f, 0.0f, 1.0f,
        3.0f, 4.0f, 0.0f, 2.0f,
        0.0f, 0.0f, 0.0f, 3.0f,
        1.0f, 2.0f, 3.0f, 3.0f;

    MatrixXf A(n, n);  //Toy Topology
    A << 1.0f, 0.0f, 1.0f, 0.0f,
        1.0f, 0.0f, 1.0f, 0.0f,
        1.0f, 0.0f, 1.0f, 0.0f,
        1.0f, 0.0f, 1.0f, 0.0f;

    VectorXf z1(n);  //Toy
    z1 << 1, -1, -1, 1;

    VectorXf z2(n);  //Toy
    z2 << 1, 1, -1, -1;

    solve(Q, A, z1, z2, n);

    return 0;
}

VectorXf solve(MatrixXf Q, MatrixXf A, VectorXf z1, VectorXf z2, int n) {
    //Input
    float pmin = 0.1f;     // minimum probability 0 < pδ < 0.5 of permutation modification
    float eta = 0.1f;      // probability decreasing rate η > 0
    float q = 0.1f;        // candidate perturbation probability q > 0
    float lambda0 = 1.0f;  // initial balancing factor λ0 > 0
    int k = 1;             // number of annealer runs k ≥ 1
    float p = 1.0f;        // probability of an element to be considered for shuffling
    int N = 10;            // Decreasing time

    //Termination Parameters
    int imax = 10;  // number of iterations
    int Nmax = 100;
    int dmin = 50;

    MatrixXf In(n, n);  //Identity matrix
    In.setIdentity();

    printf("pmin:%f, eta:%f, q:%f, lambda0: %f\n", pmin, eta, q, lambda0);
    printf("k:%d\n", k);
    printf("imax:%d, N:%d, dmin:%d\n", imax, N, dmin);

    MatrixXf P(n, n), P1(n, n), P2(n, n);
    P = In;

    int e = 0;
    int d = 0;
    float lambda = lambda0;

    P1 = g(P, n, p);
    P2 = g(P, n, p);

    MatrixXf theta1(n, n), theta2(n, n);

    theta1 = (P1.transpose() * Q * P1);
    theta1 = hadamard_product(theta1, A, n);

    theta2 = (P2.transpose() * Q * P2);
    theta2 = hadamard_product(theta2, A, n);

    // run annealer k times
    // estimate energy argmin P1^T and P2^T

    float f1 = fQ(Q, z1);
    float f2 = fQ(Q, z2);

    VectorXf z_star(n), z_first(n);
    float f_star;
    MatrixXf P_star(n, n);

    if (f1 < f2) {
        z_star = z1;
        f_star = f1;
        P_star = P1;
        z_first = z2;
    } else {
        z_star = z2;
        f_star = f2;
        P_star = P2;
        z_first = z1;
    }

    MatrixXf S(n, n);
    if ((f1 - f2) > __FLT_EPSILON__ ) {  //if(f1 != f2)
        S = kronecker_product(z_first, n) - In + diag(z_first, n);
    } else {
        S = MatrixXf::Zero(n, n);
    }

    int i = 0;
    MatrixXf Q_first(n, n);
    MatrixXf theta_first;
    float f_first;
    float random;
    do {
        Q_first = Q + lambda * S;
        if (!(i % N)) p = p - (p - pmin) * eta;

        P = g(P_star, n, p);
        theta_first = P.transpose() * Q * P;
        theta_first = hadamard_product(theta_first, A, n);

        // run annealer k times
        // estimate energy argmin P^T

        h(z_first, n, p);
        if (z_first != z_star) {
            f_first = fQ(Q, z_first);

            if (f_first < f_star) {
                swap(z_first, z_star);
                f_star = f_first;
                P_star = P;
                e = 0;
                d = 0;
                S = S + kronecker_product(z_first, n) - In + diag(z_first, n);
            } else {
                d++;
                random = (rand() % 100000000) / 100000000.0f;
                if (random <= simulated_annealing(f_first, f_star, p)) {
                    swap(z_first, z_star);
                    f_star = f_first;
                    P_star = P;
                    e = 0;
                }
            }

        } else
            e++;

        i++;
    } while (i < imax && (e + d < Nmax || d > dmin));  //Da riguardare

    return VectorXf();
}
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

VectorXf solve(MatrixXf, SparseMatrix<bool>);

int main() {
    srand(time(0));
    int n = 4;  // number of coefficients

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

    SparseMatrix<bool> A(n, n);  //Toy Topology
    vector<Triplet<bool>> t;
    t.reserve(8);
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (rand() % (2)) {
                t.push_back(Triplet<bool>(i, j, true));
                t.push_back(Triplet<bool>(j, i, true));
            }
        }
    }
    A.setFromTriplets(t.begin(), t.end());

    VectorXf sol = solve(Q, A);

    cout << "Solution: " << sol.transpose() << endl;

    return 0;
}

VectorXf solve(MatrixXf Q, SparseMatrix<bool> A) {
    //Init
    int n = Q.rows();

    VectorXf z1(n);  //Toy
    VectorXf z2(n);  //Toy
    for (int i = 0; i < n; i++) {
        if (rand() % (2))
            z1(i) = 1;
        else
            z1(i) = -1;
    }
    for (int i = 0; i < n; i++) {
        if (rand() % (2))
            z2(i) = 1;
        else
            z2(i) = -1;
    }

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
    int Nmax = 1000;
    int dmin = 50;

    int e = 0;
    int d = 0;
    float lambda = lambda0;

    //Random things
    random_device rd;
    unsigned int seed = rd() ^ rd();
    mt19937 uniform;
    uniform.seed(seed);
    uniform_real_distribution<float> real_distr(0.0, 1.0);

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
         << A << endl;

    //Algorithm
    MatrixXf Q_first(n, n);
    MatrixXf P(n, n), P1(n, n), P2(n, n), P_star(n, n);
    MatrixXf theta1(n, n), theta2(n, n), theta_first;
    VectorXf z_star(n), z_first(n);
    MatrixXf z_diag(n, n);
    MatrixXf S(n, n);  //Tabu Matrix
    float f1, f2, f_star, f_first;

    P = In;
    P1 = g(P, p);
    P2 = g(P, p);

    theta1 = P1.transpose() * Q * P1;
    theta1 = theta1.cwiseProduct(A.cast<float>());

    theta2 = P2.transpose() * Q * P2;
    theta2 = theta2.cwiseProduct(A.cast<float>());

    // run annealer k times
    // estimate energy argmin P1^T and P2^T

    f1 = fQ(Q, z1);
    f2 = fQ(Q, z2);

    if (f1 < f2) {  // f1 is better
        z_star = z1;
        f_star = f1;
        P_star = P1;
        z_first = z2;
    } else {  // f2 is better
        z_star = z2;
        f_star = f2;
        P_star = P2;
        z_first = z1;
    }

    // f1 and f2 are floats -> float comparison
    if ((f1 - f2) > __FLT_EPSILON__) {  // if(f1 != f2)
        z_diag = z_first.asDiagonal();
        S = kroneckerProduct(z1, z1.transpose()) - In + z_diag;
    } else {
        S = MatrixXf::Zero(n, n);
    }

    int i = 0;
    do {
        Q_first = Q + lambda * S;
        
        if (!(i % N)) p = p - (p - pmin) * eta;  // 0 mod N va considerato come 0?

        P = g(P_star, p);
        theta_first = P.transpose() * Q * P;
        theta_first = theta_first.cwiseProduct(A.cast<float>());

        // run annealer k times
        // estimate energy argmin P^T

        if (real_distr(uniform) <= q) h(z_first, p);  // possibly perturb the candidate

        if (z_first != z_star) {
            f_first = fQ(Q, z_first);
            if (f_first < f_star) {  // f_first is better than f_star
                swap(z_first, z_star);
                f_star = f_first;
                P_star = P;
                e = 0;
                d = 0;
                z_diag = z_first.asDiagonal();
                S = S + kroneckerProduct(z1, z1.transpose()) - In + z_diag;

            } else {
                d++;
                if (real_distr(uniform) <= simulated_annealing(f_first, f_star, p)) {
                    swap(z_first, z_star);
                    f_star = f_first;
                    P_star = P;
                    e = 0;
                }
            }

            lambda = min(lambda0, i, e);

        } else
            e++;

        i++;

    } while (i < imax && (e + d < Nmax || d > dmin));  //Da riguardare

    return z_star;
}
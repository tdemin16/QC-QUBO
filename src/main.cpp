#include <stdlib.h>
#include <sys/signal.h>

#include <chrono>
#include <cstdio>
#include <iostream>
#include <vector>

#include "../lib/Eigen/Core"
#include "../lib/Eigen/SparseCore"
#include "../lib/generators.h"
#include "../lib/lib.h"
#include "../lib/randutils.hpp"
#include "../lib/unsupported/Eigen/KroneckerProduct"

using namespace std;
using namespace Eigen;

VectorXf solve(MatrixXf);

mt19937_64 e_uniform_g;
mt19937_64 e_uniform_h;
mt19937_64 e_uniform_shuffle;
mt19937_64 e_uniform_ann;
mt19937_64 e_uniform_pert;
mt19937_64 e_uniform_vector;
uniform_real_distribution<double> d_real_uniform(0.0, 1.0);
uniform_int_distribution<unsigned long long> d_int_uniform(0, 2048);  // 2048 max number of nodes
pid_t child_pid;
int fd[4];

int main() {
    MatrixXf Q;  // QUBO Problem Matrix
    
    number_partitioning_problem(Q, "../test/input16.txt");

    VectorXf sol = solve(Q);

    cout << "z*: " << sol.transpose() << endl;

    return 0;
}

VectorXf solve(MatrixXf Q) {
    //Init
    int n = Q.outerSize();

    if (n % 8 != 0 || n == 0) {
        cout << "[Warning] n must be multiple of 8" << endl;
        exit(1);
    }

    /*---------------------------
        fd[READ] child read
        fd[WRITE] father write
        fd[READ + 2] father read
        fd[WRTIE + 2] child write
    ---------------------------*/

    if (pipe(fd) != 0) {
        cout << "[PIPE1 ERROR - CLOSING]" << endl;
    }
    if (pipe(fd + 2) != 0) {
        cout << "[PIPE2 ERROR - CLOSING]" << endl;
    }

    child_pid = fork();

    if (child_pid == 0) {
        init_child();

    } else if (child_pid == -1) {
        cout << "[FORK ERROR - CLOSING]" << endl;
        exit(4);
    }

    init_seeds();
    SparseMatrix<float> A = init_A(n);  //Chimera topology

    //Input
    double pmin = 0.2f;     // minimum probability 0 < pδ < 0.5 of permutation modification
    double eta = 0.01f;     // probability decreasing rate η > 0
    double q = 0.1f;        // candidate perturbation probability q > 0
    double lambda0 = 1.0f;  // initial balancing factor λ0 > 0
    int k = 1;              // number of annealer runs k ≥ 1
    int N = 20;             // Decreasing time

    //Termination Parameters
    int imax = 3000;  // Max number of iteration
    int Nmax = 50;    // Max number of solution equal to the best one + solution worse than the best one
    int dmin = 30;    // Number of solution that are worse than the best beyond which the best solution is not valid anymore

    MatrixXf In(n, n);  //Identity matrix
    In.setIdentity();

#ifdef SIMULATION
    cout << "Q" << endl
         << Q << endl
         << endl;
    cout << "A" << endl
         << A << endl
         << endl;
#endif

    //Algorithm
    MatrixXf Q_first(n, n);
    SparseMatrix<float> theta1(n, n), theta2(n, n), theta_first(n, n);
    VectorXf z_star(n), z_first(n), z1(n), z2(n), z_gold(n);
    MatrixXf z_diag(n, n);
    MatrixXf S(n, n);  //Tabu Matrix
    vector<int> perm(n), perm_star(n), perm1(n), perm2(n);
    double f1, f2, f_star, f_first, f_gold;
    double p = 1;             // probability of an element to be considered for shuffling
    int e = 0;                // Number of solutions that equal the best
    int d = 0;                // Number of solutions that are sub optimal
    double lambda = lambda0;  // Tabu weight
    bool perturbed;           // True when h perturbs the candidate
    bool simul_ann;
    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();

    // Initialization of perm vectors like an identity matrix of order n
    for (int i = 0; i < n; i++) {
        perm[i] = i;
        perm1[i] = i;
        perm2[i] = i;
    }

    // Hadamard product between a permuted Q and A
    theta1 = g_strong(Q, A, perm1, perm1, p);
    theta2 = g_strong(Q, A, perm2, perm2, p);

#ifdef SIMULATION
    double minimum = compute_Q(Q);  // Global minimum of Q, only for simulation pourposes

    z1 = map_back(min_energy(theta1), perm1);
    z2 = map_back(min_energy(theta2), perm2);
#else
    z1 = map_back(send_to_annealer(theta1), perm1);
    z2 = map_back(send_to_annealer(theta2), perm2);
#endif

    f1 = fQ(Q, z1);
    f2 = fQ(Q, z2);

    if (f1 < f2) {    // f1 is better
        z_star = z1;  // Best
        f_star = f1;
        perm_star = perm1;
        z_first = z2;  // Worst

        z_gold = z1;
        f_gold = f1;
    } else {  // f2 is better
        z_star = z2;
        f_star = f2;
        perm_star = perm2;
        z_first = z1;

        z_gold = z2;
        f_gold = f2;
    }

    // f1 and f2 are floats -> float comparison
    if (abs(f1 - f2) > __FLT_EPSILON__) {  // if(f1 != f2)
        z_diag = z_first.asDiagonal();     // Matrix where the diagonal is made by elemnt of z_first
        S = kroneckerProduct(z_first, z_first.transpose()) - In + z_diag;
    } else {
        S = MatrixXf::Zero(n, n);
    }

    int i = 1;
    do {
        start = chrono::steady_clock::now();
        perturbed = false;
        simul_ann = false;

        Q_first = Q + lambda * S;

        if (!(i % N))
            p = p - (p - pmin) * eta;  // 0 mod N va considerato come 0?

        theta_first = g_strong(Q_first, A, perm, perm_star, p);
#ifdef SIMULATION
        z_first = map_back(min_energy(theta_first), perm);
#else
        z_first = map_back(send_to_annealer(theta_first), perm);
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
                perm_star = perm;
                e = 0;
                d = 0;
                z_diag = z_first.asDiagonal();
                S = S + kroneckerProduct(z_first, z_first.transpose()) - In + z_diag;
            } else {
                d++;
                if (d_real_uniform(e_uniform_ann) <= simulated_annealing(f_first, f_star, p)) {
                    swap(z_first, z_star);
                    f_star = f_first;
                    perm_star = perm;
                    e = 0;
                    simul_ann = true;
                }
            }

            // Best solutions yet
            if (f_first < f_gold) {
                z_gold = z_first;
                f_gold = f_first;
            }

            lambda = min(lambda0, i - 1, e);
        } else {
            e++;
        }

#ifdef SIMULATION
        log(z_star, f_star, minimum, f_gold, lambda, p, e, d, perturbed, simul_ann, i);
#else
        log(f_star, f_gold, lambda, p, e, d, perturbed, simul_ann, i);
#endif
        end = chrono::steady_clock::now();
        chrono::duration<double> diff = end - start;
        cout << endl
             << diff.count() << "s" << endl
             << endl;
        i++;
    } while (i <= imax && (e + d < Nmax || d >= dmin));

    close(fd[READ]);
    close(fd[WRITE]);
    close(fd[READ + 2]);
    close(fd[WRITE + 2]);

#ifndef SIMULATION
    kill(child_pid, SIGKILL);
#endif

    printf("pmin:%f\teta:%f\tq:%f\tlambda0:%f\tN:%d\n", pmin, eta, q, lambda0, N);
    printf("k:%d\n", k);
    printf("imax:%d, Nmax:%d, dmin:%d\n", imax, Nmax, dmin);
    printf("e:%d\td:%d\ti:%d\n", e, d, i - 1);

    cout << endl
         << "f_gold: " << f_gold << endl
         << endl;

    return z_gold;
}
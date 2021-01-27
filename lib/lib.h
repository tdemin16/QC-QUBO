#ifndef LIB_H
#define LIB_H

#include <math.h>
#include <sys/signal.h>
#include <time.h>
#include <unistd.h>
#include <sys/wait.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>
#include <unordered_map>
#include <random>
#include <string>
#include <fstream>

#include "Eigen/Core"
#include "Eigen/IterativeLinearSolvers"
#include "Eigen/SparseCore"
#include "unsupported/Eigen/KroneckerProduct"

using namespace std;
using namespace Eigen;

#define READ 0
#define WRITE 1
#define BINARY 0
#define SPIN -1

extern mt19937_64 e_uniform_g;
extern mt19937_64 e_uniform_h;
extern mt19937_64 e_uniform_shuffle;
extern mt19937_64 e_uniform_ann;
extern mt19937_64 e_uniform_pert;
extern mt19937_64 e_uniform_vector;
extern uniform_real_distribution<double> d_real_uniform;
extern uniform_int_distribution<unsigned long long> d_int_uniform;
extern pid_t child_pid;
extern int fd[4];

VectorXf solve(MatrixXf Q, int imax=1000, int mode=BINARY, int k=4, bool logs=true);

void handle_sigint(int sig);

// change process code with python solver.py
void init_child(int mode, int k);

// init seed in order to generate random numbers
void init_seeds();

// Input: number of nodes
// Return: a SparseMatrix A containing Chimera's topology
void get_topology(unordered_map<int, int> &nodes, SparseMatrix<float> &edges, int n);

// Input: Matrix of order n Q, Vector of n integers z
// Return: z^T * Q * z
float fQ(MatrixXf Q, VectorXf x);

// Input: Matrix Q, SparseMatrix A, vector in which store the new permutation, vector from which generate the new permutation, probability of permutation
// Return: Matrix theta containing the product P^T * Q * P ○ A. permutation will contain the new permutation of indexes. Complexity O(nlogn)
SparseMatrix<float> g_strong(const MatrixXf &Q, const unordered_map<int, int> &nodes, const SparseMatrix<float> &edges, vector<int> &permutation, const vector<int> &old_permutation, double pr);

// Input: Map of n integer m
// Shuffle the matrix
void shuffle_map(map<int, int> &m);

// Input: Vector of integer v
// Shuffle the vector according to Fisher and Yates' algorithm. Used to shuffle the map
void shuffle_vector(vector<int> &v);

// Input: map of permuted indexes, vector of indexes not permuted
// Return: vector containing permuted and not permuted indexes in the correct position
vector<int> fill(const map<int, int> &m, const vector<int> &permutation);

// Input: Permutation vector
// Return: Vector in which values and indexes are swapped
vector<int> inverse(const vector<int> &permutation);

#ifndef SIMULATION
// Input: Theta matrix
// Output: z that minimize theta function
VectorXf send_to_annealer(const SparseMatrix<float> &theta, int n);
#endif

// Input: Vector of n integers z, probability pr
// With probability pr, the element zi = -zi
void h(VectorXf &z, double pr, int mode);

// Input: λ0, i, e
// Return: Minimum between λ0 and λ0/(2+i-e)
double min(double lambda0, int i, int e);

// Input: weight theta
// Return: VectorXf that minimize the weight function
VectorXf min_energy(const SparseMatrix<float> &theta, int mode);

// Input: weight function theta, VectorXf x = {-1, 1}^n
// Return: E(theta, x)
double E(const SparseMatrix<float> &theta, VectorXf x);

// Input: VectorXf of {-1, 1}^n v
// Treats the vector as a binary number but made of -1 and 1(instead of 0 and 1). Performs the increment
void increment(VectorXf &v, int mode);

VectorXf map_back(const VectorXf &z, const vector<int> &perm);

// Returns the probabilty of commuting
double simulated_annealing(double, double, double);

// Input: 2 VectorXf
// Return: true iff they are equal, false otherwise
bool comp_vectors(const VectorXf &z1, const VectorXf &z2);

void close_child();

#ifdef SIMULATION
double compute_Q(const MatrixXf &Q, int mode);
void log(const MatrixXf &Q, const VectorXf &z_star, double f_star, double min, double f_gold, double lambda, double p, int e, int d, bool perturbed, bool simul_ann, int i, string filename);
#else
void log(const VectorXf &z_star, double f_star, double f_gold, double lambda, double p, int e, int d, bool perturbed, bool simul_ann, int i, string filename);
#endif

#endif  //LIB_H
#ifndef LIB_H
#define LIB_H

#include <math.h>
#include <time.h>

#include <iostream>
#include <map>
#include <random>

#include "Eigen/Core"
#include "Eigen/IterativeLinearSolvers"
#include "Eigen/SparseCore"

using namespace std;
using namespace Eigen;

extern mt19937 e_uniform_g;
extern mt19937 e_uniform_h;
extern mt19937 e_uniform_shuffle;
extern mt19937 e_uniform_ann;
extern mt19937 e_uniform_pert;
extern unsigned int seed_g;
extern unsigned int seed_h;
extern unsigned int seed_shuffle;
extern unsigned int seed_ann;
extern unsigned int seed_pert;
extern random_device rd;
extern uniform_real_distribution<float> d_real_uniform;

// Input: Matrix of order n P, probability of permutation pr
// Return: Matrix of order n where some raws (with probability pr) are swapped
MatrixXf g(MatrixXf P, float pr);

// Input: Map of n integer m
// Shuffle the matrix
void shuffle(map<int, int> &m); //to be improved

int myrandom(int i);

// Input: Matrix of order n Q, Vector of n integers z
// Return: z^T * Q * z
float fQ(MatrixXf Q, VectorXf x);

// Input: Vector of n integers z, probability pr
// With probability pr, the element zi = -zi
void h(VectorXf &z, float pr);

// Returns the probabilty of commuting
float simulated_annealing(float, float, float);

// Input: λ0, i, e
// Return: Minimum between λ0 and λ0/(2+i-e)
float min(float lambda0, int i, int e);

// init seed in order to generate random numbers
void init_seeds();

// Input: 2 VectorXf
// Return: true iff they are equal, false otherwise
bool comp_vectors(VectorXf z1, VectorXf z2);

// Input: number of nodes
// Return: a SparseMatrix A containing CHimera's topology
SparseMatrix<float> init_A(int n);

// Input: weight theta
// Return: VectorXf that minimize the weight function
VectorXf min_energy(SparseMatrix<float> theta);

// Input: VectorXf of {-1, 1}^n v
// Treats the vector as a binary number but made of -1 and 1(instead of 0 and 1). Performs the increment
void increment(VectorXf &v);

// Input: weight function theta, VectorXf x = {-1, 1}^n
// Return: E(theta, x)
float E(SparseMatrix<float> theta, VectorXf x);

void log(VectorXf z_star, float f_star, float lambda, float p, int e, int d, bool perturbed, bool simul_ann, int i);

#endif  //LIB_H
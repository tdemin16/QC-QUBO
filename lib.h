#include <math.h>
#include <time.h>

#include <iostream>
#include <map>
#include <random>

#include "Eigen/Core"

using namespace std;
using namespace Eigen;

// Input: Matrix of order n P, probability of permutation pr
// Return: Matrix of order n where some raws (with probability pr) are swapped
MatrixXf g(MatrixXf P, float pr);

// Input: Map of n integer m
// Shuffle the matrix
void shuffle(map<int, int> &m);

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
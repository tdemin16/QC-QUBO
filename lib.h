#include <time.h>

#include <iostream>
#include <map>
#include <random>
#include "Eigen/Core"
#include <math.h>

using namespace std;
using namespace Eigen;

// Input: Matrix of order n P, dimensions of the matrix n, probability of permutation pr
// Return: Matrix of order n where some raws (with probability pr) are swapped
MatrixXf g(MatrixXf, int, float);

// Input: Map of n integer m
// Shuffle the matrix
void shuffle(map<int, int> &);

// Input: Matrix of order n Q, Vector of n integers z
// Return: z^T * Q * z
float fQ(MatrixXf, VectorXf);

// Input: Vector of n integers z, dimension of the vector n, probability pr
// With probability pr, the element zi = -zi
void h(VectorXf &, int, float);

// Returns the probabilty of commuting
float simulated_annealing(float, float, float);
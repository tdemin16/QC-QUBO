#include <time.h>

#include <iostream>
#include <map>
#include <random>

#include "Eigen/Core"

using namespace std;
using namespace Eigen;

// Input: permutation matrix(map), dimension of the matrix, probability
// Return: permuted matrix
map<vector<int>, int> g(map<vector<int>, int>, int, float);

// Input: empty matrix(map), dimension of the matrix
// The matrix become the identity matrix of order n
void make_identity(map<vector<int>, int> &, int);

// Input: matrix of order n(map), dimension of the matrix
// shuffle the matrix
void shuffle(map<vector<int>, int> &, int);

void print_map(map<vector<int>, int>, int);

// Input: matrix of order n(map), dimension of the matrix
// Return: the transposed matrix
map<vector<int>, int> transpose(map<vector<int>, int> &, int);
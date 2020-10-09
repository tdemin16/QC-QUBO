#include <cstdio>
#include <iostream>
#include <vector>

#include "Eigen/Core"
#include "lib.h"

using namespace std;
using namespace Eigen;

vector<int> solve(int);

int main() {
    int n = 6;  // number of coefficients

    solve(n);

    return 0;
}

vector<int> solve(int n) {
    //Input
    float pmin = 0.1f;     // minimum probability 0 < pδ < 0.5 of permutation modification
    float eta = 0.1f;      // probability decreasing rate η > 0
    float q = 0.1f;        // candidate perturbation probability q > 0
    float lambda0 = 1.0f;  // initial balancing factor λ0 > 0
    int k = 1;             // number of annealer runs k ≥ 1

    //Termination Parameters
    int imax = 10;  // number of iterations
    int N = 100;
    int dmin = 50;

    printf("pmin:%f, eta:%f, q:%f, lambda0: %f\n", pmin, eta, q, lambda0);
    printf("k:%d\n", k);
    printf("imax:%d, N:%d, dmin:%d\n", imax, N, dmin);

    MatrixXi P; // permutation matrix
    P.setIdentity(n, n); // set as identity matrix
    
    float p = 1; //probability of an element to be considered for shuffling
    int e = 0;
    int d = 0;

    P = g(P, n, p);
}
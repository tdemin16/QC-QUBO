#ifndef GENERATORS_H
#define GENERATORS_H

#include "lib.h"
#include "limits.h"

class NPP {
   public:
    //Input: Matrix Q, empty vector of dimension n that will contain the vector, range max of number gen
    static long long number_partitioning_problem(MatrixXd &Q, vector<int> &nums, int range);
    static long long diff(const MatrixXd &Q, const VectorXd &x, long long c);
    static void to_file(chrono::duration<double> difference, int it, double f, int n, int range, long long diff, const VectorXd &x, const vector<int> nums, string filename);
};

class QAP {
   public:
    static double quadratic_assignment_problem(MatrixXd &Q, double &max_coeff, float lambda, string file);
    static double y(const MatrixXd &Q, const VectorXd &x, double penalty);
    static void to_file(chrono::duration<double> difference, int it, int k, float lambda, double penalty, double f, string file, double y, const VectorXd &x, string filename);
};

#endif
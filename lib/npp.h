#ifndef NPP_H
#define NPP_H

#include "lib.h"

class NPP {
   public:
    //Input: Matrix Q, empty vector of dimension n that will contain the vector, range max of number gen
    static long long number_partitioning_problem(MatrixXd &Q, vector<long long> &nums, int range);
    static long long diff(const MatrixXd &Q, const VectorXd &x, long long c);
    static void to_file(chrono::duration<double> difference, int it, double f, int n, int range, long long diff, const VectorXd &x, const vector<int> nums, string filename);
};
#endif
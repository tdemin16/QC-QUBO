#ifndef GENERATORS_H
#define GENERATORS_H

#include "lib.h"

class NPP {
   public:
    //Input: Matrix Q, empty vector of dimension n that will contain the vector, range max of number gen
    static float number_partitioning_problem(MatrixXf &Q, vector<int> &nums, int range);
    static float diff(const MatrixXf &Q, const VectorXf &x, float c);
    static void to_file(chrono::duration<double> difference, int n, int range, float diff, const VectorXf &x);
};

class QAP {
   public:
    static float quadratic_assignment_problem(MatrixXf &Q, string file);
    static float y(const MatrixXf &Q, const VectorXf &x, float penalty);
};

#endif
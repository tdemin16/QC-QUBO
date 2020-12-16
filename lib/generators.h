#ifndef GENERATORS_H
#define GENERATORS_H

#include "lib.h"

class NPP {
   public:
    static float number_partitioning_problem(MatrixXf &Q, string file);
    static float diff(const MatrixXf &Q, const VectorXf &x, float c);
};

class QAP {
   public:
    static float quadratic_assignment_problem(MatrixXf &Q, string file);
    static float y(const MatrixXf &Q, const VectorXf &x, float penalty);
};

#endif
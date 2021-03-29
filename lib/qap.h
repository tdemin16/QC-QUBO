#ifndef QAP_H
#define QAP_H

#include "lib.h"

class QAP {
   public:
    static double quadratic_assignment_problem(MatrixXd &Q, double &max_coeff, float lambda, string file);
    static double y(const MatrixXd &Q, const VectorXd &x, double penalty);
    static void to_file(chrono::duration<double> difference, int it, int k, float lambda, double penalty, double f, string file, double y, const VectorXd &x, string filename);
};
#endif
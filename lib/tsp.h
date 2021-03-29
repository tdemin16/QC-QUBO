#ifndef TSP_H
#define TSP_H

#include "lib.h"

class TSP {
   public:
    static void travelling_salesman_problem(MatrixXd &Q, int n);

   private:
    static void add_cost_objective(map<pair<int,int>, double> &qubo, const MatrixXd &D, double cost_constant, int n);
    static void add_time_constraints(map<pair<int,int>, double> &qubo, double constraint_constant, int n);
    static void add_position_constraints(map<pair<int,int>, double> &qubo, double constraint_constant, int n);
    static void to_MatrixXd(MatrixXd &Q, const map<pair<int,int>, double> &qubo);
};

#endif
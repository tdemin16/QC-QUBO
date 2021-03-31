#ifndef TSP_H
#define TSP_H

#include "lib.h"

class TSP {
   public:
    static void travelling_salesman_problem(MatrixXd &Q, MatrixXd &D, int n);
    static vector<long long> decode_solution(const VectorXd &x, bool validate);
    static double cost_route(const MatrixXd &D, vector<long long> solution);
    static double tsp_brute(const MatrixXd &D, int s = 0);

   private:
    static void add_cost_objective(map<pair<long long, long long>, double> &qubo, const MatrixXd &D, double cost_constant, int n);
    static void add_time_constraints(map<pair<long long, long long>, double> &qubo, double constraint_constant, int n);
    static void add_position_constraints(map<pair<long long, long long>, double> &qubo, double constraint_constant, int n);
    static void to_MatrixXd(MatrixXd &Q, const map<pair<long long, long long>, double> &qubo);
};

#endif
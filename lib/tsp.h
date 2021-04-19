#ifndef TSP_H
#define TSP_H

#include "lib.h"

typedef struct Point {
    double x;
    double y;

    Point(double x, double y);
} Point;

class TSP {
    // We made this class from this Python version: https://github.com/BOHRTECHNOLOGY/quantum_tsp
    // Some changes occurred
   public:
    static MatrixXd build_tsp(vector<Point> &points);
    static void travelling_salesman_problem(MatrixXd &Q, const MatrixXd &D, int n);
    static vector<ll> decode_solution(const VectorXd &x, bool validate);
    static double cost_route(const MatrixXd &D, const vector<ll> &solution);
    static vector<ll> tsp_brute(const MatrixXd &D, int s = 0);
    static bool is_acceptable(const VectorXd &x);

   private:
    static void add_cost_objective(map<pair<ll, ll>, double> &qubo, const MatrixXd &D, double cost_constant, int n);
    static void add_time_constraints(map<pair<ll, ll>, double> &qubo, double constraint_constant, int n);
    static void add_position_constraints(map<pair<ll, ll>, double> &qubo, double constraint_constant, int n);
    static void to_MatrixXd(MatrixXd &Q, const map<pair<ll, ll>, double> &qubo);
};

#endif
#include "../lib/tsp.h"

void TSP::travelling_salesman_problem(MatrixXd &Q, int n) {
    MatrixXd D(n, n);  // Distance matrix
    map<pair<int, int>, double> qubo;
    double cost_constant = 10;  // tmp
    Q.resize(n * n, n * n);

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (i != j) D(i, j) = D(j, i) = rand() % 9 + 1;
        }
    }
    cout << D << endl;

    add_cost_objective(qubo, D, cost_constant, n);
    add_time_constraints(qubo, 400, n);
    add_position_constraints(qubo, 400, n);
    to_MatrixXd(Q, qubo);
    cout << Q << endl;
}

void TSP::add_cost_objective(map<pair<int, int>, double> &qubo, const MatrixXd &D, double cost_constant, int n) {
    int r, c;

    for (int t = 0; t < n; t++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    r = t * n + i;
                    c = (t + 1) % n * n + j;
                    qubo[make_pair(r, c)] = cost_constant * D(i, j);
                }
            }
        }
    }
}

void TSP::add_time_constraints(map<pair<int, int>, double> &qubo, double constraint_constant, int n) {
    int r, c;
    for (int t = 0; t < n; t++) {
        for (int i = 0; i < n; i++) {
            r = t * n + i;
            auto it = qubo.find(pair<int, int>(r, r));
            if (it != qubo.end())
                qubo[make_pair(r, r)] = -constraint_constant;
            else
                qubo[make_pair(r, r)] += -constraint_constant;

            for (int j = 0; j < n; j++) {
                c = t * n + 1;
                if (i != j) qubo[make_pair(r, c)] = 2 * constraint_constant;
            }
        }
    }
}

void TSP::add_position_constraints(map<pair<int, int>, double> &qubo, double constraint_constant, int n) {
    int r, c;
    for (int i = 0; i < n; i++) {
        for (int t1 = 0; t1 < n; t1++) {
            r = t1 * n + i;
            auto it = qubo.find(pair<int, int>(r, r));
            if (it != qubo.end())
                qubo[make_pair(r, r)] = -constraint_constant;
            else
                qubo[make_pair(r, r)] += -constraint_constant;

            for (int t2 = 0; t2 < n; t2++) {
                c = t2 * n + i;
                if (t1 != t2) qubo[make_pair(r, c)] = 2 * constraint_constant;
            }
        }
    }
}

void TSP::to_MatrixXd(MatrixXd &Q, const map<pair<int, int>, double> &qubo) {
    int r, c;
    for (auto it : qubo) {
        r = it.first.first;
        c = it.first.second;
        Q(r, c) = it.second;
    }
}
#include "../lib/tsp.h"

void TSP::travelling_salesman_problem(MatrixXd &Q, MatrixXd &D, int n) {
    map<pair<ll, ll>, double> qubo;
    Q.resize(n * n, n * n);
    D.resize(n, n);
    srand(time(0));
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (i != j) D(i, j) = D(j, i) = rand() % 9 + 1;
        }
    }
    cout << D << endl;

    const double B = 1;
    const double A = n * D.maxCoeff();

    cout << "A: " << A << endl;

    add_cost_objective(qubo, D, B, n);
    add_time_constraints(qubo, A, n);
    add_position_constraints(qubo, A, n);
    to_MatrixXd(Q, qubo);
}

void TSP::add_cost_objective(map<pair<ll, ll>, double> &qubo, const MatrixXd &D, double B, int n) {
    int r, c;

    for (int t = 0; t < n; t++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    r = t * n + i;
                    c = (t + 1) % n * n + j;
                    qubo[make_pair(r, c)] = B * D(i, j);
                }
            }
        }
    }
}

void TSP::add_time_constraints(map<pair<ll, ll>, double> &qubo, double A, int n) {
    int r, c;
    for (int t = 0; t < n; t++) {
        for (int i = 0; i < n; i++) {
            r = t * n + i;
            auto it = qubo.find(pair<ll, ll>(r, r));
            if (it != qubo.end())
                qubo[make_pair(r, r)] = -A;
            else
                qubo[make_pair(r, r)] += -A;

            for (int j = 0; j < n; j++) {
                c = t * n + j;
                if (i != j) qubo[make_pair(r, c)] = 2 * A;
            }
        }
    }
}

void TSP::add_position_constraints(map<pair<ll, ll>, double> &qubo, double A, int n) {
    int r, c;
    for (int i = 0; i < n; i++) {
        for (int t1 = 0; t1 < n; t1++) {
            r = t1 * n + i;
            auto it = qubo.find(pair<ll, ll>(r, r));
            if (it != qubo.end())
                qubo[make_pair(r, r)] = -A;
            else
                qubo[make_pair(r, r)] += -A;

            for (int t2 = 0; t2 < n; t2++) {
                c = t2 * n + i;
                if (t1 != t2) qubo[make_pair(r, c)] = 2 * A;
            }
        }
    }
}

void TSP::to_MatrixXd(MatrixXd &Q, const map<pair<ll, ll>, double> &qubo) {
    int r, c;
    for (auto it : qubo) {
        r = it.first.first;
        c = it.first.second;
        Q(r, c) = it.second;
    }
}

vector<ll> TSP::decode_solution(const VectorXd &x, bool validate) {
    ll n = sqrt(x.size());
    ll last_j = -1;
    vector<ll> solution(n, -1);
    unordered_set<ll> all;
    unordered_set<ll> ins;
    unordered_set<ll> result;
    vector<ll> sw(n, 0);

    for (ll i = 0; i < n; i++) {
        for (ll j = 0; j < n; j++) {
            if (x(n * i + j) == 1) {
                solution[i] = j;
                last_j = j;
            }
        }
        all.insert(i);
        if (last_j != -1) ins.insert(last_j);
        last_j = -1;
    }

    if (validate) {
        cout << "Solution not acceptable, validation occurred" << endl;
        for (auto it : all) {
            if (ins.find(it) == ins.end()) {
                result.insert(it);
            }
        }

        auto it = result.begin();
        for (int i = 0; i < n; i++) {
            if (solution[i] == -1) {
                solution[i] = *it;
                ins.insert(*it);
                it++;
            }
        }
        
        result.clear();
        for (auto it : all) {
            if (ins.find(it) == ins.end()) {
                result.insert(it);
            }
        }

        it = result.begin();
        if (ins.size() != all.size()) {
            for (int i = 0; i < n; i++) {
                if (sw[solution[i]] == 0)
                    sw[solution[i]]++;
                else {
                    solution[i] = *it;
                    sw[solution[i]]++;
                    it++;
                }
            }
        }
    }

    return solution;
}

double TSP::cost_route(const MatrixXd &D, vector<ll> solution) {
    int n = solution.size();
    ll a, b;
    double cost = 0;

    for (int i = 0; i < n; i++) {
        a = i;
        b = (i + 1) % n;
        cost += D(solution[a], solution[b]);
    }

    return cost;
}

double TSP::tsp_brute(const MatrixXd &D, int s) {
    // store all vertex apart from source vertex
    int n = D.outerSize();
    vector<int> vertex;
    for (int i = 0; i < n; i++)
        if (i != s)
            vertex.push_back(i);

    // store minimum weight Hamiltonian Cycle.
    int min_path = INT_MAX;
    do {
        // store current Path weight(cost)
        int current_pathweight = 0;

        // compute current path weight
        int k = s;
        for (long unsigned i = 0; i < vertex.size(); i++) {
            current_pathweight += D(k, vertex[i]);
            k = vertex[i];
        }
        current_pathweight += D(k, s);

        // update minimum
        min_path = min(min_path, current_pathweight);

    } while (next_permutation(vertex.begin(), vertex.end()));

    return min_path;
}

bool TSP::is_acceptable(const VectorXd &x) {
    ll n = sqrt(x.size());
    vector<ll> count(n, 0);
    for (ll i = 0; i < n; i++) {
        for (ll j = 0; j < n; j++) {
            if (x(n * i + j) == 1) {
                count[i]++;
            }
        }
    }

    for(auto it:count) {
        if(it != 1) return false;
    }
    return true;
}
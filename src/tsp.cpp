#include "../lib/tsp.h"

Point::Point(double x, double y) {
    this->x = x;
    this->y = y;
}

MatrixXd TSP::build_tsp(vector<Point> &points) {
    int n = points.size();
    double distance;
    MatrixXd D(n, n);

    for (int i = 0; i < n; i++) {
        D(i, i) = 0;
        for (int j = i + 1; j < n; j++) {
            distance = sqrt(pow(points[j].x - points[i].x, 2) + pow(points[j].y - points[i].y, 2));
            D(i, j) = D(j, i) = distance;
        }
    }

    return D;
}

void TSP::travelling_salesman_problem(MatrixXd &Q, const MatrixXd &D, int n) {
    map<pair<ll, ll>, double> qubo;
    Q.resize(n * n, n * n);

    const double B = 1;
    const double A = n * D.maxCoeff();

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
    // Def: an unavailable node is a node that is already stored 1 or more times in the solution
    //      that cannot be picked at random

    int rnd;
    ll n = sqrt(x.size());  // problem's dimension
    ll index;
    vector<ll> solution;  // problem solution
    vector<ll> indexes;
    set<ll> keep;            // unavailable nodes
    set<ll> all;             // all nodes
    set<ll> diff;            // nodes to insert
    vector<set<ll>> raw(n);  // solution preprocessing
                             // each subsequence can contain 0+ nodes
    double last = -1;

    if (!validate) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (x(n * i + j) == 1) last = j;
            }
            if (last != -1) solution.push_back(last);
            last = -1;
        }

    } else {
        // solution vector is initialize to -1
        solution.resize(n);
        for (lu i = 0; i < solution.size(); i++) {
            solution[i] = -1;
        }

        // foreach subsequence, its corresponding set will contain the position of all 1s
        // i.e. in 01100 with n = 5 -> {1, 2}
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (x(n * i + j) == 1) raw[i].insert(j);
            }
        }

        // solution[i] will get a value != 0 iff its corresponding set has 1 and only 1 value
        // solution[i] will get that value
        for (int i = 0; i < n; i++) {
            if (raw[i].size() == 1) {
                keep.insert(*(raw[i].begin()));  // keep keeps track of unavailable nodes
                solution[i] = *(raw[i].begin());
            }
            all.insert(i);  // all stores all nodes
        }

        // given the set of values that were unique in their subsequence
        // if the subsequence at position i present more than 1 element
        // then the solution will store a random value from that set that is not present in keep
        for (int i = 0; i < n; i++) {
            if (raw[i].size() > 1) {
                for (auto it : raw[i]) {
                    if (keep.find(it) == keep.end()) {
                        diff.insert(it);
                    }
                }
                if (diff.size() > 0) {
                    rnd = rand() % diff.size();
                    auto it = diff.begin();
                    advance(it, rnd);
                    solution[i] = *it;
                    keep.insert(*it); // update set of unavailable nodes
                    diff.clear();
                }
            }
        }

        // foreach possible node, if a node is stored multiple times in the solution
        // i.e if x = [0, 1, 0, 0,  1, 0, 0, 0,  0, 1, 0, 0] = [1, 0, 1]
        // one of the index of the solution, picked at random, will keep its value
        // the others will be set to -1
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (solution[j] == i) indexes.push_back(j); // vector of index that contain value i
            }
            if (indexes.size() > 1) {
                random_shuffle(indexes.begin(), indexes.end());
                index = indexes[0];
                for (auto it : indexes) {
                    it == index ? solution[it] = i : solution[it] = -1;
                }
                keep.insert(i); // update set of unavailable nodes
            }
            indexes.clear();
        }

        // diff = all - ins
        for (auto it : all) {
            if (keep.find(it) == keep.end()) {
                diff.insert(it);
            }
        }

        // each index of solution that has a value of -1
        // will be set to a random available value
        for (int i = 0; i < n; i++) {
            if (solution[i] == -1) {
                auto it = diff.begin();
                advance(it, rand() % diff.size());
                solution[i] = *it;
                diff.erase(it); // upadate set of available nodes
            }
        }
    }

    return solution;
}

double TSP::cost_route(const MatrixXd &D, const vector<ll> &solution) {
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

vector<ll> TSP::tsp_brute(const MatrixXd &D, int s) {
    // store all vertex apart from source vertex
    ll n = D.outerSize();
    vector<ll> vertex;
    vector<ll> best;
    for (ll i = 0; i < n; i++)
        if (i != s)
            vertex.push_back(i);

    // store minimum weight Hamiltonian Cycle.
    double min_path = __DBL_MAX__;
    double prev = __DBL_MAX__;
    do {
        // store current Path weight(cost)
        double current_pathweight = 0;

        // compute current path weight
        int k = s;
        for (long unsigned i = 0; i < vertex.size(); i++) {
            current_pathweight += D(k, vertex[i]);
            k = vertex[i];
        }
        current_pathweight += D(k, s);

        // update minimum
        min_path = min(min_path, current_pathweight);
        if (min_path < prev) best = vertex;
        prev = min_path;

    } while (next_permutation(vertex.begin(), vertex.end()));

    auto it = best.begin();
    best.insert(it, s);

    return best;
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

    for (auto it : count) {
        if (it != 1) return false;
    }
    return true;
}
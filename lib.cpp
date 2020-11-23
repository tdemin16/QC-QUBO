#include "lib.h"

void init_seeds() {
    random_device rd;
    seed_seq seed_g{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    seed_seq seed_h{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    seed_seq seed_shuffle{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    seed_seq seed_ann{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    seed_seq seed_pert{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    seed_seq seed_vector{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};

    e_uniform_g.seed(seed_g);
    e_uniform_h.seed(seed_h);
    e_uniform_ann.seed(seed_ann);
    e_uniform_pert.seed(seed_pert);
    e_uniform_shuffle.seed(seed_shuffle);
    e_uniform_vector.seed(seed_vector);
}

SparseMatrix<float> init_A(int n) {
    SparseMatrix<float> A(n, n);  //Chimera Topology
    vector<Triplet<float>> t;
    t.reserve(11 * n);             // |E| ~ 5n * 2 entries per edge + |V| = n => 10*n + n = 11n
    int block;                     // Current block (a block is made by the 8 nodes that form a bipartite subgraph)
    int row;                       // Current row of blocks
    for (int i = 0; i < n; i++) {  //For each node
        block = i / 8;
        row = i / (8 * 16);
        t.push_back(Triplet<float>(i, i, 1));  // Insert 1 on the diagonal
        if (i % 8 < 4) {                       // if i is in the "left" side of its subgraph
            for (int j = 4; j < 8; j++) {      // Insert a forward edge
                if (block * 8 + j < n) t.push_back(Triplet<float>(i, block * 8 + j, 1));
            }
            if (i + 8 * 16 < n) t.push_back(Triplet<float>(i, i + 8 * 16, 1));   // if i's corresponding in the following row exists (n is large enough) put a forward edge
            if (i - 8 * 16 >= 0) t.push_back(Triplet<float>(i, i - 8 * 16, 1));  // i i's corresponding in the previous row exists put a backward edge
        } else {                                                                 // otherwise, if i is in the "right" side of its subgraph
            for (int j = 0; j < 4; j++) {                                        // Insert a backward edge
                if (block * 8 + j < n) t.push_back(Triplet<float>(i, (block)*8 + j, 1));
            }
            if (row == (i + 8) / (8 * 16) && i + 8 < n) t.push_back(Triplet<float>(i, i + 8, 1));   // if i and its corresponding i+8 are in the same row, put a forward edge
            if (row == (i - 8) / (8 * 16) && i - 8 >= 0) t.push_back(Triplet<float>(i, i - 8, 1));  // if i and its corresponding i-8 are in the same row, put a backward edge
        }
    }
    A.setFromTriplets(t.begin(), t.end());

    return A;
}

float fQ(MatrixXf Q, VectorXf x) {
    return x.transpose() * Q * x;
}

MatrixXf g(MatrixXf P, float pr) {
    int n = P.outerSize();
    map<int, int> m;
    MatrixXf P_first(n, n);

    // with probability pr inserts i in the map with key i
    for (int i = 0; i < n; i++) {
        if (d_real_uniform(e_uniform_g) <= pr) {  // checks if the extracted numeber is equal or less than pr
            m.insert(pair<int, int>(i, i));       //if it is, inserts the number
        }
    }

    // shuffle keys
    shuffle_map(m);

    auto end = m.end();  // Used do determine whether the function find fails
    for (int i = 0; i < n; i++) {
        auto search = m.find(i);
        if (search != end) {  //if the element is found
            P_first.row(i) = P.row(m.at(i));
        } else {
            P_first.row(i) = P.row(i);
        }
    }

    //P_first is a permutation of P
    return P_first;
}

SparseMatrix<float> g_strong(MatrixXf Q, SparseMatrix<float> A, vector<int> &permutation, vector<int> old_permutation, double pr) {
    int n = Q.outerSize();
    map<int, int> m;
    SparseMatrix<float> theta(n, n);
    vector<int> inversed(n);
    vector<Triplet<float>> t;
    t.reserve(11 * n);  // A = (V, E) => |t| = |E| + |V|
    int r, c;
    double val;

    // with probability pr inserts i in the map with key i
    for (int i = 0; i < n; i++) {
        if (d_real_uniform(e_uniform_g) <= pr) {  // checks if the extracted numeber is equal or less than pr
            m.insert(pair<int, int>(i, i));       // if it is, inserts the number
        }
    }

    shuffle_map(m);

    permutation = fill(m, old_permutation);  // Generates a vector of permuted + non permuted indexes
    inversed = inverse(permutation);         // Inversed is used to know which pair row column is to be used to calculate the product Q ○ A

    for (int i = 0; i < A.outerSize(); i++) {
        for (SparseMatrix<float>::InnerIterator it(A, i); it; ++it) {
            r = it.row();
            c = it.col();
            val = Q(inversed[r], inversed[c]);
            t.push_back(Triplet<float>(r, c, val));
        }
    }
    theta.setFromTriplets(t.begin(), t.end());

    return theta;
}

void shuffle_map(map<int, int> &m) {
    vector<int> keys;  //Vector containing keys of m

    for (auto i : m) {
        keys.push_back(i.first);  //add m's keys to keys vector
    }

    // Shuffle
    //shuffle(keys.begin(), keys.end(), e_uniform_vector);
    //random_shuffle(keys.begin(), keys.end());
    shuffle_vector(keys);

    vector<int>::iterator it = keys.begin();
    //substitute old keys with new ones (shuffled)
    for (auto &i : m) {
        int ts = i.second;
        i.second = m[*it];
        m[*it] = ts;
        it++;
    }
}

// Not used
void shuffle_vector(vector<int> &v) {  // Fisher and Yates' algorithm
    int n = v.size();
    int j;

    for (int i = n - 1; i >= 0; i--) {
        j = d_int_uniform(e_uniform_vector) * i / 2047;
        swap(v[i], v[j]);
    }
}

vector<int> fill(map<int, int> m, vector<int> permutation) {
    int n = permutation.size();
    vector<int> filled(n);
    auto end = m.end();

    for (int i = 0; i < n; i++) {
        auto search = m.find(i);
        if (search != end) {                   // if index i is contained
            filled[i] = permutation[m.at(i)];  // filled[i] is equal to the permutation vector at position permutation[m[i]]
        } else {                               // if is not contained
            filled[i] = permutation[i];        // the index is the same as the old permutation
        }
    }

    return filled;
}

vector<int> inverse(vector<int> permutation) {
    int n = permutation.size();
    vector<int> inverted(n);
    for (int i = 0; i < n; i++) {
        inverted[permutation[i]] = i;
    }
    return inverted;
}

#ifndef SIMULATION
VectorXf send_to_annealer(SparseMatrix<float> theta, int *fd) {
    int n = theta.outerSize();
    char r[100];
    char c[100];
    char val[100];
    char ret[3];
    VectorXf z(n);

    memset(r, '\0', sizeof(char) * 100);
    memset(c, '\0', sizeof(char) * 100);
    memset(val, '\0', sizeof(char) * 100);  // to be adjusted, could be too small

    for (int i = 0; i < theta.outerSize(); i++) {
        for (SparseMatrix<float>::InnerIterator it(theta, i); it; ++it) {
            sprintf(r, "%ld", it.row());
            sprintf(c, "%ld", it.col());
            sprintf(val, "%lf", it.value());

            write(fd[WRITE], r, 100);  // 100 is very generous
            write(fd[WRITE], c, 100);
            write(fd[WRITE], val, 100);
        }
    }
    sprintf(val, "%s", "#\0");
    write(fd[WRITE], val, 100);

    for (int i = 0; i < n; i++) {
        read(fd[READ + 2], ret, 2);
        z(i) = atof(ret);
    }
    return z;
}
#endif

void h(VectorXf &z, double pr) {
    int n = z.size();

    for (int i = 0; i < n; i++) {
        if (d_real_uniform(e_uniform_h) <= pr) z(i) = -z(i);
    }
}

double min(double lambda0, int i, int e) {
    double lambda_first = lambda0 / (2 + i - e);

    if (lambda0 < lambda_first) return lambda0;
    return lambda_first;
}

VectorXf min_energy(SparseMatrix<float> theta) {
    int n = theta.outerSize();
    unsigned long long N = pow(2, n);  // Overflow with n > 64, not a problem since is a simulation
    VectorXf x_min(n);
    VectorXf x(n);
    double min;
    double e;
    for (int i = 0; i < n; i++) x(i) = -1;

    min = E(theta, x);
    x_min = x;
    unsigned long long i = 1;
    do {
        increment(x);
        e = E(theta, x);
        if (e < min) {
            x_min = x;
            min = e;
        }
        i++;
    } while (i < N);
    return x_min;
}

double E(SparseMatrix<float> theta, VectorXf x) {
    double e = 0;
    int r, c;
    for (int i = 0; i < theta.outerSize(); i++) {
        for (SparseMatrix<float>::InnerIterator it(theta, i); it; ++it) {
            r = it.row();
            c = it.col();
            if (r == c)
                e += it.value() * x(r);
            else
                e += it.value() * x(r) * x(c);
        }
    }
    return e;
}

void increment(VectorXf &v) {  // O(1) per l'analisi ammortizzata
    int n = v.size();
    int i = 0;
    while (i < n && v(i) == 1) {
        v(i) = -1;
        i++;
    }
    if (i < n) v(i) = 1;
}

VectorXf map_back(VectorXf z, vector<int> perm) {
    int n = perm.size();
    vector<int> inverted = inverse(perm);
    VectorXf z_ret(n);

    for (int i = 0; i < n; i++) {
        z_ret(i) = z(inverted[i]);
    }

    return z_ret;
}

double simulated_annealing(double f_first, double f_star, double p) {
    double T = -1 / log(p);
    return exp(-(f_first - f_star) / T);
}

bool comp_vectors(VectorXf z1, VectorXf z2) {
    for (int i = 0; i < z1.size(); i++) {
        if (abs(z1(i) - z2(i)) > __FLT_EPSILON__) return false;
    }
    return true;
}

#ifdef SIMULATION
SparseMatrix<float> gen_P(vector<int> perm) {
    long unsigned n = perm.size();
    SparseMatrix<float> P(n, n);
    vector<Triplet<float>> t;
    t.reserve(n);

    for (long unsigned i = 0; i < n; i++) {
        t.push_back(Triplet<float>(i, perm[i], 1));
    }
    P.setFromTriplets(t.begin(), t.end());

    return P;
}

double compute_Q(MatrixXf Q) {
    int n = Q.outerSize();
    VectorXf x_min(n);
    VectorXf x(n);
    unsigned long long N = pow(2, n);
    double min, e;
    for (int i = 0; i < n; i++) x(i) = -1;

    min = fQ(Q, x);
    x_min = x;
    unsigned long long i = 1;
    do {
        increment(x);
        e = fQ(Q, x);
        if (e < min) {
            x_min = x;
            min = e;
        }
        i++;
    } while (i < N);

    return min;
}

void log(VectorXf z_star, double f_star, double min, double f_gold, double lambda, double p, int e, int d, bool perturbed, bool simul_ann, int i) {
    cout << "---Current status at " << i << "th iteration---" << endl
         << "Best so far: f*=" << f_star << "\tz*=" << z_star.transpose() << endl
         << "To reach: min=" << min << "\tf_gold=" << f_gold << endl
         << "λ=" << lambda << "\tp=" << p << "\te=" << e << "\td=" << d;
    if (perturbed) cout << "\tperturbed";
    if (simul_ann) cout << "\tsimulated annealing";
}
#else
void log(VectorXf z_star, double f_star, double f_gold, double lambda, double p, int e, int d, bool perturbed, bool simul_ann, int i) {
    cout << "---Current status at " << i << "th iteration---" << endl
         << "Best so far: f*=" << f_star << "\tz*=" << z_star.transpose() << endl
         << "f_gold=" << f_gold << endl
         << "λ=" << lambda << "\tp=" << p << "\te=" << e << "\td=" << d;
    if (perturbed) cout << "\tperturbed";
    if (simul_ann) cout << "\tsimulated annealing";
}
#endif
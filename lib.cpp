#include "lib.h"

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
    shuffle(m);

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

void shuffle(map<int, int> &m) {
    vector<int> keys;  //Vector containing keys of m

    for (auto i : m) {
        keys.push_back(i.first);  //add m's keys to keys vector
    }

// Shuffle
    random_shuffle(keys.begin(), keys.end(), myrandom);

    vector<int>::iterator it = keys.begin();
    //substitute old keys with new ones (shuffled)
    for (auto &i : m) {
        int ts = i.second;
        i.second = m[*it];
        m[*it] = ts;
        it++;
    }
}

int myrandom (int i) { return std::rand()%i;}

float fQ(MatrixXf Q, VectorXf x) {
    return x.transpose() * Q * x;
}

void h(VectorXf &z, float pr) {
    int n = z.size();

    for (int i = 0; i < n; i++) {
        if (d_real_uniform(e_uniform_h) <= pr) z(i) = -z(i);
    }
}

float simulated_annealing(float f_first, float f_star, float p) {
    float T = -1 / log(p);
    return exp(-(f_first - f_star) / T);
}

float min(float lambda0, int i, int e) {
    float lambda_first = lambda0 / (2 + i - e);

    if (lambda0 < lambda_first) return lambda0;
    return lambda_first;
}

void init_seeds() {
#ifndef DEBUG
    seed_g = rd() ^ rd();
    seed_h = rd() ^ rd();
    seed_ann = rd() ^ rd();
    seed_pert = rd() ^ rd();
    seed_shuffle = rd() ^ rd();
#else
    seed_g = 100000;
    seed_h = 200000;
    seed_ann = 300000;
    seed_pert = 400000;
    seed_shuffle = 500000;
#endif

    e_uniform_g.seed(seed_g);
    e_uniform_h.seed(seed_h);
    e_uniform_ann.seed(seed_ann);
    e_uniform_pert.seed(seed_pert);
    e_uniform_shuffle.seed(seed_shuffle);
}

bool comp_vectors(VectorXf z1, VectorXf z2) {
    for (int i = 0; i < z1.size(); i++) {
        if (abs(z1(i) - z2(i)) > __FLT_EPSILON__) return false;
    }
    return true;
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

VectorXf min_energy(SparseMatrix<float> theta) {
    int n = theta.outerSize();
    int N = pow(2, n);
    int i;
    VectorXf x_min(n);
    VectorXf x(n);
    float min;
    float e;
    for (i = 0; i < n; i++) x(i) = -1;

    min = E(theta, x);
    x_min = x;
    i = 1;
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

void increment(VectorXf &v) {  // O(1) per l'analisi ammortizzata
    int n = v.size();
    int i = 0;
    while (i < n && v(i) == 1) {
        v(i) = -1;
        i++;
    }
    if (i < n) v(i) = 1;
}

float E(SparseMatrix<float> theta, VectorXf x) {
    float e = 0;
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

void log(VectorXf z_star, float f_star, float lambda, float p, int e, int d, bool perturbed, bool simul_ann, int i) {
    cout << "---Current status at " << i << "th iteration---" << endl
         << "Best so far: f*=" << f_star << "\tz*=" << z_star.transpose() << endl
         << "λ=" << lambda << "\tp=" << p << "\te=" << e << "\td=" << d;
    if (perturbed) cout << "\tperturbed";
    if (simul_ann) cout << "\tsimulated annealing";
    cout << endl
         << endl;
}
#include "lib.h"

void handle_sigint(int sig) {
    char send[100];
    memset(send, '\0', 100);
    sprintf(send, "%s", "END");

    write(fd[WRITE], send, 100);
    wait(NULL);
}

void init_child() {
    char first[8];
    char second[12];
    memset(first, '\0', sizeof(char) * 8);
    memset(second, '\0', sizeof(char) * 12);
    sprintf(first, "python3");
    sprintf(second, "./solver.py");
    char *args[] = {first, second, NULL};

    dup2(fd[READ], STDIN_FILENO);
    dup2(fd[WRITE + 2], STDOUT_FILENO);

    close(fd[READ]);
    close(fd[WRITE]);
    close(fd[READ + 2]);
    close(fd[WRITE + 2]);

    if (execvp(args[0], args) == -1) {
        cout << "[EXECVP ERROR - CLOSING]" << endl;
        exit(3);
    }
}

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
#ifdef SIMULATION
    int sim = 1;
#else
    int sim = 0;
#endif
    char n_nodes[6];
    char simulation[2];
    char i[5];
    char j[5];
    int r;
    int c;
    bool end = false;
    SparseMatrix<float> A(n, n);
    vector<Triplet<float>> t;

    memset(n_nodes, '\0', sizeof(char) * 6);
    memset(simulation, '\0', sizeof(char) * 2);
    memset(i, '\0', 5);
    memset(j, '\0', 5);

    sprintf(n_nodes, "%d", n);
    sprintf(simulation, "%d", sim);

    write(fd[WRITE], n_nodes, 6);
    write(fd[WRITE], simulation, 2);

    do {
        read(fd[READ + 2], i, 4);
        read(fd[READ + 2], j, 4);
        if (strncmp(i, "####", 4) != 0 && strncmp(j, "####", 4) != 0) {
            r = atoi(i);
            c = atoi(j);
            t.push_back(Triplet<float>(r, c, 1.0f));
        } else
            end = true;
    } while (!end);

    A.setFromTriplets(t.begin(), t.end());
    return A;
}

float fQ(MatrixXf Q, VectorXf x) {
    return x.transpose() * Q * x;
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
VectorXf send_to_annealer(SparseMatrix<float> theta) {
    int n = theta.outerSize();
    int curr_size = 0;
    char msg[STR_MAX_LEN];
    char val[STR_MAX_LEN];
    char tmp[100];
    VectorXf z(n);

    memset(val, '\0', sizeof(char) * STR_MAX_LEN);
    memset(tmp, '\0', 100);
    strcpy(msg, "");

    for (int i = 0; i < n; i++) {
        for (SparseMatrix<float>::InnerIterator it(theta, i); it; ++it) {
            sprintf(tmp, "%ld", it.row());
            strcat(msg, tmp);
            strcat(msg, ",");
            memset(tmp, '\0', 100);

            sprintf(tmp, "%ld", it.col());
            strcat(msg, tmp);
            strcat(msg, ",");
            memset(tmp, '\0', 100);

            sprintf(tmp, "%lf", it.value());
            strcat(msg, tmp);
            strcat(msg, ",");
            memset(tmp, '\0', 100);

            if (strlen(msg) > 3500) {
                for (int j = strlen(msg) - 1; j < STR_MAX_LEN; j++) msg[j] = '\0';
                write(fd[WRITE], msg, STR_MAX_LEN);
                strcpy(msg, "");
            }
        }
    }
    if (strlen(msg) > 0) {
        for (int j = strlen(msg) - 1; j < STR_MAX_LEN; j++) msg[j] = '\0';
        write(fd[WRITE], msg, STR_MAX_LEN);
        cout << msg << endl;
        strcpy(msg, "");
    }
    sprintf(val, "%s", "#\0");
    write(fd[WRITE], val, STR_MAX_LEN);

    while (curr_size < n) {
        read(fd[READ + 2], msg, STR_MAX_LEN);
        csv_to_vector(z, msg, curr_size);
    }

    return z;
}

void csv_to_vector(VectorXf &z, char *msg, int &count) {
    char *token = strtok(msg, ",");
    while (token != NULL) {
        z[count] = atoi(token);
        count++;
        token = strtok(NULL, ",");
    }
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
         << "f*=" << f_star << "\tz*=" << z_star.transpose() << endl
         << "To reach: min=" << min << "\tf_gold=" << f_gold << endl
         << "λ=" << lambda << "\tp=" << p << "\te=" << e << "\td=" << d;
    if (perturbed) cout << "\tperturbed";
    if (simul_ann) cout << "\tsimulated annealing";
}
#else
void log(double f_star, double f_gold, double lambda, double p, int e, int d, bool perturbed, bool simul_ann, int i) {
    cout << "---Current status at " << i << "th iteration---" << endl
         << "f*=" << f_star << endl
         << "Best so far: f_gold=" << f_gold << endl
         << "λ=" << lambda << "\tp=" << p << "\te=" << e << "\td=" << d;
    if (perturbed) cout << "\tperturbed";
    if (simul_ann) cout << "\tsimulated annealing";
}
#endif
#include "../lib/lib.h"

mt19937_64 e_uniform_g;
mt19937_64 e_uniform_h;
mt19937_64 e_uniform_shuffle;
mt19937_64 e_uniform_ann;
mt19937_64 e_uniform_pert;
mt19937_64 e_uniform_vector;
uniform_real_distribution<double> d_real_uniform;
uniform_int_distribution<unsigned long long> d_int_uniform(0, 5436);
pid_t child_pid;
int fd[4];

VectorXf solve(MatrixXf Q, int imax, int mode, int k, bool logs) {
    //Init
    int n = Q.outerSize();

    if (mode != BINARY && mode != SPIN) {
        cout << "[Warning] mode must be 0 or 1" << endl;
        exit(1);
    }

#ifdef SIMULATION
    if (n > 64) {
        cout << "[Warning] problem's dimensions are too big for a classical computation" << endl;
        exit(1);
    }
#endif

    /*---------------------------
        fd[READ] child read
        fd[WRITE] father write
        fd[READ + 2] father read
        fd[WRTIE + 2] child write
    ---------------------------*/

    if (pipe(fd) != 0) {
        cout << "[PIPE1 ERROR - CLOSING]" << endl;
    }
    if (pipe(fd + 2) != 0) {
        cout << "[PIPE2 ERROR - CLOSING]" << endl;
    }

    child_pid = fork();

    if (child_pid == 0) {
        init_child(mode, k);
    } else if (child_pid == -1) {
        cout << "[FORK ERROR - CLOSING]" << endl;
        exit(4);
    }

    init_seeds();
    unordered_map<int, int> nodes;
    SparseMatrix<float> edges;
    get_topology(nodes, edges, n);  // Pegasus topology

    //Input
    double pmin = 0.2f;     // minimum probability 0 < pδ < 0.5 of permutation modification
    double eta = 0.01f;     // probability decreasing rate η > 0
    double q = 0.1f;        // candidate perturbation probability q > 0
    double lambda0 = 1.0f;  // initial balancing factor λ0 > 0
    int N = 20;             // Decreasing time

    //Termination Parameters
    int Nmax = 100;  // Max number of solution equal to the best one + solution worse than the best one
    int dmin = 70;   // Number of solution that are worse than the best beyond which the best solution is not valid anymore

    MatrixXf In(n, n);  //Identity matrix
    In.setIdentity();

#ifdef SIMULATION
    cout << "Q" << endl
         << Q << endl
         << endl;
    cout << "edges" << endl
         << edges << endl
         << endl;
#endif

    //Algorithm
    MatrixXf Q_first(n, n);
    SparseMatrix<float> theta1(n, n), theta2(n, n), theta_first(n, n);
    VectorXf z_star(n), z_first(n), z1(n), z2(n), z_gold(n);
    MatrixXf z_diag(n, n);
    MatrixXf S(n, n);  //Tabu Matrix
    vector<int> perm(n), perm_star(n), perm1(n), perm2(n);
    double f1, f2, f_star, f_first, f_gold;
    double p = 1;             // probability of an element to be considered for shuffling
    int e = 0;                // Number of solutions that equal the best
    int d = 0;                // Number of solutions that are sub optimal
    double lambda = lambda0;  // Tabu weight
    bool perturbed;           // True when h perturbs the candidate
    bool simul_ann;
    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();
    string filename = "../out/tmp-" + to_string(time(0)) + ".txt";
    string problem = "../out/prob-" + to_string(time(0)) + ".txt";
    ofstream out_file;
    
    out_file.open(problem);
    out_file << Q;
    out_file.close();


    // Initialization of perm vectors like an identity matrix of order n
    for (int i = 0; i < n; i++) {
        perm[i] = i;
        perm1[i] = i;
        perm2[i] = i;
    }

    // Hadamard product between a permuted Q and edges
    theta1 = g_strong(Q, nodes, edges, perm1, perm1, p);
    theta2 = g_strong(Q, nodes, edges, perm2, perm2, p);

#ifdef SIMULATION
    cout << "Computing min of Q" << endl;
    double minimum = compute_Q(Q, mode);  // Global minimum of Q, only for simulation pourposes

    cout << "Computing z1" << endl;
    z1 = map_back(min_energy(theta1, mode), perm1);

    cout << "Computing z2" << endl;
    z2 = map_back(min_energy(theta2, mode), perm2);
#else
    cout << "Computing z1" << endl;
    z1 = map_back(send_to_annealer(theta1, n), perm1);

    cout << "Computing z2" << endl;
    z2 = map_back(send_to_annealer(theta2, n), perm2);
#endif

    f1 = fQ(Q, z1);
    f2 = fQ(Q, z2);

    if (f1 < f2) {    // f1 is better
        z_star = z1;  // Best
        f_star = f1;
        perm_star = perm1;
        z_first = z2;  // Worst

        z_gold = z1;
        f_gold = f1;
    } else {  // f2 is better
        z_star = z2;
        f_star = f2;
        perm_star = perm2;
        z_first = z1;

        z_gold = z2;
        f_gold = f2;
    }

    // f1 and f2 are floats -> float comparison
    if (abs(f1 - f2) > __FLT_EPSILON__) {  // if(f1 != f2)
        z_diag = z_first.asDiagonal();     // Matrix where the diagonal is made by elemnt of z_first
        S = kroneckerProduct(z_first, z_first.transpose()) - In + z_diag;
    } else {
        S = MatrixXf::Zero(n, n);
    }

    int i = 1;
    do {
        if (i == 1) cout << "Start computing the solution" << endl
                         << endl;
        start = chrono::steady_clock::now();
        perturbed = false;
        simul_ann = false;

        Q_first = Q + lambda * S;

        if (!(i % N))
            p = p - (p - pmin) * eta;  // 0 mod N va considerato come 0?

        theta_first = g_strong(Q_first, nodes, edges, perm, perm_star, p);
#ifdef SIMULATION
        z_first = map_back(min_energy(theta_first, mode), perm);
#else
        z_first = map_back(send_to_annealer(theta_first, n), perm);
#endif
        if (d_real_uniform(e_uniform_pert) <= q) {
            h(z_first, p, mode);  // possibly perturb the candidate
            perturbed = true;
        }

        if (!comp_vectors(z_first, z_star)) {
            f_first = fQ(Q, z_first);
            if (f_first < f_star) {     // f_first is better
                swap(z_first, z_star);  // z_first is better
                f_star = f_first;
                perm_star = perm;
                e = 0;
                d = 0;
                z_diag = z_first.asDiagonal();
                S = S + kroneckerProduct(z_first, z_first.transpose()) - In + z_diag;
            } else {
                d++;
                if (d_real_uniform(e_uniform_ann) <= pow(p - pmin, f_first - f_star)) {
                    swap(z_first, z_star);
                    f_star = f_first;
                    perm_star = perm;
                    e = 0;
                    simul_ann = true;
                }
            }

            // Best solutions yet
            if (f_star < f_gold) {
                z_gold = z_star;
                f_gold = f_star;
            }

            lambda = min(lambda0, i - 1, e);
        } else {
            e++;
        }

#ifdef SIMULATION
        if (logs) log(Q, z_star, f_star, minimum, f_gold, lambda, p, e, d, perturbed, simul_ann, i, filename);
#else
        if (logs) log(z_star, f_star, f_gold, lambda, p, e, d, perturbed, simul_ann, i, filename);
#endif
        end = chrono::steady_clock::now();
        chrono::duration<double> diff = end - start;
        if (logs) {
            cout << endl
                 << diff.count() << "s" << endl
                 << endl;
        }
        i++;
    } while (i <= imax && (e + d < Nmax || d >= dmin));

#ifndef SIMULATION
    close_child();
#endif

    close(fd[READ]);
    close(fd[WRITE]);
    close(fd[READ + 2]);
    close(fd[WRITE + 2]);
    
    remove(filename.c_str());

    if (logs) {
        printf("pmin:%f\teta:%f\tq:%f\tlambda0:%f\tN:%d\n", pmin, eta, q, lambda0, N);
        printf("k:%d\n", k);
        printf("imax:%d, Nmax:%d, dmin:%d\n", imax, Nmax, dmin);
        printf("e:%d\td:%d\ti:%d\n", e, d, i - 1);

        cout << endl
             << "f_gold: " << f_gold << endl
             << endl;
    }

    return z_gold;
}

void init_child(int mode, int k) {
    char first[20];
    char second[16];
    char third[3];
    char fourth[30];
    memset(first, '\0', sizeof(char) * 20);
    memset(second, '\0', sizeof(char) * 16);
    memset(third, '\0', sizeof(char) * 3);
    memset(fourth, '\0', sizeof(char) * 30);
    sprintf(first, "python3");
    sprintf(second, "../solver.py");
    sprintf(third, "%d", mode);
    sprintf(fourth, "%d", k);

    char *args[] = {first, second, third, fourth, NULL};

    dup2(fd[READ], STDIN_FILENO);        // Change child's stdin
    dup2(fd[WRITE + 2], STDOUT_FILENO);  // Change child's stdout

    // Close pipes
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

void get_topology(unordered_map<int, int> &nodes, SparseMatrix<float> &edges, int n) {
// sim tells python if it's a simulation or not
#ifdef SIMULATION
    int sim = 1;
#else
    int sim = 0;
#endif
    int k;
    char n_nodes[10];
    char simulation[2];
    int r;
    int c;
    int len_n = (to_string(n)).length() + 5;
    char i[len_n + 1];
    char j[len_n + 1];
    bool end = false;
    vector<Triplet<float>> t;
    unordered_map<int, int> swapped;

    memset(n_nodes, '\0', sizeof(char) * 10);
    memset(simulation, '\0', sizeof(char) * 2);
    memset(i, '\0', len_n + 1);
    memset(j, '\0', len_n + 1);

    sprintf(n_nodes, "%d", n);
    sprintf(simulation, "%d", sim);

    write(fd[WRITE], n_nodes, 10);    // Send number of nodes in the problem
    write(fd[WRITE], simulation, 2);  // Send if it's a simulation or not

    for (k = 0; k < n; k++) {
        read(fd[READ + 2], i, len_n);
        nodes.insert(pair<int, int>(k, atoi(i)));
    }
    edges.resize(nodes.at(k - 1) + 1, nodes.at(k - 1) + 1);

    for (auto it = nodes.begin(); it != nodes.end(); it++) {
        swapped.insert(pair<int, int>(it->second, it->first));
    }
    nodes = swapped;

    do {
        read(fd[READ + 2], i, len_n);  // Read i index
        read(fd[READ + 2], j, len_n);  // Read j index

        if (strncmp(i, "#", 1) != 0 && strncmp(j, "#", 1) != 0) {
            r = atoi(i);  // Set r as i index if i is not equal to "####"
            c = atoi(j);  // Set c as j index if j is not equal to "####"
            t.push_back(Triplet<float>(r, c, 1.0f));
        } else
            end = true;  // if "####" is recived then the matrix is totally sended
    } while (!end);

    edges.setFromTriplets(t.begin(), t.end());
}

float fQ(MatrixXf Q, VectorXf x) {
    return x.transpose() * Q * x;
}

SparseMatrix<float> g_strong(const MatrixXf &Q, const unordered_map<int, int> &nodes, const SparseMatrix<float> &edges, vector<int> &permutation, const vector<int> &old_permutation, double pr) {
    int n = Q.outerSize();
    map<int, int> m;
    SparseMatrix<float> theta(edges.outerSize(), edges.outerSize());
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

    for (int i = 0; i < edges.outerSize(); i++) {
        for (SparseMatrix<float>::InnerIterator it(edges, i); it; ++it) {
            r = it.row();
            c = it.col();
            val = Q(inversed[nodes.at(r)], inversed[nodes.at(c)]);
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
    shuffle_vector(keys);

    vector<int>::iterator it = keys.begin();
    // substitute old keys with new ones (shuffled)
    for (auto &i : m) {
        int ts = i.second;
        i.second = m[*it];
        m[*it] = ts;
        it++;
    }
}

void shuffle_vector(vector<int> &v) {  // Fisher and Yates' algorithm
    int n = v.size();
    int j;

    for (int i = n - 1; i >= 0; i--) {
        j = d_int_uniform(e_uniform_vector) * i / 5435;
        swap(v[i], v[j]);
    }
}

vector<int> fill(const map<int, int> &m, const vector<int> &permutation) {
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

vector<int> inverse(const vector<int> &permutation) {
    int n = permutation.size();
    vector<int> inverted(n);
    for (int i = 0; i < n; i++) {
        inverted[permutation[i]] = i;
    }
    return inverted;
}

#ifndef SIMULATION
VectorXf send_to_annealer(const SparseMatrix<float> &theta, int n) {
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

void h(VectorXf &z, double pr, int mode) {
    int n = z.size();

    for (int i = 0; i < n; i++) {
        if (d_real_uniform(e_uniform_h) <= pr) mode == SPIN ? z(i) = -z(i) : z(i) = (int)(z(i) + 1) % 2;
    }
}

double min(double lambda0, int i, int e) {
    double lambda_first = lambda0 / (2 + i - e);

    if (lambda0 < lambda_first) return lambda0;
    return lambda_first;
}

VectorXf min_energy(const SparseMatrix<float> &theta, int mode) {
    int n = theta.outerSize();
    unsigned long long N = pow(2, n);  // Overflow with n > 64, not a problem since is a simulation
    VectorXf x_min(n);
    VectorXf x(n);
    double min;
    double e;
    for (int i = 0; i < n; i++) x(i) = mode;

    min = E(theta, x);
    x_min = x;
    unsigned long long i = 1;
    do {
        increment(x, mode);
        e = E(theta, x);
        if (e < min) {
            x_min = x;
            min = e;
        }
        i++;
    } while (i < N);
    return x_min;
}

double E(const SparseMatrix<float> &theta, VectorXf x) {
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

void increment(VectorXf &v, int mode) {  // O(1) per l'analisi ammortizzata
    int n = v.size();
    int i = 0;
    while (i < n && v(i) == 1) {
        v(i) = mode;
        i++;
    }
    if (i < n) v(i) = 1;
}

VectorXf map_back(const VectorXf &z, const vector<int> &perm) {
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

bool comp_vectors(const VectorXf &z1, const VectorXf &z2) {
    for (int i = 0; i < z1.size(); i++) {
        if (abs(z1(i) - z2(i)) > __FLT_EPSILON__) return false;
    }
    return true;
}

void close_child() {
    char send[100];
    memset(send, '\0', 100);
    sprintf(send, "%s", "END");

    write(fd[WRITE], send, 100);
    wait(NULL);
}

#ifdef SIMULATION
double compute_Q(const MatrixXf &Q, int mode) {
    int n = Q.outerSize();
    VectorXf x_min(n);
    VectorXf x(n);
    unsigned long long N = pow(2, n);
    double min, e;
    for (int i = 0; i < n; i++) x(i) = mode;

    min = fQ(Q, x);
    x_min = x;
    unsigned long long i = 1;
    do {
        increment(x, mode);
        e = fQ(Q, x);
        if (e <= min) {
            x_min = x;
            min = e;
        }
        i++;
    } while (i < N);
    return min;
}

void log(const MatrixXf &Q, const VectorXf &z_star, double f_star, double min, double f_gold, double lambda, double p, int e, int d, bool perturbed, bool simul_ann, int i, string filename) {
    ofstream out_file;
    stringstream ss;
    ss << z_star.transpose();
    string out = "---Current status at " + to_string(i) + "th iteration---\n" + "f*=" + to_string(f_star) + "\tz*=" + ss.str() + "\n" + "To reach: min=" + to_string(min) + "\tf_gold=" + to_string(f_gold) + "\n" + "λ=" + to_string(lambda) + "\tp=" + to_string(p) + "\te=" + to_string(e) + "\td=" + to_string(d);

    if (perturbed) out += "\tperturbed";
    if (simul_ann) out += "\tsimulated annealing";

    cout << out;

    ss.str("");
    ss << Q;
    out += "\n" + ss.str();
    out_file.open(filename);
    out_file << out;
    out_file.close();
}
#else
void log(const VectorXf &z_star, double f_star, double f_gold, double lambda, double p, int e, int d, bool perturbed, bool simul_ann, int i, string filename) {
    ofstream out_file;
    stringstream ss;
    string out = "---Current status at " + to_string(i) + "th iteration---\n" + "f*=" + to_string(f_star) + "\tf_gold=" + to_string(f_gold) + "\n" + "λ=" + to_string(lambda) + "\tp=" + to_string(p) + "\te=" + to_string(e) + "\td=" + to_string(d);

    if (perturbed) out += "\tperturbed";
    if (simul_ann) out += "\tsimulated annealing";

    cout << out;

    ss << z_star.transpose();
    out += "\n" + ss.str();
    out_file.open(filename);
    out_file << out;
    out_file.close();
}
#endif

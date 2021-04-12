#include "../lib/lib.h"

mt19937_64 e_uniform_g;
mt19937_64 e_uniform_h;
mt19937_64 e_uniform_shuffle;
mt19937_64 e_uniform_ann;
mt19937_64 e_uniform_pert;
mt19937_64 e_uniform_vector;
uniform_real_distribution<double> d_real_uniform;
pid_t child_pid;
int fd[4];

VectorXd solve(MatrixXd Q, int imax, int k, short logs, string filename) {
    //Init
    int n = Q.outerSize();  // Get number of vaiables

#ifdef SIMULATION
    if (n > 64) {
        cout << "[Warning] problem's dimensions are too big for a classical computation" << endl;
        exit(1);
    }
#endif

    //logs
    bool log_console;
    bool log_file;
    logs % 2 == 0 ? log_console = false : log_console = true;
    (logs >> 1) % 2 == 0 ? log_file = false : log_file = true;

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
        init_child(k);
    } else if (child_pid == -1) {
        cout << "[FORK ERROR - CLOSING]" << endl;
        exit(4);
    }

    init_seeds(filename, log_file);
    unordered_map<int, int> nodes;  // map from topology nodes to ideal ones
    SparseMatrix<double> edges;     // SparseMatrix of edges
    get_topology(nodes, edges, n);  // Dwave topology

    //Input
    const double pmin = 0.1;       // minimum probability 0 < pδ < 0.5 of permutation modification
    const double eta = 0.01f;      // probability decreasing rate η > 0
    const double q = 0.2f;         // candidate perturbation probability q > 0
    const double lambda0 = 3 / 2;  // initial balancing factor λ0 > 0
    const int N = 10;              // Decreasing time

    //Termination Parameters
    const int Nmax = 100;  // Max number of solution equal to the best one + solution worse than the best one
    const int dmin = 70;   // Number of solution that are worse than the best beyond which the best solution is not valid anymore

    MatrixXd In(n, n);  //Identity matrix
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
    MatrixXd Q_first(n, n);
    SparseMatrix<double> theta1(n, n), theta2(n, n), theta_first(n, n);
    VectorXd z_star(n), z_first(n), z1(n), z2(n), z_gold(n);
    MatrixXd z_diag(n, n);
    MatrixXd S(n, n);  //Tabu Matrix
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
    string problem = "../out/prob-" + filename + ".txt";  // Temporary problem file
    string file_log = "../out/tmp-" + filename + ".txt";
    ofstream out_file;

    if (log_file) {
        // With large n could take some time (approx from 30s to 1 min)
        out_file.open(problem);
        out_file << Q;
        out_file.close();
    }

    // Initialization of perm vectors like an identity matrix of order n according to thesis
    for (int i = 0; i < n; i++) {
        perm[i] = i;
        perm1[i] = i;
        perm2[i] = i;
    }

    // Construction of matrix theta considering DWave's topology
    theta1 = g_strong(Q, nodes, edges, perm1, perm1, p);
    theta2 = g_strong(Q, nodes, edges, perm2, perm2, p);

#ifdef SIMULATION
    cout << "Computing min of Q" << endl;
    double minimum = compute_Q(Q);  // Global minimum of Q, only for simulation pourposes

    cout << "Computing z1" << endl;
    z1 = map_back(min_energy(theta1), perm1);

    cout << "Computing z2" << endl;
    z2 = map_back(min_energy(theta2), perm2);
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
        S = MatrixXd::Zero(n, n);
    }

    int i = 1;
    /*int count_acc = 0;
    string acc_prob = "../out/acc-" + filename + ".txt";
    out_file.open(acc_prob);
    stringstream ss;
    ss << Q;
    out_file << ss.str() << endl;
    ss.str("");*/
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
        z_first = map_back(min_energy(theta_first), perm);
#else
        z_first = map_back(send_to_annealer(theta_first, n), perm);
#endif
        if (d_real_uniform(e_uniform_pert) <= q) {
            h(z_first, p);  // possibly perturb the candidate
            perturbed = true;
        }

        //TMP
        /*if (is_acceptable(z_first)) {
            out_file << "ACC:" << endl;
            count_acc++;
            out_file << g_strong(Q, nodes, edges, perm, perm, 0);
            out_file << endl;
            for (auto it : perm) out_file << it << " ";
            out_file << endl;
        } else if (i == (int)imax / 4) {
            out_file << "NA:" << endl;
            out_file << g_strong(Q, nodes, edges, perm, perm, 0);
            out_file << endl;
            for (auto it : perm) out_file << it << " ";
            out_file << endl;
        }*/

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
        log(Q, z_star, f_star, minimum, f_gold, lambda, p, e, d, perturbed, simul_ann, i, imax, filename, logs);
#else
        log(z_star, f_star, f_gold, lambda, p, e, d, perturbed, simul_ann, i, imax, filename, logs);
#endif
        end = chrono::steady_clock::now();
        chrono::duration<double> diff = end - start;
        if (log_console) {
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

    if (log_file) {
        remove(file_log.c_str());
        remove(problem.c_str());
    }
    if (log_console) {
        printf("pmin:%f\teta:%f\tq:%f\tlambda0:%f\tN:%d\n", pmin, eta, q, lambda0, N);
        printf("k:%d\n", k);
        printf("imax:%d, Nmax:%d, dmin:%d\n", imax, Nmax, dmin);
        printf("e:%d\td:%d\ti:%d\n", e, d, i - 1);

        cout << endl
             << "f_gold: " << f_gold << endl;
    }
#ifdef SIMULATION
    if (f_gold == minimum) {
        cout << "minimum reached" << endl;
    } else {
        cout << "minimum not reached" << endl;
    }
#endif

    //cout << "Count Acc: " << count_acc << endl;

    return z_gold;
}

void init_child(int k) {
    char first[20];
    char second[16];
    char third[30];
    memset(first, '\0', sizeof(char) * 20);
    memset(second, '\0', sizeof(char) * 16);
    memset(third, '\0', sizeof(char) * 30);
    sprintf(first, "python3");
    sprintf(second, "../solver.py");
    sprintf(third, "%d", k);

    char *args[] = {first, second, third, NULL};

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

void init_seeds(string filename, bool log) {
    random_device rd;
    seed_seq seed_g{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    seed_seq seed_h{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    seed_seq seed_shuffle{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    seed_seq seed_ann{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    seed_seq seed_pert{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    seed_seq seed_vector{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    if (log) {
        ofstream of;
        of.open("../out/seed-" + filename + ".txt");
        of << "g=";
        seed_g.param(ostream_iterator<int>(of, " "));
        of << endl
           << "h=";
        seed_h.param(ostream_iterator<int>(of, " "));
        of << endl
           << "shuffle=";
        seed_shuffle.param(ostream_iterator<int>(of, " "));
        of << endl
           << "ann=";
        seed_ann.param(ostream_iterator<int>(of, " "));
        of << endl
           << "pert=";
        seed_pert.param(ostream_iterator<int>(of, " "));
        of << endl
           << "vector=";
        seed_pert.param(ostream_iterator<int>(of, " "));
        of.close();
    }

    e_uniform_g.seed(seed_g);
    e_uniform_h.seed(seed_h);
    e_uniform_ann.seed(seed_ann);
    e_uniform_pert.seed(seed_pert);
    e_uniform_shuffle.seed(seed_shuffle);
    e_uniform_vector.seed(seed_vector);
}

void get_topology(unordered_map<int, int> &nodes, SparseMatrix<double> &edges, int n) {
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
    vector<Triplet<double>> t;
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
            t.push_back(Triplet<double>(r, c, 1.0f));
        } else
            end = true;  // if "####" is recived then the matrix is totally sended
    } while (!end);

    edges.setFromTriplets(t.begin(), t.end());
}

double fQ(MatrixXd Q, VectorXd x) {
    return x.transpose() * Q * x;
}

SparseMatrix<double> g_strong(const MatrixXd &Q, const unordered_map<int, int> &nodes, const SparseMatrix<double> &edges, vector<int> &permutation, const vector<int> &old_permutation, double pr) {
    int n = Q.outerSize();
    map<int, int> m;
    SparseMatrix<double> theta(edges.outerSize(), edges.outerSize());
    vector<int> inversed(n);
    vector<Triplet<double>> t;
    t.reserve(16 * n);  // A = (V, E) => |t| = |E| + |V|
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
        for (SparseMatrix<double>::InnerIterator it(edges, i); it; ++it) {
            r = it.row();
            c = it.col();
            val = Q(inversed[nodes.at(r)], inversed[nodes.at(c)]);
            t.push_back(Triplet<double>(r, c, val));
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
    uniform_int_distribution<int> d_int_uniform(0, n - 1);

    for (int i = n - 1; i >= 0; i--) {
        j = d_int_uniform(e_uniform_vector);
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
VectorXd send_to_annealer(const SparseMatrix<double> &theta, int n) {
    char r[100];
    char c[100];
    char val[100];
    char ret[3];
    VectorXd z(n);

    memset(r, '\0', sizeof(char) * 100);
    memset(c, '\0', sizeof(char) * 100);
    memset(val, '\0', sizeof(char) * 100);  // to be adjusted, could be too small

    for (int i = 0; i < theta.outerSize(); i++) {
        for (SparseMatrix<double>::InnerIterator it(theta, i); it; ++it) {
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

void h(VectorXd &z, double pr) {
    for (int i = 0; i < z.size(); i++) {
        if (d_real_uniform(e_uniform_h) <= pr) z(i) = (int)(z(i) + 1) % 2;
    }
}

double min(double lambda0, int i, int e) {
    double lambda_first = lambda0 / (2 + i - e);

    if (lambda0 < lambda_first) return lambda0;
    return lambda_first;
}

VectorXd min_energy(const SparseMatrix<double> &theta) {
    int n = theta.outerSize();
    unsigned long long N = pow(2, n);  // Overflow with n > 64, not a problem since is a simulation
    VectorXd x_min(n);
    VectorXd x(n);
    double min;
    double e;
    for (int i = 0; i < n; i++) x(i) = 0;

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

double E(const SparseMatrix<double> &theta, VectorXd x) {
    double e = 0;
    int r, c;
    for (int i = 0; i < theta.outerSize(); i++) {
        for (SparseMatrix<double>::InnerIterator it(theta, i); it; ++it) {
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

void increment(VectorXd &v) {  // O(1) per l'analisi ammortizzata
    int n = v.size();
    int i = 0;
    while (i < n && v(i) == 1) {
        v(i) = 0;
        i++;
    }
    if (i < n) v(i) = 1;
}

VectorXd map_back(const VectorXd &z, const vector<int> &perm) {
    int n = perm.size();
    vector<int> inverted = inverse(perm);
    VectorXd z_ret(n);

    for (int i = 0; i < n; i++) {
        z_ret(i) = z(inverted[i]);
    }

    return z_ret;
}

// Unused
double simulated_annealing(double f_first, double f_star, double p) {
    double T = -1 / log(p);
    return exp(-(f_first - f_star) / T);
}

bool comp_vectors(const VectorXd &z1, const VectorXd &z2) {
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
double compute_Q(const MatrixXd &Q) {
    int n = Q.outerSize();
    VectorXd x_min(n);
    VectorXd x(n);
    unsigned long long N = pow(2, n);
    double min, e;
    for (int i = 0; i < n; i++) x(i) = 0;

    min = fQ(Q, x);
    x_min = x;
    unsigned long long i = 1;
    do {
        increment(x);
        e = fQ(Q, x);
        if (e <= min) {
            x_min = x;
            min = e;
        }
        i++;
    } while (i < N);
    return min;
}

VectorXd compute_Q_vector(const MatrixXd &Q) {
    int n = Q.outerSize();
    VectorXd x_min(n);
    VectorXd x(n);
    unsigned long long N = pow(2, n);
    double min, e;
    for (int i = 0; i < n; i++) x(i) = 0;

    min = fQ(Q, x);
    x_min = x;
    unsigned long long i = 1;
    do {
        increment(x);
        e = fQ(Q, x);
        if (e <= min) {
            x_min = x;
            min = e;
        }
        i++;
    } while (i < N);
    return x_min;
}

void log(const MatrixXd &Q, const VectorXd &z_star, double f_star, double min, double f_gold, double lambda, double p, int e, int d, bool perturbed, bool simul_ann, int i, int imax, string filename, short logs) {
    ofstream out_file;
    stringstream ss;
    string out;
    string tmp_file;
    string file_log;

    if (logs) {
        ss << z_star.transpose();
        out = "---Current status at " + to_string(i) + "th iteration---\n" + "f*=" + to_string(f_star) + "\tz*=" + ss.str() + "\n" + "To reach: min=" + to_string(min) + "\tf_gold=" + to_string(f_gold) + "\n" + "λ=" + to_string(lambda) + "\tp=" + to_string(p) + "\te=" + to_string(e) + "\td=" + to_string(d);
        if (perturbed) out += "\tperturbed";
        if (simul_ann) out += "\tsimulated annealing";
    }

    if (logs % 2 != 0) {
        cout << out;
    }

    if ((logs >> 1) % 2 != 0) {
        tmp_file = "../out/tmp-" + filename + ".txt";
        file_log = "../out/" + filename + ".txt";
        ss.str("");
        ss << Q;
        out += "\n" + ss.str();
        out_file.open(tmp_file);
        out_file << out;
        out_file.close();
        int par = imax / 4;
        if (i % par == 0) {
            out_file.open(file_log, ios_base::app);
            out_file << out + "\n";
            out_file.close();
        }
    }
}
#else
void log(const VectorXd &z_star, double f_star, double f_gold, double lambda, double p, int e, int d, bool perturbed, bool simul_ann, int i, int imax, string filename, short logs) {
    ofstream out_file;
    stringstream ss;
    string out;
    string tmp_file;
    string file_log;

    if (logs) {
        out = "---Current status at " + to_string(i) + "th iteration---\n" + "f*=" + to_string(f_star) + "\tf_gold=" + to_string(f_gold) + "\n" + "λ=" + to_string(lambda) + "\tp=" + to_string(p) + "\te=" + to_string(e) + "\td=" + to_string(d);
        if (perturbed) out += "\tperturbed";
        if (simul_ann) out += "\tsimulated annealing";
    }

    if (logs % 2 != 0) {
        cout << out;
    }

    if ((logs >> 1) % 2 != 0) {
        tmp_file = "../out/tmp-" + filename + ".txt";
        file_log = "../out/" + filename + ".txt";
        ss << z_star.transpose();
        out += "\n" + ss.str();
        out_file.open(tmp_file);
        out_file << out;
        out_file.close();

        int par = imax / 4;
        if (i % par == 0) {
            out_file.open(file_log, ios_base::app);
            out_file << out + "\n";
            out_file.close();
        }
    }
}
#endif

bool is_acceptable(const VectorXd &x) {
    long long n = sqrt(x.size());
    vector<long long> count(n, 0);
    for (long long i = 0; i < n; i++) {
        for (long long j = 0; j < n; j++) {
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
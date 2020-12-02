#include "../lib/lib.h"

mt19937_64 e_uniform_g;
mt19937_64 e_uniform_h;
mt19937_64 e_uniform_shuffle;
mt19937_64 e_uniform_ann;
mt19937_64 e_uniform_pert;
mt19937_64 e_uniform_vector;
uniform_real_distribution<double> d_real_uniform;
uniform_int_distribution<unsigned long long> d_int_uniform(0, 2048);
pid_t child_pid;
int fd[4];

VectorXf solve(MatrixXf Q) {
    //Init
    int n = Q.outerSize();

    if (n % 8 != 0 || n == 0) {
        cout << "[Warning] n must be multiple of 8" << endl;
        exit(1);
    }

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
        init_child();

    } else if (child_pid == -1) {
        cout << "[FORK ERROR - CLOSING]" << endl;
        exit(4);
    }

    init_seeds();
    SparseMatrix<float> A = init_A(n);  //Chimera topology

    //Input
    double pmin = 0.2f;     // minimum probability 0 < pδ < 0.5 of permutation modification
    double eta = 0.01f;     // probability decreasing rate η > 0
    double q = 0.1f;        // candidate perturbation probability q > 0
    double lambda0 = 1.0f;  // initial balancing factor λ0 > 0
    int k = 1;              // number of annealer runs k ≥ 1
    int N = 20;             // Decreasing time

    //Termination Parameters
    int imax = 3000;  // Max number of iteration
    int Nmax = 50;    // Max number of solution equal to the best one + solution worse than the best one
    int dmin = 30;    // Number of solution that are worse than the best beyond which the best solution is not valid anymore

    MatrixXf In(n, n);  //Identity matrix
    In.setIdentity();

#ifdef SIMULATION
    cout << "Q" << endl
         << Q << endl
         << endl;
    cout << "A" << endl
         << A << endl
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

    // Initialization of perm vectors like an identity matrix of order n
    for (int i = 0; i < n; i++) {
        perm[i] = i;
        perm1[i] = i;
        perm2[i] = i;
    }

    // Hadamard product between a permuted Q and A
    theta1 = g_strong(Q, A, perm1, perm1, p);
    theta2 = g_strong(Q, A, perm2, perm2, p);

#ifdef SIMULATION
    double minimum = compute_Q(Q);  // Global minimum of Q, only for simulation pourposes

    z1 = map_back(min_energy(theta1), perm1);
    z2 = map_back(min_energy(theta2), perm2);
#else
    z1 = map_back(send_to_annealer(theta1), perm1);
    z2 = map_back(send_to_annealer(theta2), perm2);
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
        start = chrono::steady_clock::now();
        perturbed = false;
        simul_ann = false;

        Q_first = Q + lambda * S;

        if (!(i % N))
            p = p - (p - pmin) * eta;  // 0 mod N va considerato come 0?

        theta_first = g_strong(Q_first, A, perm, perm_star, p);
#ifdef SIMULATION
        z_first = map_back(min_energy(theta_first), perm);
#else
        z_first = map_back(send_to_annealer(theta_first), perm);
#endif
        if (d_real_uniform(e_uniform_pert) <= q) {
            h(z_first, p);  // possibly perturb the candidate
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
                if (d_real_uniform(e_uniform_ann) <= simulated_annealing(f_first, f_star, p)) {
                    swap(z_first, z_star);
                    f_star = f_first;
                    perm_star = perm;
                    e = 0;
                    simul_ann = true;
                }
            }

            // Best solutions yet
            if (f_first < f_gold) {
                z_gold = z_first;
                f_gold = f_first;
            }

            lambda = min(lambda0, i - 1, e);
        } else {
            e++;
        }

#ifdef SIMULATION
        log(z_star, f_star, minimum, f_gold, lambda, p, e, d, perturbed, simul_ann, i);
#else
        log(f_star, f_gold, lambda, p, e, d, perturbed, simul_ann, i);
#endif
        end = chrono::steady_clock::now();
        chrono::duration<double> diff = end - start;
        cout << endl
             << diff.count() << "s" << endl
             << endl;
        i++;
    } while (i <= imax && (e + d < Nmax || d >= dmin));

    close(fd[READ]);
    close(fd[WRITE]);
    close(fd[READ + 2]);
    close(fd[WRITE + 2]);

#ifndef SIMULATION
    kill(child_pid, SIGKILL);
#endif

    printf("pmin:%f\teta:%f\tq:%f\tlambda0:%f\tN:%d\n", pmin, eta, q, lambda0, N);
    printf("k:%d\n", k);
    printf("imax:%d, Nmax:%d, dmin:%d\n", imax, Nmax, dmin);
    printf("e:%d\td:%d\ti:%d\n", e, d, i - 1);

    cout << endl
         << "f_gold: " << f_gold << endl
         << endl;

    return z_gold;
}

void handle_sigint(int sig) {
    char send[100];
    memset(send, '\0', 100);
    sprintf(send, "%s", "END");

    write(fd[WRITE], send, 100);
    wait(NULL);
}

void init_child() {
    char first[20];
    char second[16];
    memset(first, '\0', sizeof(char) * 20);
    memset(second, '\0', sizeof(char) * 16);
    sprintf(first, "python3");
    sprintf(second, "../solver.py");
    char *args[] = {first, second, NULL};

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

SparseMatrix<float> init_A(int n) {
// sim tells python if it's a simulation or not
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

    write(fd[WRITE], n_nodes, 6);     // Send number of nodes in the problem
    write(fd[WRITE], simulation, 2);  // Send if it's a simulation or not

    do {
        read(fd[READ + 2], i, 4); // Read i index
        read(fd[READ + 2], j, 4); // Read j index
        if (strncmp(i, "####", 4) != 0 && strncmp(j, "####", 4) != 0) {
            r = atoi(i); // Set r as i index if i is not equal to "####"
            c = atoi(j); // Set c as j index if j is not equal to "####"
            t.push_back(Triplet<float>(r, c, 1.0f));
        } else
            end = true; // if "####" is recived then the matrix is totally sended
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
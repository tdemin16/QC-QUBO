#include "../lib/generators.h"

long long NPP::number_partitioning_problem(MatrixXd &Q, vector<int> &nums, int range) {
    if (range < 1) {
        cout << "Range must be >= 1" << endl;
        exit(1);
    }

    int n = nums.size();
    mt19937_64 e_uniform;
    uniform_int_distribution<long> d_int_uniform(1, range);
    random_device rd;
    seed_seq seed{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    long long c = 0;
    long long prod;

    Q.resize(n, n);
    e_uniform.seed(seed);

    for (int i = 0; i < n; i++) {
        c += nums[i] = d_int_uniform(e_uniform);
    }

    for (long long i = 0; i < n; i++) {
        for (long long j = i; j < n; j++) {
            if (i != j) {
                prod = nums[i] * nums[j];
                Q(i, j) = prod;
                Q(j, i) = prod;
            } else {
                Q(i, i) = nums[i] * (nums[i] - c);
            }
        }
    }

    return c;
}

double NPP::diff(const MatrixXd &Q, const VectorXd &x, long long c) {
    return sqrt(pow(c, 2) + 4 * fQ(Q, x));
}

void NPP::to_file(chrono::duration<double> difference, int it, double f, int n, int range, double diff, const VectorXd &x, const vector<int> nums, string filename) {
    long unsigned len = 0;
    ofstream of;
    of.open("../out/" + filename + ".txt");

    of << difference.count() << "s" << endl
       << "imax=" << it << endl
       << "f=" << f << endl
       << "DIM=" << n << endl
       << "RANGE=[1, " << range << "]" << endl
       << "DIFF=" << diff << endl
       << "Solution vector:" << endl
       << x.transpose() << endl
       << "Number sequence:" << endl;

    for (auto it : nums) {
        len++;
        of << to_string(it);
        if (len < nums.size()) of << ", ";
    }

    of.close();
}

double QAP::quadratic_assignment_problem(MatrixXd &Q, double &max_coeff, float lambda, string file) {
    int n;
    long long N;
    double pen;
    ifstream in;

    in.open(file);
    if (in.fail()) {
        cout << "File doesn't exist" << endl;
        exit(1);
    }

    in >> n;
    MatrixXd f(n, n);
    MatrixXd d(n, n);
    N = n * n;
    Q.resize(N, N);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            in >> f(i, j);
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            in >> d(i, j);
        }
    }

    Q = kroneckerProduct(f, d);
    max_coeff = Q.maxCoeff();

    pen = (lambda * max_coeff);

    long long k;
    for (long long i = 0; i < N; i++) {
        k = (i / n) * n;
        for (long long j = k; j < k + n; j++) {
            if (i != j)
                Q(i, j) += pen;
            else
                Q(i, j) -= pen * 2;
        }
    }

    for (long long i = 0; i < N; i++) {
        k = i % n;
        for (long long j = k; j < N; j += n) {
            if (i != j) Q(i, j) += pen;
        }
    }

    return pen;
}

double QAP::y(const MatrixXd &Q, const VectorXd &x, double penalty) {
    return penalty * sqrt(Q.innerSize()) * 2 + fQ(Q, x);
}
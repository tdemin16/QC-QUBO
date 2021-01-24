#include "../lib/generators.h"

float NPP::number_partitioning_problem(MatrixXf &Q, vector<int> &nums, int range) {
    if (range < 1) {
        cout << "Range must be >= 1" << endl;
        exit(1);
    }

    int n = nums.size();
    mt19937_64 e_uniform;
    uniform_int_distribution<long> d_int_uniform(1, range);
    random_device rd;
    seed_seq seed{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    float c = 0;
    int prod;

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

float NPP::diff(const MatrixXf &Q, const VectorXf &x, float c) {
    return sqrt(pow(c, 2) + 4 * fQ(Q, x));
}

void NPP::to_file(chrono::duration<double> difference, int it, float f, int n, int range, float diff, const VectorXf &x, const vector<int> nums, const MatrixXf &Q) {
    long unsigned len = 0;
    pid_t pid = time(0);
    ofstream of;
    of.open("../out/" + to_string(pid) + ".txt");

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
        len ++;
        of << to_string(it);
        if(len < nums.size()) of << ", ";
    }
    of << endl << Q;

    of.close();
}

float QAP::quadratic_assignment_problem(MatrixXf &Q, string file) {
    int n;
    long long N;
    float pen;
    float std_dev = 0;
    float mean = 0;
    long long count = 0;
    ifstream in;

    in.open(file);
    if (in.fail()) {
        cout << "File doesn't exist" << endl;
        exit(1);
    }

    in >> n;
    MatrixXf f(n, n);
    MatrixXf d(n, n);
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

    for (long long i = 0; i < N; i++) {
        for (long long j = i + 1; j < N; j++) {
            if (Q(i, j) != 0) {
                mean += Q(i, j);
                count++;
            }
        }
    }

    mean /= count;
    mean *= 2;

    count = 0;
    for (long long i = 0; i < N; i++) {
        for (long long j = i + 1; j < N; j++) {
            if (Q(i, j) != 0) {
                std_dev += pow(Q(i, j) - mean, 2);
                count++;
            }
        }
    }

    std_dev /= count;
    std_dev = sqrt(std_dev);

    pen = (2 * Q.maxCoeff() + std_dev);
    cout << pen << endl;

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

float QAP::y(const MatrixXf &Q, const VectorXf &x, float penalty) {
    return penalty * sqrt(Q.innerSize()) * 2 + fQ(Q, x);
}
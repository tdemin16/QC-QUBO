#include "../lib/generators.h"

float NPP::number_partitioning_problem(MatrixXf &Q, string file) {
    vector<double> nums;
    long long n;
    float c = 0;
    int prod;
    ifstream in;

    in.open(file);
    if (in.fail()) {
        cout << "File doesn't exist" << endl;
        exit(1);
    }

    in >> n;
    nums.reserve(n);
    Q.resize(n, n);

    for (long long i = 0; i < n; i++) {
        in >> nums[i];
        c += nums[i];
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

float NPP::diff(const MatrixXf &Q, const VectorXf &x, float c){
    return sqrt(pow(c, 2) + 4 * fQ(Q, x));
}

float QAP::quadratic_assignment_problem(MatrixXf &Q, string file) {
    int n;
    long long N;
    float pen;
    float std_dev = 0;
    float mean;
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

    mean = Q.mean();
    for (long long i = 0; i < N; i++) {
        for (long long j = i + 1; j < N; j++) {
            if (Q(i, j) != 0) {
                std_dev += pow(Q(i, j) - mean, 2);
                count++;
            }
        }
    }

    std_dev /= count - 1;
    pen = (Q.maxCoeff() + sqrt(std_dev)) * 2;
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
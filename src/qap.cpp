#include "../lib/qap.h"

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

void QAP::to_file(chrono::duration<double> difference, int it, int k, float lambda, double penalty, double f, string file, double y, const VectorXd &x, string filename) {
    ofstream of;
    of.open("../out/" + filename + ".txt");

    of << difference.count() << "s" << endl
       << "imax=" << it << endl
       << "k=" << k << endl
       << "lambda=" << lambda << endl
       << "penalty=" << penalty << endl
       << "f=" << f << endl
       << "file=" + file << endl
       << "y=" << y << endl
       << "SOLUTION_VECTOR:" << endl
       << x.transpose() << endl;

    of.close();
}
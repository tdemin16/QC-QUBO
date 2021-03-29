#include "../lib/npp.h"

long long NPP::number_partitioning_problem(MatrixXd &Q, vector<long long> &nums, int range) {
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

long long NPP::diff(const MatrixXd &Q, const VectorXd &x, long long c) {
    return sqrt(pow(c, 2) + 4 * fQ(Q, x));
}

void NPP::to_file(chrono::duration<double> difference, int it, double f, int n, int range, long long diff, const VectorXd &x, const vector<int> nums, string filename) {
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
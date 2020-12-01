#include "../lib/generators.h"

int number_partitioning_problem(MatrixXf &Q, string file) {
    vector<double> nums;
    int n;
    int c = 0;
    int prod;
    ifstream in;
    in.open(file);

    in >> n;
    nums.reserve(n);
    Q.resize(n, n);

    for (int i = 0; i < n; i++) {
        in >> nums[i];
        c += nums[i];
    }

    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            if(i != j) {
                prod = nums[i] *nums[j];
                Q(i, j) = prod;
                Q(j, i) = prod;
            } else {
                Q(i, i) = nums[i] * (nums[i] - c);
            }
        }
    }

    return c;
}
#include "../lib/lib.h"
#define N 2500

//g++ -std=c++2a -o3 -m64 -fopenmp -Wall -Wextra diff.cpp ../src/lib.cpp -o diff

int main() {
    ifstream in;
    VectorXd z(N);
    long long c = 0;
    vector<long long> num(N);
    double diff2;
    long long prod;
    string s;
    MatrixXd Q(N, N);
    in.open("100000.txt");
    if (!in.fail()) {
        cout << "ex" << endl;

        for (int i = 0; i < N; i++) {
            in >> z(i);
        }

        for (int i = 0; i < N; i++) {
            in >> num[i];
            c += num[i];
        }
        for (long long i = 0; i < N; i++) {
            for (long long j = i; j < N; j++) {
                if (i != j) {
                    prod = num[i] * num[j];
                    Q(i, j) = prod;
                    Q(j, i) = prod;
                    if (prod < 0) cout << "war" << endl;
                } else {
                    Q(i, i) = num[i] * (num[i] - c);
                }
            }
        }

        diff2 = pow(c, 2) + 4 * fQ(Q, z);
        cout.precision(20);
        cout << fQ(Q, z) << endl;
        cout << pow(c, 2) << endl;
        cout << diff2 << endl;
        cout << sqrt(diff2) << endl;

        in.close();
    }

    return 0;
}
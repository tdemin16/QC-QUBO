#include "../lib/lib.h"

//g++ -std=c++2a -o3 -m64 -fopenmp -Wall -Wextra diff.cpp ../src/lib.cpp -o diff

int main() {
    ifstream in;
    VectorXd z(5000);
    long long c = 0;
    vector<int> num(5000);
    double diff2;
    int prod;
    string s;
    MatrixXd Q(5000, 5000);
    in.open("range_10000.txt");
    if (!in.fail()) {
        cout << "ex" << endl;

        for (int i = 0; i < 5000; i++) {
            in >> z(i);
        }

        for (int i = 0; i < 5000; i++) {
            in >> num[i];
            c += num[i];
        }

        for (long long i = 0; i < 5000; i++) {
            for (long long j = i; j < 5000; j++) {
                if (i != j) {
                    prod = num[i] * num[j];
                    Q(i, j) = prod;
                    Q(j, i) = prod;
                } else {
                    Q(i, i) = num[i] * (num[i] - c);
                }
            }
        }

        diff2 = pow(c, 2) + 4 * fQ(Q, z);
        cout.precision(20);
        cout << fQ(Q, z) <<endl;
        cout << pow(c, 2) << endl;
        cout << diff2 << endl;
        cout << sqrt(diff2) << endl;

        in.close();
    }

    return 0;
}
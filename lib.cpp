#include "lib.h"

MatrixXf g(MatrixXf P, int n, float pr) {
    float prob;
    srand(time(0));
    map<int, int> m;
    MatrixXf P_first(n, n);

    // with probability pr inserts i in the map with key i
    for (int i = 0; i < n; i++) {
        prob = (rand() % (10000000)) / 10000000.0f;  // generate a random number from 0 to 1
        if (prob <= pr) {                            // checks if the extracted numeber is equal or less than pr
            m.insert(pair<int, int>(i, i));          //if it is, inserts the number
        }
    }

    // shuffle the keys
    shuffle(m);

    auto end = m.end();
    for (int i = 0; i < n; i++) {
        auto search = m.find(i);
        if (search != end) {
            for (int j = 0; j < n; j++) {
                P_first(i, j) = P(m.at(i), j); //if m contains i, then take the vector Pi where i = m.at(i)
            }
        } else {
            for (int j = 0; j < n; j++) {
                P_first(i, j) = P(i, j); //if m doesn't conatain i, then take the vector Pi
            }
        }
    }

    //P_first is a permutation of P
    return P_first;
}

void shuffle(map<int, int> &m) {
    srand(time(0));
    vector<int> keys; //Vector containing keys of m

    for (auto i : m) {
        keys.push_back(i.first); //add m's keys to keys vector
    }
    random_shuffle(keys.begin(), keys.end()); //to be improved

    vector<int>::iterator it = keys.begin();
    //substitute old keys with new ones (shuffled)
    for (auto &i : m) {
        int ts = i.second;
        i.second = m[*it];
        m[*it] = ts;
        it++;
    }
}

MatrixXf hadamard_product(MatrixXf P, MatrixXf Q, int n) {
    MatrixXf R(n, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            R(i, j) = P(i, j) * Q(i, j);
        }
    }

    return R;
}

float fQ(MatrixXf Q, VectorXf x) {
    return x.transpose() * Q * x;
}

MatrixXf kronecker_product(VectorXf z, int n) {
    MatrixXf M(n, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            M(i, j) = z(i) * z(j);
        }
    }

    return M;
}

MatrixXf diag(VectorXf z, int n) {
    MatrixXf M(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                M(i, j) = z(i);
            } else {
                M(i, j) = 0;
            }
        }
    }

    return M;
}

void h(VectorXf &z, int n, float pr) {
    float random;
    for(int i = 0; i< n; i++) {
        random = (rand() % 10000000) / 10000000.0f;
        if(random <= pr) {
            z(i) = -z(i);
        }
    }
}

float simulated_annealing(float f_first, float f_star, float p) {
    float T = -1 / log(p);
    return exp(-(f_first - f_star) / T);
}
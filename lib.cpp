#include "lib.h"

MatrixXf g(MatrixXf P, float pr) {
    int n = P.rows();
    map<int, int> m;  // Room for improvements
    MatrixXf P_first(n, n);

    random_device rd;
    unsigned int seed = rd() ^ rd();
    mt19937 uniform;
    uniform.seed(seed);
    uniform_real_distribution<float> real_distr(0.0, 1.0);

    // with probability pr inserts i in the map with key i
    for (int i = 0; i < n; i++) {
        if (real_distr(uniform) <= pr) {     // checks if the extracted numeber is equal or less than pr
            m.insert(pair<int, int>(i, i));  //if it is, inserts the number
        }
    }

    // shuffle the keys
    shuffle(m);

    auto end = m.end();
    for (int i = 0; i < n; i++) {
        auto search = m.find(i);
        if (search != end) {
            for (int j = 0; j < n; j++) {
                P_first(i, j) = P(m.at(i), j);  //if m contains i, then take the vector Pi where i = m.at(i)
            }
        } else {
            for (int j = 0; j < n; j++) {
                P_first(i, j) = P(i, j);  //if m doesn't conatain i, then take the vector Pi
            }
        }
    }

    //P_first is a permutation of P
    return P_first;
}

void shuffle(map<int, int> &m) {
    vector<int> keys;  //Vector containing keys of m

    for (auto i : m) {
        keys.push_back(i.first);  //add m's keys to keys vector
    }
    random_shuffle(keys.begin(), keys.end());  //to be improved

    vector<int>::iterator it = keys.begin();
    //substitute old keys with new ones (shuffled)
    for (auto &i : m) {
        int ts = i.second;
        i.second = m[*it];
        m[*it] = ts;
        it++;
    }
}

float fQ(MatrixXf Q, VectorXf x) {
    return x.transpose() * Q * x;
}

void h(VectorXf &z, float pr) {
    int n = z.cols();

    random_device rd;
    unsigned int seed = rd() ^ rd();
    mt19937 uniform;
    uniform.seed(seed);
    uniform_real_distribution<float> real_distr(0.0, 1.0);

    for (int i = 0; i < n; i++) {
        if (real_distr(uniform) <= pr) z(i) = -z(i);
    }
}

float simulated_annealing(float f_first, float f_star, float p) {
    float T = -1 / log(p);
    return exp(-(f_first - f_star) / T);
}

float min(float lambda0, int i, int e) {
    float lambda_first = lambda0 / (2 + i - e);
    
    if(lambda0 < lambda_first) return lambda0;
    return lambda_first;
}
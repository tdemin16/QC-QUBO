#include "lib.h"

MatrixXi g(MatrixXi P, int n, float pr) {
    float prob;
    srand(time(0));
    map<int, int> m;

    for (int i = 0; i < n; i++) {
        prob = (rand() % (10000000)) / 10000000.0f; // generate a random number
        if(prob <= pr) {
            m.insert(i, i);
        }
    }
}
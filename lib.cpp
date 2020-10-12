#include "lib.h"

void make_identity(map<vector<int>, int> P, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j)
                P.insert(pair<vector<int>, int>({i, j}, 1));
            else
                P.insert(pair<vector<int>, int>({i, j}, 0));
        }
    }
    /* Prints the map
    for (auto it : P) {
        cout << "{";
        for (auto vec : it.first) {
            cout << " " << vec;
        }
        cout << " }, ";
        cout << it.second << endl;
    }*/
}

map<vector<int>, int> g(map<vector<int>, int> P, int n, float pr) {
    float prob;
    srand(time(0));
    map<int, int> m;
    map<vector<int>, int> P_first;

    for (auto it : P) {
        cout << "{";
        for (auto vec : it.first) {
            cout << " " << vec;
        }
        cout << " }, ";
        cout << it.second << endl;
    }

    for (int i = 0; i < n; i++) {
        prob = (rand() % (10000000)) / 10000000.0f;  // generate a random number
        if (prob <= pr) {
            m.insert(pair<int, int>(i, i));
        }
    }
    
    auto end = m.end();
    for (int i = 0; i < n; i++) {
        auto search = m.find(i);
        if (search != end) {
            for (int j = 0; j < n; j++) {
                //cout << P.at({i, j}) << endl;
                //P_first.insert(pair<vector<int>, int>({i, j}, P.at(vector<int>(m.at(i), j))));
            }
        } else {
            for (int j = 0; j < n; j++) {
                //P_first.insert(pair<vector<int>, int>({i, j}, P.at(vector<int>(i, j))));
            }
        }
    }

    return P_first;
}
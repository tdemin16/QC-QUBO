#include "lib.h"

void make_identity(map<vector<int>, int> &P, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j)
                P.insert(pair<vector<int>, int>({i, j}, 1));  // Inserts 1 on the diagonal
            else
                P.insert(pair<vector<int>, int>({i, j}, 0));  // Inserts 0 on the other parts
        }
    }
}

map<vector<int>, int> g(map<vector<int>, int> P, int n, float pr) {
    float prob;
    srand(time(0));
    map<int, int> m;
    map<vector<int>, int> P_first;

    // with probability pr inserts i in the map with key i
    for (int i = 0; i < n; i++) {
        prob = (rand() % (10000000)) / 10000000.0f;  // generate a random number from 0 to 1
        if (prob <= pr) {                            // checks if the extracted numeber is equal or less than pr
            m.insert(pair<int, int>(i, i));          //if it is, inserts the number
        }
    }

    // shuffle the matrix P
    shuffle(P, n);

    auto end = m.end();
    for (int i = 0; i < n; i++) {
        auto search = m.find(i);
        if (search != end) {
            for (int j = 0; j < n; j++) {
                P_first.insert(pair<vector<int>, int>({i, j}, P.at({m.at(i), j}))); //if i is a key of m, inserts the line vector P m[i] in Pi'
            }
        } else {
            for (int j = 0; j < n; j++) {
                P_first.insert(pair<vector<int>, int>({i, j}, P.at({i, j}))); //if it isn't,  inserts the line vector Pi in Pi'
            }
        }
    }

    return P_first;
}

void shuffle(map<vector<int>, int> &P, int n) {
    //print_map(P, n);
    srand(time(0));
    vector<vector<int>> keys;

    for (auto i : P) {
        keys.push_back(i.first);
    }
    random_shuffle(keys.begin(), keys.end());

    vector<vector<int>>::iterator it = keys.begin();
    for (auto &i : P) {
        int ts = i.second;
        i.second = P[*it];
        P[*it] = ts;
        it++;
    }
    //cout << endl <<endl;
    //print_map(P, n);
}

void print_map(map<vector<int>, int> P, int n) {
    for (auto it : P) {
        cout << "{";
        for (auto vec : it.first) {
            cout << " " << vec;
        }
        cout << " }, ";
        cout << it.second << endl;
    }
}

map<vector<int>, int> transpose(map<vector<int>, int> &P, int n) {
    //print_map(P, n);
    int temp;

    // swap each element with its transposed 
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            temp = P.at({i, j});
            P[{i, j}] = P.at({j, i});
            P.at({j, i}) = temp;
        }
    }
    //cout << endl;
    //print_map(P, n);
}
#include "bilayer.h"

dMatrix solution(const vector<double> &w){
    dMatrix a(w.size(), vector<double> (w.size(), 0));
    vector<int> index;
    int temp;
    
    for (unsigned i = 0; i < w.size(); ++i) {
        index.push_back(i);
    }
    for (unsigned i = 0; i < w.size()-1; ++i) {
        for (unsigned j = 0; j < w.size()-1-i ; ++j) {
            if (w[index[w.size()-1-j]] > w[index[w.size()-2-j]]) {
                temp = index[w.size()-1-j];
                index[w.size()-1-j] = index[w.size()-2-j];
                index[w.size()-2-j] = temp;
            }
        }
    }
    if (w[index[0]] > (accumulate(w.begin(), w.end(), 0.0) - w[index[0]])) {
        a[index[0]][index[0]] = (w[index[0]] - w[index[1]]-w[index[2]]-w[index[3]])/w[index[0]];
        a[index[0]][index[1]] = w[index[1]]/w[index[0]];
        a[index[0]][index[2]] = w[index[2]]/w[index[0]];
        a[index[0]][index[3]] = w[index[3]]/w[index[0]];
        a[index[1]][index[0]] = 1;
        a[index[2]][index[0]] = 1;
        a[index[3]][index[0]] = 1;
    }
    else{
        a[index[0]][index[1]] = (w[index[0]]+ w[index[1]] - w[index[2]] -w[index[3]])/2/w[index[0]];
        a[index[0]][index[2]] = (w[index[0]]- w[index[1]] + w[index[2]] -w[index[3]])/2/w[index[0]];
        a[index[0]][index[3]] = (w[index[3]])/w[index[0]];
        a[index[1]][index[0]] = (w[index[0]]+ w[index[1]] - w[index[2]] -w[index[3]])/2/w[index[1]];
        a[index[1]][index[2]] = (-w[index[0]]+ w[index[1]] + w[index[2]] +w[index[3]])/2/w[index[1]];
        a[index[2]][index[0]] = (w[index[0]]- w[index[1]] + w[index[2]] -w[index[3]])/2/w[index[2]];
        a[index[2]][index[1]] = (-w[index[0]]+ w[index[1]] + w[index[2]] +w[index[3]])/2/w[index[2]];
        a[index[3]][index[0]] = 1;
        
    }
    return a;
}
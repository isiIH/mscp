#ifndef SCP_H
#define SCP_H

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#include <config.h>
#include <Set.h>

using namespace std;

class SCP {
public:
    Set X;
    vector<vector<int>> F;
    vector<Set> bF;

    ulong n, m, nWX, nWF;

    SCP();
    SCP(const string filename);
    SCP(const SCP& other);
    SCP& operator=(const SCP& other);

    void readFile(const string filename);
    void readFileScp(const string filename);
    void readFilePartition(const string filename);
    void analyzeF();
    void printSubsets(const vector<Set> &C);

};

#endif
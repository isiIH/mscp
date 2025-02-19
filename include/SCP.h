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

    vector<int> uniqueSets;
    Set excludedSets;

    ulong n, m, nWX, nWF;

    SCP();
    SCP(string filename);

    void readFile(string filename);
    void readFileScp(string filename);
    void readFilePartition(string filename);
    void analyzeF();
    void printSubsets(const vector<Set> &C);

};

#endif
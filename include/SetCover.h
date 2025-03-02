#ifndef SET_COVER_H
#define SET_COVER_H

#include <iostream>
#include <algorithm>

#include <config.h>
#include <RowCovering.h>
#include <SCP.h>

using namespace std;

class SetCover {
public:
    vector<int> solution;
    Set U;
    vector<RowCovering> rowMap;
    SCP scp;

    vector<int> uniqueSets;
    Set excludedSets;

    SetCover();
    SetCover(SCP &scp);

    void preprocess();
    void rowReduction();
    void columnDomination();
    void updateRowMap(const int setIndex);
    void push_back(const int s);
    void erase(const int s);
    Set unionSets(const int ignoreSet = -1);
    bool isCovered(const int ignoreSet = -1);
    int size();
    void printRowMap();
};

#endif
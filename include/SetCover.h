#ifndef SET_COVER_H
#define SET_COVER_H

#include <iostream>
#include <algorithm>
#include <execution>

#include <config.h>
#include <RowCovering.h>
#include <SCP.h>
#include <Group.h>
#include <UnionFind.h>

using namespace std;

class SetCover {
public:
    vector<int> solution;
    Set U;
    vector<RowCovering> rowMap;
    SCP scp;
    // Group g;
    UnionFind g;

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
    bool isCovered(const Set& X, const int ignoreSet = -1);
    int size();
    void printRowMap();
};

#endif
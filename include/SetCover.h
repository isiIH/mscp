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
    vector<RowCovering> rowMap;
    SCP scp;
    vector<int> uniqueSets;
    Set excludedSets;

    SetCover();
    SetCover(SCP &scp);

    void preprocess();
    void rowReduction();
    void columnDomination();
    void updateRowMap(int setIndex);
    Set unionSets();
    bool isCovered();
    int size();
    void printRowMap();
};

#endif
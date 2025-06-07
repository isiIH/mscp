#ifndef SET_COVER_H
#define SET_COVER_H

#include <iostream>
#include <algorithm>
#include <execution>
#include <chrono>

#include <config.h>
#include <Edge.h>
#include <RowCovering.h>
#include <SCP.h>

using namespace std;

class SetCover {
public:
    vector<int> solution;
    Set X;
    Set U;
    vector<RowCovering> rowMap;
    SCP* scp;

    vector<int> uniqueSets;
    Set excludedSets;

    SetCover();
    SetCover(SCP &scp);
    SetCover(const SetCover& other);
    SetCover& operator=(const SetCover& other);

    void preprocess();
    void rowReduction();
    void columnDomination();
    void updateRowMap(const int setIndex);
    void push_back(const int s);
    void erase(const int s);
    Set unionSets();
    bool isCovered();
    void redundantSets();
    int size();
    void printRowMap();
};

#endif
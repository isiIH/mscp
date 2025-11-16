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
    vector<int> solution; // Stores the indices of the subsets included in the solution.
    Set X; // Universe
    Set U; // Residual Universe (uncovered elements)
    vector<RowCovering> rowMap; // Maps each element to the subsets that cover it
    SCP* scp; // Pointer to the SCP problem instance handler

    // Preprocess
    vector<int> uniqueSets; // Stores essential sets (those covering grade-1 elements)
    Set excludedSets; // Stores subsets found to be dominated by others

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

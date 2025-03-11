#ifndef GRASP_H
#define GRASP_H

#include <iostream>
#include <vector>
#include <set>
#include <numeric>
#include <cmath>
#include <chrono>
#include <omp.h>
#include <assert.h>

#include <config.h>
#include <Set.h>
#include <SCP.h>
#include <SetCover.h>

using namespace std;

class Grasp {
public:
    SCP scp;
    SetCover bestSol;

    Grasp(SCP &scp);

    SetCover search();
    void searchPerGroup(SetCover& solution);
    void updateSolution(SetCover& solution, const vector<RowCovering>& rowMap, bool &improve);
    void randSuccintSC(SetCover &C, const bool& improve);
};

#endif
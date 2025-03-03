#ifndef GRASP_H
#define GRASP_H

#include <iostream>
#include <vector>
#include <set>
#include <numeric>
#include <cmath>
#include <chrono>

#include <config.h>
#include <Set.h>
#include <SCP.h>
#include <SetCover.h>

using namespace std;

class Grasp {
public:
    SCP scp;

    int function;
    bool improve;
    int last_visited;

    SetCover bestSol;

    Grasp(SCP &scp);

    SetCover search();
    void updateSolution(SetCover& solution, const vector<RowCovering>& rowMap);
    void randSuccintSC(SetCover &C);
};

#endif
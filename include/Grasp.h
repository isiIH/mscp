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
#include <random>

#include <config.h>
#include <Set.h>
#include <SCP.h>
#include <SetCover.h>
#include <Group.h>

using namespace std;

class Grasp {
public:
    SCP scp;
    SetCover bestSol;
    Group g;

    Grasp(SCP &scp);

    SetCover search();
    void searchPerGroup(SetCover& solution);
    void updateSolution(SetCover& solution, const vector<RowCovering>& rowMap, bool &improve, mt19937& gen);
    void randSuccintSC(SetCover &C, const bool& improve, mt19937& gen);
};

#endif
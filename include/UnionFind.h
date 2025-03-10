
#ifndef UNION_FIND_H
#define UNION_FIND_H

#include <vector>
#include <iostream>
#include <omp.h>
#include <unordered_map>

#include <Set.h>

using namespace std;

class UnionFind {
private:
    vector<int> parent, rank;
public:
    vector<vector<int>> groups;

    UnionFind();
    UnionFind(const int m);

    void unite(const int u, const int v);
    int find(const int u);
    void findGroups(const Set& excluded);
    int sizeGroups();
    void printGroups();
};

#endif
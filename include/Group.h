
#ifndef UNION_FIND_H
#define UNION_FIND_H

#include <vector>
#include <iostream>
#include <omp.h>
#include <execution>
#include <algorithm>
#include <limits>
#include <map>
#include <chrono>

#include <Edge.h>
#include <Set.h>
#include <RowCovering.h>

using namespace std;

class UnionFind {
public:
    vector<int> parent, rank;
    UnionFind();
    UnionFind(const int n);
    int find(const int u);
    bool unite(const int u, const int v);
};

class Group {
private:
    int n, nW;
    bool type; // 0:UNION-FIND, 1:MST
    UnionFind uf;
    vector<Edge> mst;
    vector<int> subtreeW;
    vector<int> edgeW;
    map<int, vector<pair<int, int>>> adj;
    int totalWeight = 0;
public:
    vector<Set> U;
    vector<vector<int>> groups;

    Group();
    Group(const bool type, const int n, const int nW);
    
    void buildMST(vector<Edge>& edges);
    void dfs(const int node, const int parent);
    void calcBestCut();
    void collectGroup(const int u, Set& visited, const int groupId, const int bestId);
    void distributeSubsets();

    void findGroups(vector<Edge>& edges, const vector<RowCovering>& rowMap);
    int sizeGroups();
    void printGroups();

};

#endif
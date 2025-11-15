
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
    SCP* scp;
    bool type; // 0:UNION-FIND, 1:MST
    UnionFind uf;
    vector<Edge> edges;
    vector<Edge> mst; 
    vector<int> subtreeW; // Cumulative weight of each edge if cut
    vector<int> edgeW; // Cumulative weight of each subtree including the edge weight
    map<int, vector<pair<int, int>>> adj;
    int totalWeight = 0;
    vector<int> elemToGroup;
public:
    vector<Set> U;
    vector<vector<int>> groups;
    vector<int> subsetToGroup;
    vector<vector<RowCovering>> groupMap;

    Group();
    Group(const int n, SCP& scp);
    
    void buildMST(const int n, vector<Edge>& edges);
    void dfs(const int node, const int parent);
    void calcBestCut();
    void collectGroup(const int u, Set& visited, const int groupId, const int bestId);

    void findGroups(const vector<RowCovering>& rowMap);

    int sizeGroups();
    void printGroups();

    void createGraph(const vector<RowCovering>& rowMap);
    void distributeSubsets(const Set& excludedSets, const vector<RowCovering>& rowMap);

};

#endif
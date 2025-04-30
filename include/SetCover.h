#ifndef SET_COVER_H
#define SET_COVER_H

#include <iostream>
#include <algorithm>
#include <execution>
#include <chrono>
#include <map>

#include <config.h>
#include <RowCovering.h>
#include <SCP.h>
#include <UnionFind.h>

using namespace std;

class Edge {
public:
    int u, v, w; // node u, node v, weight w
    Edge() {}
    Edge(int u, int v, int w) : u(u), v(v), w(w) {}
    bool operator<(const Edge &other) const {
        return w > other.w;
    }
};

class SetCover {
public:
    vector<int> solution;
    Set U;
    vector<RowCovering> rowMap;
    SCP scp;
    UnionFind g;

    vector<Edge> edges;
    vector<Edge> mst;
    vector<int> subtreeW;
    vector<int> childW;
    map<int, vector<pair<int, int>>> adj;
    int totalWeight = 0;

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

    void buildMST();
    void dfs(int node, int parent);
    void calcBestCut();
    void collectGroup(int u, Set& visited, vector<int>& group, const pair<int,int> &cut);
};

#endif
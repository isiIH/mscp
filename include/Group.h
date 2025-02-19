#ifndef GROUP_H
#define GROUP_H

#include <iostream>
#include <vector>
#include <map>

#include <BasicCDS.h>

using namespace std;
using namespace cds;

class Group {
    vector<ulong*> list_groups;
    ulong* visited;
    int m, nWF;

public:
    map<int, vector<int>> graph; //Store adjacent nodes

    Group(int m);
    Group(int m, int node);
    ~Group();

    void add_edge(int u, int v);
    void create_groups();
    void dfs(int node, ulong* group);
    int groups();
    void print();
    void printGraph();
};

#endif
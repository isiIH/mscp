#ifndef GROUP_H
#define GROUP_H

#include <iostream>
#include <vector>
#include <map>

#include <BasicCDS.h>
#include <Set.h>

using namespace std;
using namespace cds;

class Group {
    vector<Set> list_groups;
    Set visited;
    int m, nWF;

public:
    map<int, vector<int>> graph; //Store adjacent nodes

    Group();
    Group(const int m);

    void add_edge(const int u, const int v);
    void create_groups(const Set& excluded);
    void dfs(const int node, Set& group);
    int groups();
    void print();
    void printGraph();
};

#endif
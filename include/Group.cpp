#include <Group.h>

Group::Group(int m) {
    this->m = m;
    nWF = m/(sizeof(ulong)*8); 
    if(m%(sizeof(ulong)*8) > 0) nWF++;

    visited = new ulong[nWF];
    for(int i=0; i<nWF; i++) visited[i] = 0;
}

Group::Group(int m, int node) : Group(m) {
    setBit64(visited, node);
}

Group::~Group() {
    for (auto group : list_groups) {
        delete[] group;
    }
    delete[] visited;
}

void Group::add_edge(int u, int v) {
    // Set both subset edges (undirected graph)
    graph[u].push_back(v);
    graph[v].push_back(u);
}

void Group::printGraph() {
    cout << "----------------" << endl;
    cout << "Adjacency List" << endl;
    cout << "----------------" << endl;

    for(auto u : graph) {
        cout << "S" << u.first << " -> "; 
        for(int v : u.second) {
            cout << v << " ";
        }
        cout << endl;
    }
}

void Group::create_groups() {
    for(int u=0; u<m; u++) {
        if(!checkBit(visited, u)) { // If node has not been visited yet create a new group
            ulong *group = new ulong[nWF];
            for(int i=0; i<nWF; i++) group[i] = 0;
            dfs(u, group); // Deph-First Search
            list_groups.push_back(group);
        }
    }
}

void Group::dfs(int node, ulong* group) {
    // Set node as visited and add to the created group
    setBit64(visited, node); 
    setBit64(group, node);

    // Check the adjacent nodes
    for(int v : graph[node]) {
        if(!checkBit(visited, v))
            dfs(v, group);
    }

    // for(int v=(node+1); v<m; v++) {
    //     if(!checkBit(visited, v) && graph.check_edge(node, v))
    //         dfs(v, group);
    // }
}

void Group::print() {
    cout << "------" << endl;
    cout << "Groups" << endl;
    cout << "------" << endl;
    
    int g_idx = 0;
    for(ulong* group : list_groups) {
        g_idx++;
        cout << "Group " << g_idx << ": ";
        for(int i=0; i<nWF; i++) {
            printBitsUlong(group[i]);
        }
        cout << endl;
    }
}

int Group::groups() {
    return list_groups.size();
}
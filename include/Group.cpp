#include <Group.h>

Group::Group() {}

Group::Group(const int m) {
    this->m = m;
    nWF = m/(sizeof(ulong)*8); 
    if(m%(sizeof(ulong)*8) > 0) nWF++;
}

void Group::add_edge(const int u, const int v) {
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

void Group::create_groups(const Set& excluded) {
    // copy the excluded sets to ignore the edges
    visited = Set(excluded);

    for(int u=0; u<m; u++) {
        if(!visited.check(u)) { // If node has not been visited yet create a new group
            Set group(nWF);
            dfs(u, group); // Deph-First Search
            list_groups.push_back(group);
        }
    }
}

void Group::dfs(const int node, Set& group) {
    // Set node as visited and add to the created group
    visited.push_back(node);
    group.push_back(node);

    // Check the adjacent nodes
    for(int v : graph[node]) {
        if(!visited.check(v))
            dfs(v, group);
    }
}

void Group::print() {
    cout << "------" << endl;
    cout << "Groups" << endl;
    cout << "------" << endl;
    
    int g_idx = 0;
    for(const Set& group : list_groups) {
        g_idx++;
        cout << "Group " << g_idx << ": ";
        group.print();
        cout << endl;
    }
}

int Group::groups() {
    return list_groups.size();
}
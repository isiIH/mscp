#include <iostream>
#include "../include/BasicCDS.h"
#include <vector>

using namespace std;
using namespace cds;

class Graph {
public:
    map<int, vector<int>> adj_list; //Store adjacent nodes

    Graph() {}

    void add_edge(int u, int v) {
        // Set both subset edges (undirected graph)
        adj_list[u].push_back(v);
        adj_list[v].push_back(u);
    }

    void print() {
        cout << "----------------" << endl;
        cout << "Adjacency List" << endl;
        cout << "----------------" << endl;

        for(auto u : adj_list) {
            cout << "S" << u.first << " -> "; 
            for(int v : u.second) {
                cout << v << " ";
            }
            cout << endl;
        }
    }
};


// class Graph {
//     vector<ulong*> adj_matrix; //Store edges between subsets
//     int row, col;

// public:
//     Graph() {}

//     Graph(int m) {
//         adj_matrix = vector<ulong*>(m);

//         // Bitwise representation of edges
//         int nWF = m/(sizeof(ulong)*8); 
//         if(m%(sizeof(ulong)*8) > 0) nWF++;

//         row = m;
//         col = nWF;

//         // Initialize graph m x nWF
//         for(int i=0; i<m; i++) {
//             adj_matrix[i] = new ulong[nWF];
//             for(int j=0; j<nWF; j++) adj_matrix[i][j] = 0;
//         }
//     }

//     void add_edge(int u, int v) {
//         // Set both subset edges to 1 (undirected graph)
//         setBit64(adj_matrix[u], v);
//         setBit64(adj_matrix[v], u);
//     }

//     bool check_edge(int u, int v) {
//         return checkBit(adj_matrix[u], v);
//     }

//     void print() {
//         cout << "----------------" << endl;
//         cout << "Adjacency Matrix" << endl;
//         cout << "----------------" << endl;

//         for(int u=0; u<row; u++) {
//             cout << "S" << u << " -> ";
//             for(int v=0; v<row; v++) {
//                 if(checkBit(adj_matrix[u], v)) {
//                     cout << v << " ";
//                 }
//             }
//             cout << endl;
//         }
//     }

//     int rows() {
//         return row;
//     }

//     int cols() {
//         return col;
//     }
// };

class Group {
    vector<ulong*> list_groups;
    Graph graph;
    ulong* visited;
    int m, nWF;

public:
    Group(Graph g, int m) {

        graph = g;
        
        // m = graph.rows();
        // nWF = graph.cols();

        this->m = m;
        nWF = m/(sizeof(ulong)*8); 
        if(m%(sizeof(ulong)*8) > 0) nWF++;

        visited = new ulong[nWF];
        for(int i=0; i<nWF; i++) visited[i] = 0;
    }

    Group(Graph g, int m, int node) : Group(g, m) {
        setBit64(visited, node);
    }

    void create_groups() {
        for(int u=0; u<m; u++) {
            if(!checkBit(visited, u)) { // If node has not been visited yet create a new group
                ulong *group = new ulong[nWF];
                for(int i=0; i<nWF; i++) group[i] = 0;
                dfs(u, group); // Deph-First Search
                list_groups.push_back(group);
            }
        }
    }

    void dfs(int node, ulong* group) {
        // Set node as visited and add to the created group
        setBit64(visited, node); 
        setBit64(group, node);

        // Check the adjacent nodes
        for(int v : graph.adj_list[node]) {
            if(!checkBit(visited, v))
                dfs(v, group);
        }

        // for(int v=(node+1); v<m; v++) {
        //     if(!checkBit(visited, v) && graph.check_edge(node, v))
        //         dfs(v, group);
        // }
    }

    void print() {
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

    int groups() {
        return list_groups.size();
    }
};
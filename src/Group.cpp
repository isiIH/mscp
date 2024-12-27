#include <iostream>
#include "../include/BasicCDS.h"
#include <vector>

using namespace std;
using namespace cds;

class Graph {
    vector<ulong*> adj_matrix; //Store edges between subsets
    int row, col;

public:
    Graph() {}

    Graph(int m) {
        adj_matrix = vector<ulong*>(m);

        // Bitwise representation of edges
        int nWF = m/(sizeof(ulong)*8); 
        if(m%(sizeof(ulong)*8) > 0) nWF++;

        row = m;
        col = nWF;

        // Initialize graph m x nWF
        for(int i=0; i<m; i++) {
            adj_matrix[i] = new ulong[nWF];
            for(int j=0; j<nWF; j++) {
                adj_matrix[i][j] = 0;
            }
        }
    }

    void add_edge(int u, int v) {
        // Set both subset edges to 1 (undirected graph)
        setBit64(adj_matrix[u], v);
        setBit64(adj_matrix[v], u);
    }

    bool check_edge(int u, int v) {
        return checkBit(adj_matrix[u], v);
    }

    void print() {
        cout << "----------------" << endl;
        cout << "Adjacency Matrix" << endl;
        cout << "----------------" << endl;

        for(int u=0; u<row; u++) {
            cout << "S" << u << " -> ";
            for(int v=0; v<row; v++) {
                if(checkBit(adj_matrix[u], v)) {
                    cout << v << " ";
                }
            }
            cout << endl;
        }
    }

    int rows() {
        return row;
    }

    int cols() {
        return col;
    }
};

class Group {
    vector<ulong*> list_groups;
    Graph g;
    ulong* visited;
    int m, nWF;

public:
    Group(Graph g) {
        this->g = g;
        m = g.rows();
        nWF = g.cols();

        visited = new ulong[nWF];
        for(int i=0; i<nWF; i++) visited[i] = 0;

        for(int u=0; u<m; u++) {
            if(!checkBit(visited, u)) {
                ulong *group = new ulong[nWF];
                for(int i=0; i<nWF; i++) group[i] = 0;
                dfs(u, group);
                list_groups.push_back(group);
            }
        }
    }

    void dfs(int node, ulong* group) {
        setBit64(visited, node); 
        setBit64(group, node);

        for(int v=0; v<m; v++) {
            if(g.check_edge(node, v) && !checkBit(visited, v)) {
                dfs(v, group);
            }
        }
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
};

int main()
{
    int m = 8;
    Graph g(m);

    g.add_edge(0, 3);
    g.add_edge(0, 5);
    g.add_edge(1, 2);
    g.add_edge(1, 4);
    g.add_edge(4, 2);
    g.add_edge(6, 7);

    g.print();

    Group c(g);
    c.print();

    return 0;
}
#ifndef EDGE_H
#define EDGE_H

class Edge {
public:
    int u, v, w;
    Edge() {}
    Edge(int u, int v, int w) : u(u), v(v), w(w) {}

    bool operator<(const Edge &other) const {
        return w > other.w;
    }
};

#endif
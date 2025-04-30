#include <UnionFind.h>

UnionFind::UnionFind() {}

UnionFind::UnionFind(const int m) {
    parent.resize(m);
    rank.resize(m, 0);

    #pragma omp parallel for
    for(int i=0; i<m; i++)
        parent[i] = i;
}

bool UnionFind::unite(const int u, const int v) {
    int uRoot = find(u);
    int vRoot = find(v);

    // if they are in the same group already
    if(uRoot == vRoot) return false;

    // unite them by rank
    if(rank[uRoot] < rank[vRoot]) {
        parent[uRoot] = vRoot;
    } else if(rank[vRoot] < rank[uRoot]) {
        parent[vRoot] = uRoot;
    } else {
        parent[vRoot] = uRoot;
        rank[uRoot]++;
    }
    return true;
}

int UnionFind::find(const int u) {
    int root = parent[u];

    // if it is not the representative of u
    if(parent[root] != root) {
        return parent[u] = find(parent[root]);
    }

    return root;
}

void UnionFind::findGroups(const Set& excluded) {
    groups.clear();
    unordered_map<int, vector<int>> group_map;

    for (int i = 0; i < parent.size(); i++) {
        if(!excluded.check(i)) {
            int root = find(i);
            // printf("root: %d, i: %d\n", root, i);
            group_map[root].push_back(i);
        }
    }

    for (auto &entry : group_map) {
        groups.push_back(entry.second);
    }
}

int UnionFind::sizeGroups() {
    return groups.size();
}

void UnionFind::printGroups() {
    for(int i=0; i<groups.size(); i++) {
        printf("Group %d: ", (i+1));
        for (const int e : groups[i])
            printf("%d ", e);
        printf("\n");
    }
}
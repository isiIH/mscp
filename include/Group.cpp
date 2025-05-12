#include <Group.h>

/*----------UnionFind----------*/
UnionFind::UnionFind() {}

UnionFind::UnionFind(const int n) {
    parent.resize(n);
    rank.resize(n, 0);

    #pragma omp parallel for
    for(int i=0; i<n; i++)
        parent[i] = i;
}

int UnionFind::find(const int u) {
    // if it is not the representative of u
    if (parent[u] != u)
        parent[u] = find(parent[u]); // Path compression
    return parent[u];
}

bool UnionFind::unite(const int u, const int v) {
    int uRoot = find(u), vRoot = find(v);

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

/*----------Group----------*/
Group::Group() {}

Group::Group(const int n, const int nW) : n(n), nW(nW), type(SEG_TYPE) /*0:UNION-FIND, 1:MST*/ {
    uf = UnionFind(n);
    if(type) { // MST
        subtreeW.resize(n);
        elemToGroup.resize(n);
    }
}

void Group::createGraph(const vector<RowCovering>& rowMap) {
    #pragma omp parallel for
    for(int i=0; i<rowMap.size() - 1; i++) {
        vector<Edge> local_edges;
        for(int j=i+1; j<rowMap.size(); j++) {
            int sharedSubsets = rowMap[i].countIntersection(rowMap[j].col_covering);
            if(sharedSubsets > 0) {
                // Add the edges to the group
                local_edges.push_back(Edge(rowMap[i].row, rowMap[j].row, sharedSubsets));
            }
        }  
        #pragma omp critical
        {
            edges.insert(edges.end(), local_edges.begin(), local_edges.end());
        }
    }
}

void Group::findGroups(const vector<RowCovering> rowMap) {
    createGraph(rowMap);
    if(type) { // MST
        buildMST(edges);
        dfs(mst[0].u, -1);  
        printf("Total edges: %ld\n", edges.size());
        printf("Total Weight: %d\n", totalWeight);
        calcBestCut();
        
    } else { // UNION-FIND
        #pragma omp parallel for
        for(const Edge& e : edges) {
            // printf("Edge: %d %d %d\n", e.u, e.v, e.w);
            uf.unite(e.u, e.v);
        }

        unordered_map<int, vector<int>> group_map;

        for (const RowCovering& e : rowMap) {
            int root = uf.find(e.row);
            // printf("root: %d, i: %d\n", root, i);
            group_map[root].push_back(e.row);
        }

        for (auto &entry : group_map) {
            Set u(nW);
            for(int e : entry.second) {
                u.push_back(e);
                
            }
            U.push_back(u);
        }
    }

    groups.resize(U.size());
}

void Group::distributeSubsets(const vector<Set>& bF, const Set& excludedSets, const vector<RowCovering>& rowMap) {
    vector<Set> covered = U;
    int m = bF.size();

    subsetToGroup.resize(m, -1);
    vector<bool> visited(m, false);
    for(const RowCovering& row : rowMap) {
        int bestSet = -1;
        int bestCover = 0;
        for(int ss : row.col_covering) {
            if(visited[ss] && subsetToGroup[ss] == elemToGroup[row.row]) {
                bestSet = -1;
                break;
            }
            if(visited[ss]) continue;

            int elemsCover = U[elemToGroup[row.row]].intersectionLength(bF[ss]);
            if(elemsCover > bestCover) {
                bestCover = elemsCover;
                bestSet = ss;
            }
        }

        if(bestSet != -1) {
            groups[elemToGroup[row.row]].push_back(bestSet);
            subsetToGroup[bestSet] = elemToGroup[row.row];
            covered[elemToGroup[row.row]].substract(bF[bestSet]);
            visited[bestSet] = true;
        }
    }
    
    for(RowCovering row : rowMap) {
        if(covered[elemToGroup[row.row]].check(row.row)) {
           
            printf("Row %d: Group: %d\n", row.row, elemToGroup[row.row]);
            for(int ss : row.col_covering) {
                printf("Subset %d Group %d\n", ss, subsetToGroup[ss]);
                int elemsCover = U[subsetToGroup[ss]].intersectionLength(bF[ss]);
                printf("Elems Cover: %d\n", elemsCover);
            }
        }
    }
}

int Group::sizeGroups() {
    return U.size();
}

void Group::printGroups() {
    for(Set& u : U) {
        u.print();
    }
    for(int i=0; i<groups.size(); i++) {
        // printf("Group %d: (%ld)", (i+1), groups[i].size());
        for (const int e : groups[i])
            printf("%d ", e);
        printf("\n");
    }
}

void Group::buildMST(vector<Edge>& edges) {
    // Sort edges in descending order
    sort(execution::par, edges.begin(), edges.end());

    int edgeId = 0;
    for(const Edge& e : edges) {
        if(uf.unite(e.u, e.v)) {
            mst.push_back(e);
            totalWeight += e.w;
            adj[e.u].emplace_back(e.v, edgeId);
            adj[e.v].emplace_back(e.u, edgeId);
            if (mst.size() == (n - 1)) break;
            edgeId++;
        }
    }

    // for(const Edge& e : mst) {
    //     printf("(%d <- (%d) -> %d)\n", e.u+1, e.w, e.v+1);
    // }

    printf("MST size: %ld\n", mst.size());

    edgeW.resize(mst.size());
}

void Group::dfs(int node, int parent) {
    subtreeW[node] = 0;
    for(auto& [neighbor, uid] : adj[node]) {
        if (neighbor == parent) continue;
        dfs(neighbor, node);
        edgeW[uid] = subtreeW[neighbor]; // Cumulative weight if the edge is cut
        subtreeW[node] += edgeW[uid] + mst[uid].w; // Cumulative weight of the node
    }
}

void Group::calcBestCut() {
    int bestId, diff, bestDiff = numeric_limits<int>::max();

    for(int i=0; i<mst.size(); i++) {
        diff = abs(2 * edgeW[i] - totalWeight + mst[i].w); // difference in weight of each subtree if the edge is cut
        if(diff < bestDiff) {
            bestDiff = diff;
            bestId = i;
        }
    }
    printf("Best Diff: %d\n", bestDiff);
    printf("Best Id: %d\n", bestId);

    Set visited(nW);
    U.resize(2, Set(nW));
    #pragma omp parallel sections
    {
        #pragma omp section
        collectGroup(mst[bestId].u, visited, 0, bestId);

        #pragma omp section
        collectGroup(mst[bestId].v, visited, 1, bestId);
    }
}

void Group::collectGroup(const int u, Set& visited, const int groupId, const int bestId) {
    visited.push_back(u);
    U[groupId].push_back(u);
    elemToGroup[u] = groupId;
    for (auto& [v, eid] : adj[u]) {
        if (visited.check(v)) continue;
        Edge& e = mst[eid];
        Edge& cutEdge = mst[bestId];
        if ((e.u == cutEdge.u && e.v == cutEdge.v) || (e.u == cutEdge.v && e.v == cutEdge.u))
            continue;
        collectGroup(v, visited, groupId, bestId);
    }
}
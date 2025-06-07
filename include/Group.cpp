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

Group::Group(const int n, const int nW) : nW(nW), type(SEG_TYPE) /*0:UNION-FIND, 1:MST*/ {
    uf = UnionFind(n);
    elemToGroup.resize(n);
    if(type) { // MST
        subtreeW.resize(n);
    }
}

void Group::createGraph(const vector<RowCovering>& rowMap) {
    int rowSize = rowMap.size();
    #pragma omp parallel for
    for(int i=0; i<rowSize - 1; i++) {
        vector<Edge> local_edges;
        int sharedSubsets;
        for(int j=i+1; j<rowSize; j++) {
            sharedSubsets = rowMap[i].countIntersection(rowMap[j].col_covering);
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

void Group::findGroups(const vector<RowCovering>& rowMap) {
    if(PRINT) {
        cout << "------------------------" << endl;
        cout << "Executing FindGroups..." << endl;
        cout << "------------------------" << endl;
    }
    createGraph(rowMap);
    if(type) { // MST
        buildMST(rowMap.size(), edges);
        dfs(mst[0].u, -1);  
        if(CHECK) {
            printf("Total edges: %ld\n", edges.size());
            printf("Total Weight: %d\n", totalWeight);
        }
        calcBestCut();
        
    } else { // UNION-FIND
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

        int groupId = 0;
        for (auto &entry : group_map) {
            Set u(nW);
            for(int e : entry.second) {
                u.push_back(e);
                elemToGroup[e] = groupId;
            }
            U.push_back(u);
            groupId++;
        }
    }

    groups.resize(U.size());
    groupMap.resize(U.size());
}

void Group::distributeSubsets(const vector<Set>& bF, const Set& excludedSets, const vector<RowCovering>& rowMap) {
    int m = bF.size();
    subsetToGroup.resize(m, -1);
    vector<bool> visited(m, false);

    int rowGroup, bestSet, bestGroup, bestCover, elemsCover;
    // Cover each row with the best subset
    for(const RowCovering& row : rowMap) {
        rowGroup = elemToGroup[row.row];
        groupMap[rowGroup].push_back(row);
        bestSet = -1;
        bestCover = 0;
        for(const int ss : row.col_covering) {
            if(subsetToGroup[ss] == rowGroup) {
                bestSet = -1;
                break;
            }
            if(visited[ss]) continue;

            elemsCover = U[rowGroup].intersectionLength(bF[ss]);
            if(elemsCover > bestCover) {
                bestCover = elemsCover;
                bestSet = ss;
            }
        }

        if(bestSet != -1) {
            groups[rowGroup].push_back(bestSet);
            subsetToGroup[bestSet] = rowGroup;
            visited[bestSet] = true;
        }
    }

    // Verify if any subset is not assigned to a group
    for(int i=0; i<m; i++) {
        if(visited[i] || excludedSets.check(i)) continue;

        bestGroup = -1;
        bestCover = 0;
        for(int j=0; j<groups.size(); j++) {
            elemsCover = U[j].intersectionLength(bF[i]);
            if(elemsCover > bestCover) {
                bestCover = elemsCover;
                bestGroup = j;
            }
        }

        if(bestGroup != -1) {
            groups[bestGroup].push_back(i);
            subsetToGroup[i] = bestGroup;
        }
    }

    // Erase subsets that are not assigned to any group
    #pragma omp parallel for
    for(vector<RowCovering>& rowMap : groupMap) {
        for(RowCovering& row : rowMap) {
            rowGroup = elemToGroup[row.row];
            vector<int>& cols = row.col_covering;
            cols.erase(
                remove_if(cols.begin(), cols.end(), [&](int subsetIdx) {
                    return subsetToGroup[subsetIdx] != rowGroup;
                }),
                cols.end()
            );
        }
    }
}

int Group::sizeGroups() {
    return U.size();
}

void Group::printGroups() {
    for(int i=0; i<groups.size(); i++) {
        printf("Group %d: \n", i);
        U[i].print();

        printf("Subsets: ");
        for (const int e : groups[i])
            printf("%d ", e);
        printf("\n");
    }
}

void Group::buildMST(const int n, vector<Edge>& edges) {
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

    if(CHECK) printf("MST size: %ld\n", mst.size());

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
    if(CHECK) {
        printf("Best Diff: %d\n", bestDiff);
        printf("Best Id: %d\n", bestId);
    }

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
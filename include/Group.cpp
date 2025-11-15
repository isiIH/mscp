#include <Group.h>

/*----------UnionFind----------*/
UnionFind::UnionFind() {}

UnionFind::UnionFind(const int n) {
    // Initialize parent and rank arrays
    parent.resize(n);
    rank.resize(n, 0);

    // Each parent points to itself at first
    #pragma omp parallel for
    for(int i=0; i<n; i++)
        parent[i] = i;
}

// Find the representative of the set that contains u
int UnionFind::find(const int u) {
    // if it is not the representative of u
    if (parent[u] != u)
        parent[u] = find(parent[u]); // Path compression
    return parent[u];
}

// Unite the sets that contain u and v by rank
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

Group::Group(const int n, SCP &scp) : scp(&scp), type(SEG_TYPE) /*0:UNION-FIND, 1:MST*/ {
    uf = UnionFind(n); // Initialize Union-Find
    elemToGroup.resize(n, -1); // Map elements to the group they belong
    if(type) { // MST
        subtreeW.resize(n);
    }
}

// Build the Graph that maps the relationships between rows using countIntersection
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

// Find groups using UNION-FIND or MST
void Group::findGroups(const vector<RowCovering>& rowMap) {
    if(PRINT) {
        cout << "------------------------" << endl;
        cout << "Executing FindGroups..." << endl;
        cout << "------------------------" << endl;
    }

    auto start_time = chrono::high_resolution_clock::now();
    createGraph(rowMap);
    auto end_time = chrono::high_resolution_clock::now();
    if(PRINT) printf("Time createGraph: %f\n", chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0);

    if(type) { // MST
        buildMST(rowMap.size(), edges); // Build Maximum Spanning Tree
        dfs(mst[0].u, -1); // Calculate subtree weights
        if(CHECK) {
            printf("Total edges: %ld\n", edges.size());
            printf("Total Weight: %d\n", totalWeight);
        }
        calcBestCut(); // Calculate the best cut to form groups
        
    } else { // UNION-FIND
        auto start_time = chrono::high_resolution_clock::now();
        // Unite nodes based on the edges
        for(const Edge& e : edges) {
            // printf("Edge: %d %d %d\n", e.u, e.v, e.w);
            uf.unite(e.u, e.v);
        }

        unordered_map<int, vector<int>> group_map;

        // Map each row to its group representative
        for (const RowCovering& e : rowMap) {
            int root = uf.find(e.row); 
            // printf("root: %d, i: %d\n", root, i);
            group_map[root].push_back(e.row);
        }

        // Represent the universe of each group and map elements to their group
        int groupId = 0;
        for (auto &entry : group_map) {
            Set u(scp->nWX);
            for(int e : entry.second) {
                u.push_back(e);
                elemToGroup[e] = groupId;
            }
            U.push_back(u);
            groupId++;
        }
        auto end_time = chrono::high_resolution_clock::now();
        if(PRINT) printf("Time union-find: %f\n", chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0);
    }

    groups.resize(U.size());
    groupMap.resize(U.size());
}

// Distribute subsets to the groups already formed
void Group::distributeSubsets(const Set& excludedSets, const vector<RowCovering>& rowMap) {
    if(type) { // MST
        subsetToGroup.resize(scp->m, -1);
        vector<bool> visited(scp->m, false);

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

                elemsCover = U[rowGroup].intersectionLength(scp->bF[ss]);
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
        for(int i=0; i<scp->m; i++) {
            if(visited[i] || excludedSets.check(i)) continue;

            bestGroup = -1;
            bestCover = 0;
            for(int j=0; j<groups.size(); j++) {
                elemsCover = U[j].intersectionLength(scp->bF[i]);
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
                int rowGroup = elemToGroup[row.row];
                vector<int>& cols = row.col_covering;
                cols.erase(
                    remove_if(cols.begin(), cols.end(), [&](int subsetIdx) {
                        return subsetToGroup[subsetIdx] != rowGroup;
                    }),
                    cols.end()
                );
            }
        }
    } else { // UNION-FIND
        // Distribute subsets based on the elemToGroup mapping
        for(int i=0; i<scp->m; i++) {
            if(excludedSets.check(i)) continue;
            int rowGroup = -1;
            for(int j=0; j<scp->F[i].size(); j++) {
                if(elemToGroup[scp->F[i][j] - 1] != -1) {
                    rowGroup = elemToGroup[scp->F[i][j] - 1];
                    break;
                }
            }
            groups[rowGroup].push_back(i);
        }
        // Map each row to its group representative
        for(const RowCovering& row : rowMap) {
            groupMap[elemToGroup[row.row]].push_back(row);
        }
    }
}

// Return the number of groups found
int Group::sizeGroups() {
    return U.size();
}

// Print groups found
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

// Build Maximum Spanning Tree using Union-Find
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

// Depth-First Search to calculate subtree weights
void Group::dfs(int node, int parent) {
    subtreeW[node] = 0;
    for(auto& [neighbor, uid] : adj[node]) {
        if (neighbor == parent) continue;
        dfs(neighbor, node);
        edgeW[uid] = subtreeW[neighbor]; // Cumulative weight if the edge is cut
        subtreeW[node] += edgeW[uid] + mst[uid].w; // Cumulative weight of the node
    }
}

// Calculate the best cut in the MST to form two groups
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

    // Identify the universe of each group after the best cut
    Set visited(scp->nWX);
    U.resize(2, Set(scp->nWX));
    #pragma omp parallel sections
    {
        #pragma omp section
        collectGroup(mst[bestId].u, visited, 0, bestId);

        #pragma omp section
        collectGroup(mst[bestId].v, visited, 1, bestId);
    }
}

// Collect elements of a group using DFS
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
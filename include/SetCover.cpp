#include <SetCover.h>

SetCover::SetCover() {};

SetCover::SetCover(SCP &scp) : scp(scp) {
    if(GROUP_SEG) {
        g = UnionFind(scp.m);
    }
    excludedSets = Set(scp.nWF);
    U = Set(scp.X);
    rowMap.resize(scp.n);
    preprocess();
}

void SetCover::preprocess() {
    if(PRINT) {
        cout << "------------------------" << endl;
        cout << "Executing PreSetCover..." << endl;
        cout << "------------------------" << endl;
    }

    // Column Domination
    auto start_time = chrono::high_resolution_clock::now();
    columnDomination();
    auto end_time = chrono::high_resolution_clock::now();
    printf("Time Column Domination: %f\n", chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0);

    // Create map structure
    start_time = chrono::high_resolution_clock::now();
    if(PRINT) printf("Creating Row Map...\n");
    #pragma omp parallel for schedule(dynamic, 1)
    for(int i=0; i<scp.n; i++)
        rowMap[i] = RowCovering(excludedSets, scp.bF, i);

    sort(execution::par, rowMap.begin(), rowMap.end(), [&](RowCovering a, RowCovering b){return a.n_columns < b.n_columns;});
    end_time = chrono::high_resolution_clock::now();
    printf("Time Create Map: %f\n", chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0);

    printRowMap();
    // Add uniques elements
    start_time = chrono::high_resolution_clock::now();
    rowReduction();
    end_time = chrono::high_resolution_clock::now();
    printf("Time RowReduction: %f\n", chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0);

    printRowMap();

    // Universe Segmentation
    start_time = chrono::high_resolution_clock::now();
    // if(GROUP_SEG) {
    //     #pragma omp parallel for
    //     for(const Edge& e : edges) {
    //         // printf("Edge: %d %d %d\n", e.u, e.v, e.w);
    //         if(!excludedSets.check(e.u) && !excludedSets.check(e.v)) g.unite(e.u, e.v);
    //     }
    //     g.findGroups(excludedSets);
    //     printf("Groups: %d\n", g.sizeGroups());
    // }
    // g.printGroups();
    if(GROUP_SEG) {
        buildMST();
        subtreeW.resize(scp.m);
        childW.resize(mst.size());
        dfs(mst[0].u, -1);  
        // print adjacency list
        // printf("Adjacency List:\n");
        // for(int i=0; i<adj.size(); i++) {
        //     printf("%d: ", i);
        //     for(pair<int, int> &e : adj[i]) {
        //         printf("(%d, %d) ", e.first, e.second);
        //     }
        //     printf("\n");
        // }
        // print subtreeW();
        // printf("SubtreeW: ");  
        // for(int i=0; i<subtreeW.size(); i++) {
        //     printf("%d ", subtreeW[i]);
        // }
        // printf("\n");
        // printf("ChildW: ");
        // for(int i=0; i<childW.size(); i++) {
        //     printf("%d ", childW[i]);
        // }
        // printf("\n");
        printf("Total edges: %ld\n", edges.size());
        printf("Total Weight: %d\n", totalWeight);
        calcBestCut();
    }
    end_time = chrono::high_resolution_clock::now();
    printf("Time Segmentation: %f\n", chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0);

    if(1) {
        cout << "Added " << uniqueSets.size() << " subsets" << endl; 
        cout << "Excluded " << excludedSets.size() - uniqueSets.size() << " subsets" << endl;
        cout << "|X| = " << rowMap.size() << endl;
    }
}

void SetCover::rowReduction() {
    if(PRINT) printf("Executing Row Reduction...\n");
    int setIndex;

    // Check if the sorted rowMap has unique elements
    while(!rowMap.empty() && rowMap[0].n_columns == 1) {
        setIndex = rowMap[0].col_covering[0];

        // Add subset of grade 1
        uniqueSets.push_back(setIndex);
        // push_back(setIndex);
        updateRowMap(setIndex);

        excludedSets.push_back(setIndex);
    }
}

void SetCover::columnDomination() {
    if(PRINT) printf("Executing Column Domination...\n");

    // sort the subsets in ascending order
    vector<pair<int, int>> indexedSubsets(scp.m);

    #pragma omp parallel for
    for(int i=0; i<scp.m; i++)
        indexedSubsets[i] = {i, scp.bF[i].size()};

    sort(execution::par, indexedSubsets.begin(), indexedSubsets.end(), [](pair<int, int> a, pair<int, int> b) {
        return a.second < b.second;
    });

    #pragma omp parallel shared(indexedSubsets) 
    {

        vector<Edge> local_edges;
        #pragma omp for schedule(dynamic, 1) nowait
        for(int i=0; i < scp.m-1; i++) {
            bool ignore = false;
            int nIntersect;
            pair<int, int> setA = indexedSubsets[i];
            int setB;

            for(int j=i+1; j < scp.m; j++) {
                setB = indexedSubsets[j].first;

                nIntersect = scp.bF[setA.first].intersectionLength(scp.bF[setB]);

                // If the intersection is the same size of the smallest subset 
                if(nIntersect == setA.second) {
                    #pragma omp critical
                    {
                        excludedSets.push_back(setA.first);
                    }
                    ignore = true;
                    break;
                }

                // Save the neightbors if the subset A is not excluded
                if(GROUP_SEG && nIntersect) local_edges.push_back(Edge(setA.first, setB, nIntersect));
            }
            
            // if(ignore) {
            //     local_edges.clear();
            // }
        }

        #pragma omp critical
        {
            edges.insert(edges.end(), local_edges.begin(), local_edges.end());
        }

    }
}

void SetCover::updateRowMap(const int setIndex) {
    int nRows = rowMap.size();
    int aux = 0, aux2;
    // Remove all the subset's elements from the rowMap
    for(int i=0; i<nRows; i++) {
        aux2 = i - aux;
        if(scp.bF[setIndex].check(rowMap[aux2].row)) {
            rowMap.erase(rowMap.begin() + aux2);
            aux++;
        }
    }

    U.substract(scp.bF[setIndex]);
}

void SetCover::push_back(const int s) {
    solution.push_back(s);
    updateRowMap(s);
}

void SetCover::erase(const int s) {
    solution.erase(solution.begin() + s);
}

Set SetCover::unionSets(const int ignoreSet) {
    Set C(scp.nWX);
    for(const int idS : solution) {
        if(idS != ignoreSet) {
            for(int i=0; i<scp.nWX; i++)
                C.S[i] |= scp.bF[idS].S[i];
        }
    }
    return C;
}

bool SetCover::isCovered(const Set& X, const int ignoreSet) {
    Set coveredElements = unionSets(ignoreSet);
    
    for (int i = 0; i < scp.nWX; i++) if ((coveredElements.S[i] & X.S[i]) != X.S[i]) {
        return false;
    }
    return true;
}

int SetCover::size() {
    return solution.size();
}

void SetCover::printRowMap() {
    if(CHECK) {
        for(RowCovering row : rowMap) {
            printf("(%d) |%d| => ", row.row, row.n_columns);
            for (int index : row.col_covering)
                printf("%d ", index);
            printf("\n");
        }
    }
}

void SetCover::buildMST() {
    UnionFind uf = UnionFind(scp.m);

    // Sort edges in descending order
    sort(edges.begin(), edges.end());

    int edgeId = 0;
    for(const Edge& e : edges) {
        if(!excludedSets.check(e.u) && !excludedSets.check(e.v) && uf.unite(e.u, e.v)) {
            mst.push_back(e);
            totalWeight += e.w;
            adj[e.u].emplace_back(e.v, edgeId);
            adj[e.v].emplace_back(e.u, edgeId);
            edgeId++;
            if (mst.size() == scp.m - 1) break;
        }
    }

    // for(const Edge& e : mst) {
    //     printf("(%d <- (%d) -> %d)\n", e.u+1, e.w, e.v+1);
    // }

    printf("MST size: %ld\n", mst.size());
}

void SetCover::dfs(int node, int parent) {
    subtreeW[node] = 0;
    for(pair<int, int> &e : adj[node]) {
        if (e.first == parent) continue;
        dfs(e.first, node);
        // peso de la rama que "cuelga" por e.eid:
        childW[e.second] = subtreeW[e.first];
        // acumula ese peso en el subárbol de v
        subtreeW[node] += childW[e.second] + mst[e.second].w;
    }
}

void SetCover::calcBestCut() {
    int bestDiff = numeric_limits<int>::max();;
    int bestId = -1;
    int diff;
    for(int i=0; i<mst.size(); i++) {
        diff = llabs(2 * childW[i] - totalWeight + mst[i].w);
        if(diff < bestDiff) {
            bestDiff = diff;
            bestId = i;
        }
    }
    printf("Best Diff: %d\n", bestDiff);
    printf("Best Id: %d\n", bestId);

    Set visited(scp.m);
    vector<vector<int>> groups(2);
    collectGroup(mst[bestId].u, visited, groups[0], {mst[bestId].u, mst[bestId].v});
    collectGroup(mst[bestId].v, visited, groups[1], {mst[bestId].u, mst[bestId].v});

    // print groups
    // printf("Group 1: ");
    // for(int i=0; i<groups[0].size(); i++) {
    //     printf("%d ", groups[0][i]);
    // }
    // printf("\n");
    // printf("Group 2: ");
    // for(int i=0; i<groups[1].size(); i++) {
    //     printf("%d ", groups[1][i]);
    // }
    // printf("\n");
}

void SetCover::collectGroup(int u, Set& visited, vector<int>& group, const pair<int,int> &cut) {
    visited.push_back(u);
    group.push_back(u);
    for (auto &[v, eid] : adj[u]) {
        if (visited.check(v)) continue;
        Edge e = mst[eid];
        if ((e.u == cut.first && e.v == cut.second) || (e.u == cut.second && e.v == cut.first))
            continue;
        collectGroup(v, visited, group, cut);
    }
}
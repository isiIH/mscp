#include <SetCover.h>

SetCover::SetCover() {};

SetCover::SetCover(SCP &scp) : scp(scp) {
    if(GROUP_SEG) {
        g = Group(SEG_TYPE, scp.n, scp.nWX);
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
    if(GROUP_SEG) {
        edges.clear();
        int sharedSubsets;
        for(int i=0; i<rowMap.size() - 1; i++) {
            for(int j=i+1; j<rowMap.size(); j++) {
                sharedSubsets = rowMap[i].countIntersection(rowMap[j].col_covering);
                if(sharedSubsets > 0) {
                    // Add the edges to the group
                    edges.push_back(Edge(rowMap[i].row, rowMap[j].row, sharedSubsets));
                }
            }
            
        }
        g.findGroups(edges, rowMap);
        printf("Groups: %d\n", g.sizeGroups());
    }
    g.printGroups();
    end_time = chrono::high_resolution_clock::now();
    printf("Time Segmentation: %f\n", chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0);

    if(PRINT) {
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

        #pragma omp for schedule(dynamic, 1) nowait
        for(int i=0; i < scp.m-1; i++) {
            vector<Edge> local_edges;
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
            
            if(ignore) {
                local_edges.clear();
            } else {
                #pragma omp critical
                {
                    edges.insert(edges.end(), local_edges.begin(), local_edges.end());
                }
            }
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
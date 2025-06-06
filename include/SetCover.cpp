#include <SetCover.h>

SetCover::SetCover() {};

SetCover::SetCover(SCP &scp) : scp(scp) {
    excludedSets = Set(scp.nWF);
    X = Set(scp.X);
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
    if(PRINT) printf("Time Column Domination: %f\n", chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0);

    // Create map structure
    start_time = chrono::high_resolution_clock::now();
    if(PRINT) printf("Creating Row Map...\n");
    #pragma omp parallel for schedule(dynamic, 1)
    for(int i=0; i<scp.n; i++)
        rowMap[i] = RowCovering(excludedSets, scp.bF, i);

    sort(execution::par, rowMap.begin(), rowMap.end(), [&](RowCovering a, RowCovering b){return a.n_columns < b.n_columns;});
    end_time = chrono::high_resolution_clock::now();
    if(PRINT) printf("Time Create Map: %f\n", chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0);

    printRowMap();
    // Add uniques elements
    start_time = chrono::high_resolution_clock::now();
    rowReduction();
    end_time = chrono::high_resolution_clock::now();
    if(PRINT) printf("Time RowReduction: %f\n", chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0);

    printRowMap();

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
                    break;
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

Set SetCover::unionSets() {
    Set C(scp.nWX);
    for(const int idS : solution) {
        for(int i=0; i<scp.nWX; i++)
            C.S[i] |= scp.bF[idS].S[i];
    }
    return C;
}

bool SetCover::isCovered() {
    Set coveredElements = unionSets();
    
    for (int i = 0; i < scp.nWX; i++) if ((coveredElements.S[i] & X.S[i]) != X.S[i]) {
        return false;
    }
    return true;
}

void SetCover::redundantSets() {
    vector<int> elemCover(scp.n, 0);
    for(int ss : solution) {
        for(int e : scp.F[ss]) {
            if(X.check((e-1))) elemCover[(e-1)]++;
        }
    }

    bool isRedundant;
    int i = 0;
    while(i < solution.size()) {
        isRedundant = true;
        for(int e : scp.F[solution[i]]) {
            if(X.check((e-1)) && elemCover[(e-1)] == 1) {
                isRedundant = false;
                break;
            }
        }

        if(isRedundant) {
            if(CHECK) printf("Redundant subset erased: %d\n", solution[i]);
            for(int e : scp.F[solution[i]]) {
                if(X.check((e-1))) elemCover[(e-1)]--;
            }
            erase(i);
        } else i++;
    }
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
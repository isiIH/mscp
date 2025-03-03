#include <SetCover.h>

SetCover::SetCover() {};

SetCover::SetCover(SCP &scp) : scp(scp) {
    g = Group(scp.m);
    excludedSets = Set(scp.nWF);
    U = Set(scp.X);
    preprocess();
}

void SetCover::preprocess() {
    if(PRINT) {
        cout << "------------------------" << endl;
        cout << "Executing PreSetCover..." << endl;
        cout << "------------------------" << endl;
    }

    // Column Domination
    columnDomination();

    // Create map structure
    for(int i=0; i<scp.n; i++)
        rowMap.push_back(RowCovering(excludedSets, scp.bF, i));

    sort(rowMap.begin(), rowMap.end(), [&](RowCovering a, RowCovering b){return a.n_columns < b.n_columns;});

    // Add uniques elements
    rowReduction();

    printRowMap();

    // // Universe Segmentation
    g.create_groups(excludedSets);
    // cout << "Original Groups: " << g.groups() << endl;
    // g.print();
    // excludedSets.print();

    if(PRINT) {
        cout << "Added " << uniqueSets.size() << " subsets" << endl; 
        cout << "Excluded " << excludedSets.size() << " subsets" << endl;
        cout << "|X| = " << scp.X.size() << endl;
    }
}

void SetCover::rowReduction() {
    int setIndex, p = 0;

    // Check if the sorted rowMap has unique elements
    while(rowMap[p].n_columns == 1) {
        setIndex = rowMap[p].col_covering[0];

        // Add subset of grade 1
        uniqueSets.push_back(setIndex);
        // push_back(setIndex);
        updateRowMap(setIndex);

        excludedSets.push_back(setIndex);

        p++;
    }
}

void SetCover::columnDomination() {
    // sort the subsets in ascending order
    vector<pair<int, int>> indexedSubsets;
    for(int i=0; i<scp.m; i++) {
        indexedSubsets.push_back({i, scp.bF[i].size()});
    }
    sort(indexedSubsets.begin(), indexedSubsets.end(), [](pair<int, int> &a, pair<int, int> &b){
        return a.second < b.second;
    });

    int nIntersect;
    pair<int, int> setA, setB;

    for(int i=0; i < scp.m-1; i++) {
        for(int j=i+1; j < scp.m; j++) {
            setA = indexedSubsets[i];
            setB = indexedSubsets[j];

            nIntersect = scp.bF[setA.first].intersectionLength(scp.bF[setB.first]);

            // If the intersection is the same size of the smallest subset 
            if(nIntersect == setA.second) {
                excludedSets.push_back(setA.first);
                break;
            }
                    
            // if the set is not excluded and has an intersection > 0, add to the graph
            if(nIntersect) 
                g.add_edge(setA.first, setB.first);
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
            cout << "(" << row.row << ") |" << row.n_columns << "| => ";
            for (int index : row.col_covering) {
                cout << index << " ";
            }
            cout << endl;
        }
    }
}
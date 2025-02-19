#include <SetCover.h>

SetCover::SetCover() {};

SetCover::SetCover(SCP &scp) : scp(scp) {
    excludedSets = Set(scp.nWF);
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
    // Group g(scp.m);
    // for(int i=0; i<scp.m-1; i++) {
    //         for(int j=i+1; j<scp.m; j++) {
    //             if(scp.intersectionLength(scp.bF[i], scp.bF[j]) != 0)
    //                 g.add_edge(i, j);
    //         }
    // }

    // g.create_groups();
    // cout << "Original Groups: " << g.groups() << endl;

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

        updateRowMap(setIndex);

        // Add subset of grade 1
        uniqueSets.push_back(setIndex);
        solution.push_back(setIndex);

        p++;
    }
}

void SetCover::columnDomination() {
    // sort the subsets in ascending order
    vector<pair<int, Set>> indexedSubsets;
    for(int i=0; i<scp.m; i++) {
        indexedSubsets.push_back({i, scp.bF[i]});
    }
    sort(indexedSubsets.begin(), indexedSubsets.end(), [](pair<int, Set> &a, pair<int, Set> &b){
        return a.second.size() < b.second.size();
    });

    Set setA, setB;
    for(int i=0; i < scp.m-1; i++) {
        for(int j=i+1; j < scp.m; j++) {
            setA = indexedSubsets[i].second;
            setB = indexedSubsets[j].second;

            // If the intersection is the same size of the smallest subset 
            if(setA.intersectionLength(setB) == setA.size()) {
                // cout << "SetA: ";
                // for(int e : scp.F[indexedSubsets[i].first]) cout << e << " ";
                // cout << endl;
                // cout << "SetB: ";
                // for(int e : scp.F[indexedSubsets[j].first]) cout << e << " ";
                // cout << endl;
                // cout << endl;

                excludedSets.push_back(indexedSubsets[i].first);
                break;
            }
        }
    }
}

void SetCover::updateRowMap(int setIndex) {
    int nRows, aux, aux2;
    for(int e : scp.F[setIndex]) {
        // Update the covered elements
        // scp.X.erase(e-1);

        // Remove all the subset's elements from the rowMap
        nRows = rowMap.size();
        aux = 0;
        for(int i=0; i<nRows; i++) {
            aux2 = i - aux;
            if(scp.bF[setIndex].check(rowMap[aux2].row)) {
                rowMap.erase(rowMap.begin() + aux2);
                aux++;
            }
        }
    }
}

Set SetCover::unionSets() {
    Set C(scp.nWX);
    for(const int idS : solution) for(int i=0; i<scp.nWX; i++) C.S[i] |= scp.bF[idS].S[i];
    return C;
}

bool SetCover::isCovered() {
    Set coveredElements = unionSets();
    
    for (int i = 0; i < scp.nWX; i++) if ((coveredElements.S[i] & scp.X.S[i]) != scp.X.S[i]) {
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
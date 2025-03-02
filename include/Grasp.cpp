#include <Grasp.h>

Grasp::Grasp(SCP &scp): scp(scp) {}

SetCover Grasp::search() {
    // Preprocess
    SetCover initialSol(scp);

    SetCover newSol, bestSol = initialSol;
    Set unionSC;
    int i, col, nRemove;
    vector<int> setsRemoved;
    improve = false;

    // Initial solution
    randSuccintSC(bestSol);

    if(PRINT) cout << "Initial Sol. Cardinality: " << bestSol.size() << endl;

    for(int iter=0; iter< MAX_ITER; iter++){
        // Perturbation
        newSol = bestSol;
        nRemove = rand() % (int)ceil((newSol.size()-newSol.uniqueSets.size()) * RCL) + 1;

        if(PRINT) {
            cout << "--------------------------------------------" << endl;
            cout << "IT: " << (iter+1) << endl;
            cout << nRemove << " subsets deleted" << endl;
        }

        if(CHECK) {
            cout << "SOLUTION = { ";
            for(int i=0; i<newSol.size(); i++) cout << "S" << newSol.solution[i] << " ";
            cout << "}" << endl;
            cout << "DEL = { ";
        }

        for(int i=0; i<nRemove; i++) {
            col = rand()%(newSol.size()-newSol.uniqueSets.size()) + newSol.uniqueSets.size();

            if(CHECK) {
                cout << newSol.solution[col] << " ";
            }
                
            newSol.erase(col);
        }

        if(CHECK) cout << "}" << endl;

        // Update RowMap & U
        unionSC = newSol.unionSets();
        for(RowCovering row : initialSol.rowMap)  {
            if(!unionSC.check(row.row)) {
                newSol.rowMap.push_back(row);
                newSol.U.push_back(row.row);
            }
        }

        // New solution
        randSuccintSC(newSol);

        // Delete redundant subsets
        i=newSol.uniqueSets.size();
        while(i < newSol.size()){
            if(newSol.isCovered(newSol.solution[i])) {
                if(CHECK) cout << "Redundant subset erased: " << newSol.solution[i] << endl;
                newSol.erase(i);
            }
            else i++;
        }

        // Evaluate and Upgrade solution
        if(newSol.size() < bestSol.size()) {
            bestSol = newSol;
            improve = true;
        } else improve = false;

        if(PRINT) {
            cout << endl;
            cout << "Sol. Cardinality: " << newSol.size() << endl;
            cout << "Best Cardinality: " << bestSol.size() << endl;
        }
    }

    return bestSol;
}

void Grasp::randSuccintSC(SetCover &C) {
    if(CHECK) {
        cout << "--------------------------------------------" << endl;
    }
    // function = rand() % 3;
    function = 0;
    set<int> subsets;
    // vector<int> subsets;
    double coverage, bestCoverage;
    int grade, p, bestSet;

    double rand_subset, total = 0;
    vector<pair<int, int>> subsets_coverage;

    if(PRINT) {
        switch(function) {
            case 0: cout << "Using function (1/rowsCovered)" << endl; break;
            case 1: cout << "Using function (1/sqrt(rowsCovered))" << endl; break;
            case 2: cout << "Using function (1/log(1 + rowsCovered))" << endl; break;
            case 3: cout << "Using function (1/rowsCovered²)" << endl; break;
            default: break;
        }
    }

    while( C.U.size() > 0 ) { // Iterate until U is empty
        p = 0;
        grade = C.rowMap[p].n_columns;
        bestCoverage = numeric_limits<double>::max();
        
        // Collect all the element's subsets of grade K
        // subsets = C.rowMap[p].col_covering;
        while(p < C.rowMap.size() && C.rowMap[p].n_columns == grade) {
            for(int ss : C.rowMap[p].col_covering) subsets.insert(ss);
            p++;
        }

        if(CHECK) {
            cout << "grade: " << grade << endl;
            for(int ss : subsets) cout << ss << " ";
            cout << endl;
        }

        // Evaluate each subset with a coverage function
        for(int ss : subsets) {
            coverage = C.U.intersectionLength(scp.bF[ss]);
            switch(function) {
                case 0: coverage = 1/coverage; break;
                case 1: coverage = 1/sqrt(coverage); break;
                case 2: coverage = 1/log(1 + coverage); break;
                case 3: coverage = 1/(coverage * coverage); break;
                default: break;
            }

            if(coverage < bestCoverage) {
                bestCoverage = coverage;
                bestSet = ss;
            }
            if(!improve) {
                total += coverage;
                subsets_coverage.push_back(make_pair(ss, total));
            }
        }
        
        if(!improve && rand() % 25 == 0) {
            if(CHECK) cout << "random set" << endl;
            rand_subset = ((double) rand()) / RAND_MAX;
            for(int i=0; i<subsets_coverage.size(); i++) {
                if(rand_subset <= subsets_coverage[i].second / total) {
                    bestSet = subsets_coverage[i].first;
                    break;
                }
            }
            subsets_coverage.clear();
            total = 0;
        }

        // Erase element covered by the best candidate and add to the solution
        C.push_back(bestSet);

        if(CHECK) {
            cout << "Best Coverage: " << bestCoverage << endl;
            cout << "Pos. Subset: " << bestSet << endl;
            cout << "|U|: " << C.U.size() << endl;
            for(int e : scp.F[bestSet]) {
                cout << e << " ";
            }
            cout << endl;
            C.printRowMap();
        }

        subsets.clear();
    }
}

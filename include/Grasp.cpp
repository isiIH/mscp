#include <Grasp.h>

Grasp::Grasp(SCP &scp): scp(scp) {}

SetCover Grasp::search() {
    // Preprocess
    SetCover initialSol(scp);

    int i;
    Set U(scp.X);
    SetCover newSol, bestSol = initialSol;
    Set unionSC;
    int col;
    int nRemove;
    vector<int> setsRemoved;
    improve = false;

    // Initial solution
    randSuccintSC(U, bestSol);

    improve = true;

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
            cout << "{ ";
            for(int i=0; i<newSol.size(); i++) cout << "S" << newSol.solution[i] << " ";
            cout << "}" << endl;
        }

        for(int i=0; i<nRemove; i++) {
            col = rand()%(newSol.size()-newSol.uniqueSets.size()) + newSol.uniqueSets.size();
            setsRemoved.push_back(newSol.solution[col]);

            if(CHECK) {
                cout << newSol.solution[col] << endl;
            }
            
            
            newSol.solution.erase(newSol.solution.begin() + col);
        }

        // Update RowMap & U
        unionSC = newSol.unionSets();
        for(RowCovering row : initialSol.rowMap)  {
            if(!unionSC.check(row.row)) {
                newSol.rowMap.push_back(row);
                U.push_back(row.row);
            }
        }
        setsRemoved.clear();


        // New solution
        randSuccintSC(U, newSol);


        // Eliminar subsets redundantes (que no agregan elementos nuevos)
        // i=newSol.uniqueSets.size();
        // while(i < newSol.size()){
        //     SetCover sol = newSol;
        //     sol.solution.erase(sol.solution.begin() + i);
        //     if(sol.isCovered()) {
        //         if(CHECK) cout << "Redundant subset erased: " << newSol.solution[i] << endl;
        //         newSol.solution.erase(newSol.solution.begin() + i);
        //     }
        //     else i++;
        // }

        if(newSol.size() < bestSol.size()) {
            bestSol = newSol;
            // improve = true;

            //Penalizar columnas repetidas en la solución anterior
            // for(int ss : new_sol)  {
            //     if(find(sol.begin(), sol.end(), ss) != sol.end()) {
            //         scp.n_columns_colums[ss]++;
            //         scp.worst_columns[ss] = 1.1;
            //         if(CHECK) cout << "subset " << ss << " repeated" << endl;
            //     } else {
            //         scp.worst_columns[ss] = 0.8;
            //         scp.n_columns_colums[ss] = 0;
            //     }
            // }
        }
        // } else improve = false;

        if(PRINT) {
            cout << endl;
            cout << "Sol. Cardinality: " << newSol.size() << endl;
            // printSubsets(new_sol);
            cout << "Best Cardinality: " << bestSol.size() << endl;
        }
    }

    return bestSol;
    
}

void Grasp::randSuccintSC(Set &U, SetCover &C) {
    // scp.function = rand() % 4;
    function = 0;
    set<int> subsets;
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

    while( C.rowMap.size() > 0 ) { // Iterate until rowMap is empty
        p = 0;
        grade = C.rowMap[p].n_columns;
        bestCoverage = numeric_limits<double>::max();
        
        // Collect all the elements'subsets of grade K
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
            coverage = U.intersectionLength(scp.bF[ss]);
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
                subsets_coverage.push_back(make_pair(bestSet, total));
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
        U.substract(scp.bF[bestSet]);
        C.solution.push_back(bestSet);

        // Delete the elements from the best candidate
        C.updateRowMap(bestSet);

        if(CHECK) {
            cout << "Best Coverage: " << bestCoverage << endl;
            cout << "Pos. Subset: " << bestSet << endl;
            cout << "|U|: " << U.size() << endl;
            // scp.printSubset(U);
        }

        subsets.clear();
    }
}

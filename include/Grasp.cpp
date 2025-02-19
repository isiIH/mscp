#include <Grasp.h>

Grasp::Grasp(SCP &scp): scp(scp) {}

SetCover Grasp::search() {
    int i;
    Set U(scp.X);
    SetCover sol(scp);
    SetCover new_sol = sol;
    ulong* unionSC;
    int col;
    int nRemove;
    vector<int> setsRemoved;
    improve = false;


    // Initial solution
    // randSuccintSC(U, sol);

    // if(PRINT) cout << "Initial Sol. Cardinality: " << sol.size() << endl;

    // for(int iter=0; iter< MAX_ITER; iter++){
    //     // Perturbation
    //     new_sol = sol;
    //     nRemove = rand() % (int)ceil((new_sol.size()-scp.uniqueSets.size()) * RCL) + 1;

    //     if(PRINT) {
    //         cout << "--------------------------------------------" << endl;
    //         cout << "IT: " << (iter+1) << endl;
    //         cout << nRemove << " subsets deleted" << endl;
    //     }

    //     if(CHECK) {
    //         cout << "{ ";
    //         for(int i=0; i<new_sol.size(); i++) cout << "S" << new_sol[i] << " ";
    //         cout << "}" << endl;
    //     }

    //     for(int i=0; i<nRemove; i++) {
    //         col = rand()%(new_sol.size()-scp.uniqueSets.size()) + scp.uniqueSets.size();
    //         setsRemoved.push_back(new_sol[col]);
            
    //         new_sol.erase(new_sol.begin() + col);
    //     }

    //     // Update U
    //     unionSC = scp.unionSets(new_sol);
    //     for(int ss : setsRemoved)  {
    //         for(int e : scp.F[ss]) {
    //             if(!checkBit(unionSC, (e-1)))
    //                 setBit64(U, (e-1));
    //         }
    //     }

    //     setsRemoved.clear();

    //     // sort(C.rowMap.begin(), C.rowMap.end(), [&](RowCovering a, RowCovering b){return a.n_columns < b.n_columns;});
        
    //     // Nueva solución
    //     new_sol = randSuccintSC(U, new_sol);

    //     // Eliminar subsets redundantes (que no agregan elementos nuevos)
    //     i=scp.uniqueSets.size();
    //     while(i < new_sol.size()){
    //         vector<int> sol = new_sol;
    //         sol.erase(sol.begin() + i);
    //         if(scp.isCovered(sol)) {
    //             if(CHECK) cout << "Redundant subset erased: " << new_sol[i] << endl;
    //             new_sol.erase(new_sol.begin() + i);
    //         }
    //         else i++;
    //     }

    //     if(new_sol.size() < sol.size()) {
    //         sol = new_sol;
    //         par->improve = true;

    //         //Penalizar columnas repetidas en la solución anterior
    //         // for(int ss : new_sol)  {
    //         //     if(find(sol.begin(), sol.end(), ss) != sol.end()) {
    //         //         scp.n_columns_colums[ss]++;
    //         //         scp.worst_columns[ss] = 1.1;
    //         //         if(CHECK) cout << "subset " << ss << " repeated" << endl;
    //         //     } else {
    //         //         scp.worst_columns[ss] = 0.8;
    //         //         scp.n_columns_colums[ss] = 0;
    //         //     }
    //         // }
    //     } else par->improve = false;

    //     if(PRINT) {
    //         cout << endl;
    //         cout << "Sol. Cardinality: " << new_sol.size() << endl;
    //         // printSubsets(new_sol);
    //         cout << "Best Cardinality: " << sol.size() << endl;
    //     }
    // }

    return sol;
    
}

void Grasp::randSuccintSC(Set U, SetCover &C) {
    // scp.function = rand() % 4;
    function = 0;
    int posSet;
    set<int> subsets;
    double coverage;
    double bestCoverage = numeric_limits<double>::max();
    int grade;
    int p;
    // ulong* Ux = new ulong[scp.nWX];

    double total = 0;
    vector<pair<int, int>> subsets_coverage;
    double rand_subset;

    if(PRINT) {
        switch(function) {
            case 0: cout << "Using function (1/rowsCovered)" << endl; break;
            case 1: cout << "Using function (1/sqrt(rowsCovered))" << endl; break;
            case 2: cout << "Using function (1/log(1 + rowsCovered))" << endl; break;
            case 3: cout << "Using function (1/rowsCovered²)" << endl; break;
            default: break;
        }
    }

    while( U.size() > 0 ) {
        p = last_visited;
        
        while(p < C.rowMap.size() && !U.check(C.rowMap[p].row-1)) p++;
        grade = C.rowMap[p].n_columns;
        // for(int i=0; i<scp.nWX; i++) Ux[i] = 0;
        while(p < C.rowMap.size() && C.rowMap[p].n_columns == grade) {
            if(U.check(C.rowMap[p].row-1))
                for(int ss : C.rowMap[p].col_covering) subsets.insert(ss);
            // setBit64(Ux, scp.elem_pos[C.rowMap[p].row]);
            p++;
        }
        if(CHECK) {
            for(int ss : subsets) cout << ss << " ";
            cout << endl;
        }

        for(int ss : subsets) {
            coverage = U.intersectionLength(scp.bF[ss]);
            switch(function) {
                case 0: coverage = 1/coverage; break;
                case 1: coverage = 1/sqrt(coverage); break;
                case 2: coverage = 1/log(1 + coverage); break;
                case 3: coverage = 1/(coverage * coverage); break;
                default: break;
            }

            // coverage *= scp.worst_columns[ss];

            if(coverage < bestCoverage) {
                bestCoverage = coverage;
                posSet = ss;
            }
            if(!improve) {
                total += coverage;
                subsets_coverage.push_back(make_pair(posSet, total));
            }
        }
        if(!improve && rand() % 25 == 0) {
            if(CHECK) cout << "random set" << endl;
            rand_subset = ((double) rand()) / RAND_MAX;
            for(int i=0; i<subsets_coverage.size(); i++) {
                if(rand_subset <= subsets_coverage[i].second / total) {
                    posSet = subsets_coverage[i].first;
                    break;
                }
            }
            subsets_coverage.clear();
            total = 0;
        }

        U.substract(scp.bF[posSet]);
        C.solution.push_back(posSet);

        // for(int e : scp.F[posSet]) {
        //     C.rowMap.erase(remove_if(C.rowMap.begin(), C.rowMap.end(), [e](const item& mp) {return mp.row == e;}), C.rowMap.end());
        // }

        if(CHECK) {
            cout << "Best Coverage: " << bestCoverage << endl;
            cout << "Pos. Subset: " << posSet << endl;
            cout << "|U|: " << U.size() << endl;
            // scp.printSubset(U);
        }
        bestCoverage = numeric_limits<double>::max();;
        subsets.clear();
    }
}

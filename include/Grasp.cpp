#include <Grasp.h>

Grasp::Grasp(SCP &scp): scp(scp) {}

SetCover Grasp::search() {
    // Preprocess
    SetCover solution(scp);
    
    if(GROUP_SEG) searchPerGroup(solution);
    else {
        auto start_time = chrono::high_resolution_clock::now();
        // Initial solution
        SetCover newSol = solution;
        randSuccintSC(newSol, false);
        
        if(PRINT) printf("Initial Sol. Cardinality: %d\n", newSol.size());
        
        bool improve = true;
        for(int iter=0; iter< MAX_ITER; iter++) {
            if(PRINT) {
                printf("--------------------------------------------\n");
                printf("IT: %d\n", (iter + 1));
            }
            updateSolution(newSol, solution.rowMap, improve);
        }
        
        solution = newSol;
        auto end_time = chrono::high_resolution_clock::now();
        if(PRINT) printf("Time: %f\n", chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0);
    }

    // Add unique sets (grade 1) to the solution
    solution.solution.insert(solution.solution.end(), solution.uniqueSets.begin(), solution.uniqueSets.end());

    return solution;
}

void Grasp::searchPerGroup(SetCover& solution) {
    // Universe Segmentation
    auto start_time = chrono::high_resolution_clock::now();
    g = Group(scp.n, scp.nWX);
    g.findGroups(solution.rowMap);
    g.distributeSubsets(scp.bF, solution.excludedSets, solution.rowMap);
    if(PRINT) printf("Groups: %d\n", g.sizeGroups());
    if(CHECK) g.printGroups();
    auto end_time = chrono::high_resolution_clock::now();
    if(PRINT) printf("Time Segmentation: %f\n", chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0);

    int numGroups = g.sizeGroups();
    vector<vector<int>> groupSolutions(numGroups);

    #pragma omp parallel for default(none) shared(groupSolutions, solution, numGroups, g)
    for(int ng=0; ng < numGroups; ng++) {
        SetCover groupSol = solution;
        vector<int> group = g.groups[ng];
        // if the group is empty or with one subset, skip it
        if(group.size() < 2) {
            groupSolutions[ng] = group;
            continue;
        }

        // Universe of the group
        groupSol.X = g.U[ng];
        groupSol.U = g.U[ng];
        // Erase rows that are not in the group
        vector<RowCovering>& rowMap = groupSol.rowMap;

        rowMap.erase(
            remove_if(rowMap.begin(), rowMap.end(), [&](const RowCovering& row) {
                return !groupSol.X.check(row.row);
            }),
            rowMap.end()
        );
        // Erase the subsets that doesn't belong to the group
        for (RowCovering& row : rowMap) {
            vector<int>& cols = row.col_covering;
            cols.erase(
                remove_if(cols.begin(), cols.end(), [&](int subsetIdx) {
                    return g.subsetToGroup[subsetIdx] != ng;
                }),
                cols.end()
            );
        }

        if(CHECK) {
            printf("|U| = %d\n", groupSol.X.size());
            printf("|rowMap| = %ld\n", groupSol.rowMap.size());
        }

        // Initial solution
        auto start_time = chrono::high_resolution_clock::now();
        SetCover sol = groupSol;
        randSuccintSC(sol, false);
        
        if(PRINT) printf("Initial Sol. Cardinality: %d\n", sol.size());
        
            
        bool improve = true;
        for(int iter=0; iter< MAX_ITER; iter++){
            if(PRINT) {
                printf("--------------------------------------------\n");
                printf("Group %d IT: %d\n", (ng + 1), (iter + 1));
            }
            updateSolution(sol, groupSol.rowMap, improve);

        }
        
        if(PRINT) printf("Group %d (%d, %ld) size: %ld\n", (ng + 1), groupSol.X.size(), group.size(), sol.solution.size());
        
        groupSolutions[ng] = sol.solution;

        auto end_time = chrono::high_resolution_clock::now();

        if(PRINT) printf("Time: %f\n", chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0);

    }

    // copy every group sol to the original solution
    for(vector<int>& groupSol : groupSolutions) {
        solution.solution.insert(solution.solution.end(), groupSol.begin(), groupSol.end());
    }
    
    // Eliminar subconjuntos redundantes
    solution.redundantSets();
}

void Grasp::updateSolution(SetCover& solution, const vector<RowCovering>& rowMap, bool& improve) {
    // Perturbation
    int i, col;
    Set unionSC;
    SetCover newSol = solution;
    int nRemove = rand() % (int)ceil((newSol.size()) * MAX_RM) + 1;

    if(PRINT) printf("%d subsets deleted\n", nRemove);

    if(CHECK) {
        printf("Solution = { ");
        for(i=0; i<newSol.size(); i++) printf("S%d ", newSol.solution[i]);
        printf("}\nDeleted = { ");
    }

    for(i=0; i<nRemove; i++) {
        col = rand()%newSol.size();

        if(CHECK) printf("%d ", newSol.solution[col]);
            
        newSol.erase(col);
    }

    if(CHECK) printf("}\n");

    // Update RowMap & U
    unionSC = newSol.unionSets();
    for(RowCovering row : rowMap)  {
        if(!unionSC.check(row.row)) {
            newSol.rowMap.push_back(row);
            newSol.U.push_back(row.row);
        }
    }

    // New solution
    randSuccintSC(newSol, improve);

    // Delete redundant subsets
    newSol.redundantSets();

    // Evaluate and Upgrade solution
    if(newSol.size() < solution.size()) {
        solution = newSol;
        improve = true;
    } else improve = false;

    if(PRINT) {
        printf("\nSol.Cardinality: %d\nBest Cardinality: %d\n", newSol.size(), solution.size());
    }
}

void Grasp::randSuccintSC(SetCover &C, const bool& improve) {
    if(CHECK) printf("------------------------------------\n");
    int function;
    set<int> subsets;
    double coverage, bestCoverage;
    int grade, p, bestSet;

    double rand_subset, total = 0;
    vector<pair<int, double>> subsets_coverage;

    while( !C.rowMap.empty() ) { // Iterate until rowMap is empty
        function = rand() % 4;
        if(CHECK) {
            switch(function) {
                case 0: printf("Using function (1/rowsCovered)\n"); break;
                case 1: printf("Using function (1/sqrt(rowsCovered))\n");; break;
                case 2: printf("Using function (1/log(1 + rowsCovered))\n"); break;
                case 3: printf("Using function (1/rowsCovered²)\n"); break;
                default: break;
            }
        }
        p = 0;
        grade = C.rowMap[p].n_columns;
        bestCoverage = numeric_limits<double>::max();
        
        // Collect all the element's subsets of grade K
        while(p < C.rowMap.size() && C.rowMap[p].n_columns == grade) {
            for(int ss : C.rowMap[p].col_covering) subsets.insert(ss);
            p++;
        }

        if(CHECK) {
            printf("grade: %d\n", grade);
            for(int ss : subsets) printf("%d ", ss);
            printf("\n");
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

            if(CHECK) printf("Subset %d: %f\n", ss, coverage);

            if(coverage < bestCoverage) {
                bestCoverage = coverage;
                bestSet = ss;
            }
            if(!improve) {
                subsets_coverage.push_back(make_pair(ss, (1 - coverage)));
                total += (1 - coverage);
            }
        }
        
        if(!improve && rand() % 25 == 0) {
            if(CHECK) printf("random set\n");
            rand_subset = ((double) rand()) / RAND_MAX;

            for(int i=0; i<subsets_coverage.size(); i++) {
                subsets_coverage[i].second = subsets_coverage[i].second / total;
                // printf("subset %d coverage: %f\n", subsets_coverage[i].first, subsets_coverage[i].second);
            }
            total = 0;
            for(int i=0; i<subsets_coverage.size(); i++) {
                total += subsets_coverage[i].second;
                subsets_coverage[i].second = total;
            }
            // printf("total prob: %f \n", total);
            for(int i=0; i<subsets_coverage.size(); i++) {
                if(rand_subset <= subsets_coverage[i].second) {
                    bestSet = subsets_coverage[i].first;
                    break;
                }
            }
        }

        // Erase element covered by the best candidate and add to the solution
        C.push_back(bestSet);

        if(CHECK) {
            printf("Best Coverage: %f\n", bestCoverage);
            printf("Pos. Subset: %d\n", bestSet);
            printf("|U|: %d\n", C.U.size());
            for(int e : scp.F[bestSet]) 
                printf("%d ", e);
            printf("\n");
        }

        subsets.clear();
        subsets_coverage.clear();
        total = 0;
    }
}

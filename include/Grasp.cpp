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
        printf("Time: %f\n", chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0);
    }

    // Add unique sets (grade 1) to the solution
    solution.solution.insert(solution.solution.end(), solution.uniqueSets.begin(), solution.uniqueSets.end());

    return solution;
}

void Grasp::searchPerGroup(SetCover& solution) {
    int numGroups = solution.g.sizeGroups();
    vector<SetCover> groupSolutions(numGroups);

    return;

    // for(int ng=0; ng < numGroups; ng++) {
    //     groupSolutions[ng] = solution;
    //     vector<int> group = solution.g.groups[ng];
    //     // if the group is empty or with one subset, skip it
    //     if(group.size() < 2) {
    //         groupSolutions[ng].solution = group;
    //         groupSolutions[ng].U.clear();
    //         continue;
    //     }

    //     // Compute the union of the subsets
    //     groupSolutions[ng].U.clear();
    //     for(int ss : group)
    //         groupSolutions[ng].U.add(scp.bF[ss]);
    //     // Intersect with the original U to delete unique elements
    //     for(int i=0; i<scp.nWX; i++)
    //         groupSolutions[ng].U.S[i] &= solution.U.S[i];

    //     // Compute the new rowMap
    // }

    // // Check if there is no intersection between the groups
    // Set mergeU(groupSolutions[0].U);
    // for(int ng=1; ng < numGroups; ng++)
    //     for(int i=0; i<scp.nWX; i++) mergeU.S[i] = mergeU.S[i] & groupSolutions[ng].U.S[i];

    // // If there is intersection, distribute the elements
    // if(mergeU.size() != 0) {
    //     for(int e=0; e<scp.n; e++) {
    //         if(!mergeU.check(e)) continue;

    //         vector<int> subsets = solution.rowMap[e].col_covering;
    //         // choose one subset randomly
    //         int randSubset = rand() % subsets.size();
    //         int groupId = solution.g.subsetToGroup[subsets[randSubset]];
    //         printf("GroupId %d\n", groupId);
    //         // delete the element from the other groups
    //         for(int ng=0; ng < numGroups; ng++) {
    //             if(ng == groupId) {
                    
    //             };
                
    //             groupSolutions[ng].U.erase(e);

    //         }
    //         // int i = 0;
    //         // for(int ss : subsets) {
    //         //     int groupDel = solution.g.subsetToGroup[ss];
    //         //     printf("Element %d Group %d, subset %d\n", e+1, groupDel, ss);
    //         //     if(groupDel != groupId)  {
    //         //         RowCovering& rc = groupSolutions[groupDel].rowMap[e];
    //         //         rc.col_covering.erase(rc.col_covering.begin() + i);
    //         //         rc.n_columns--;
    //         //         groupSolutions[groupDel].U.erase(e);

    //         //         printf("U size: %d\n", groupSolutions[groupDel].U.size());
    //         //         printf("RowMap size: %ld %d\n", rc.col_covering.size(), rc.n_columns);
    //         //     } else i++;

    //         // }
    //     }
    // }

    // if(1) {
    //     for(int ng=0; ng < numGroups; ng++) {
    //         printf("Group %d (%d) size: %ld\n", (ng + 1), groupSolutions[ng].U.size(), groupSolutions[ng].solution.size());
    //         // for(int ss : groupSolutions[ng].solution) {
    //         //     printf("%d ", ss);
    //         // }
    //         // printf("\n");
    //         // printf("RowMap: ");
    //         // for(RowCovering row : groupSolutions[ng].rowMap) {
    //         //     printf("(%d) |%d| => ", row.row, row.n_columns);
    //         //     for (int index : row.col_covering)
    //         //         printf("%d ", index);
    //         // }
    //         // printf("\n");
    //         // printf("U: ");
    //         // for(int i=0; i<scp.nWX; i++) {
    //         //     if(groupSolutions[ng].U.check(i)) {
    //         //         printf("%d ", i);
    //         //     }
    //         // }
    //         // printf("\n");
    //         printf("U size: %d\n", groupSolutions[ng].U.size());
    //         printf("RowMap size: %ld\n", groupSolutions[ng].rowMap.size());
    //     }
    // }

    // return;

    // for(int ng=0; ng < numGroups; ng++) {
        
    // }
    
    // #pragma omp parallel for default(none) shared(groupSolutions, solution, numGroups)
    // for(int ng=0; ng < numGroups; ng++) {
        
    //     // Define U as the union of the subsets
    //     SetCover groupSol = solution;
    //     vector<int> group = solution.g.groups[ng];
    //     if(group.size() == 1) {
    //         groupSolutions[ng] = group;
    //         continue;
    //     }

    //     groupSol.U.clear();
    //     for(int ss : group) {
    //         if(CHECK) printf("%d ", ss);
    //         groupSol.U.add(scp.bF[ss]);
    //     }
        
    //     // if unique sets were found, intersect with the original U
    //     for(int i=0; i<scp.nWX; i++) {
    //         groupSol.U.S[i] &= solution.U.S[i];
    //     }

    //     if(groupSol.U.size() == 0) continue;
        
    //     // upgrade rowMap
    //     int i=0;
    //     while(i < groupSol.rowMap.size()) {
    //         if(!groupSol.U.check(groupSol.rowMap[i].row)) {
    //             groupSol.rowMap.erase(groupSol.rowMap.begin() + i);
    //         }
    //         else i++;
    //     }
        
    //     int universe = groupSol.U.size();
        
    //     if(CHECK) {
    //         printf("|U| = %d\n", groupSol.U.size());
    //         printf("|rowMap| = %ld\n", groupSol.rowMap.size());
    //     }
        
    //     // Initial solution
    //     auto start_time = chrono::high_resolution_clock::now();
    //     SetCover sol = groupSol;
    //     randSuccintSC(sol, false);
        
    //     if(PRINT) printf("Initial Sol. Cardinality: %d\n", sol.size());
        
    //     bool improve = true;
    //     for(int iter=0; iter< MAX_ITER; iter++){
    //         if(PRINT) {
    //             printf("--------------------------------------------\n");
    //             printf("Group %d IT: %d\n", (ng + 1), (iter + 1));
    //         }
    //         assert(universe == groupSol.rowMap.size());
    //         updateSolution(sol, groupSol.rowMap, improve);

    //     }
        
    //     printf("Group %d (%d, %ld) size: %ld\n", (ng + 1), universe, group.size(), sol.solution.size());
        
    //     groupSolutions[ng] = sol.solution;

    //     auto end_time = chrono::high_resolution_clock::now();

    //     printf("Time: %f\n", chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0);
    // }

    // // copy every group sol to the original solution
    // for(vector<int> groupSol : groupSolutions) {
    //     solution.solution.insert(solution.solution.end(), groupSol.begin(), groupSol.end());
    // }
    
    // // Eliminar subconjuntos redundantes
    // int i=0;
    // while(i < solution.size()){
    //     if(solution.isCovered(solution.U, solution.solution[i])) {
    //         printf("Redundant subset erased: %d\n", solution.solution[i]);
    //         solution.erase(i);
    //     }
    //     else i++;
    // }
    // printf("Best Sol: %d\n", solution.size());
}

void Grasp::updateSolution(SetCover& solution, const vector<RowCovering>& rowMap, bool& improve) {
    // Perturbation
    int i, col;
    Set unionSC;
    SetCover newSol = solution;
    int nRemove = rand() % (int)ceil((newSol.size()) * RCL) + 1;

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
    i=0;
    while(i < newSol.size()){
        if(newSol.isCovered(scp.X, newSol.solution[i])) {
            if(CHECK) printf("Redundant subset erased: %d\n", newSol.solution[i]);
            newSol.erase(i);
        }
        else i++;
    }

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
    // int function = rand() % 3;
    int function = 0;
    set<int> subsets;
    // vector<int> subsets;
    double coverage, bestCoverage;
    int grade, p, bestSet;

    double rand_subset, total = 0;
    vector<pair<int, int>> subsets_coverage;

    if(PRINT) {
        switch(function) {
            case 0: printf("Using function (1/rowsCovered)\n"); break;
            case 1: printf("Using function (1/sqrt(rowsCovered))\n");; break;
            case 2: printf("Using function (1/log(1 + rowsCovered))\n"); break;
            case 3: printf("Using function (1/rowsCovered²)\n"); break;
            default: break;
        }
    }

    while( !C.rowMap.empty() ) { // Iterate until rowMap is empty
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
            if(CHECK) printf("random set");;
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
            printf("Best Coverage: %f\n", bestCoverage);
            printf("Pos. Subset: %d\n", bestSet);
            printf("|U|: %d\n", C.U.size());
            for(int e : scp.F[bestSet]) 
                printf("%d ", e);
            printf("\n");
        }

        subsets.clear();
    }
}

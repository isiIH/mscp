#include <Grasp.h>

Grasp::Grasp(SCP &scp): scp(scp) {}

SetCover Grasp::search() {
    // Preprocess
    SetCover solution(scp);
    
    if(GROUP_SEG) searchPerGroup(solution);
    else {
        // Initial solution
        SetCover newSol = solution;
        randSuccintSC(newSol, false);

        if(PRINT) printf("Initial Sol. Cardinality: %d\n", newSol.size());

        bool improve = true;
        for(int iter=0; iter< MAX_ITER; iter++)
            updateSolution(newSol, solution.rowMap, improve);

        solution = newSol;
    }

    // Add unique sets (grade 1) to the solution
    solution.solution.insert(solution.solution.end(), solution.uniqueSets.begin(), solution.uniqueSets.end());

    return solution;
}

void Grasp::searchPerGroup(SetCover& solution) {
    int numGroups = solution.g.sizeGroups();
    vector<vector<int>> groupSolutions(numGroups);
    
    #pragma omp parallel for default(none) shared(groupSolutions, solution, numGroups)
    for(int ng=0; ng < numGroups; ng++) {
        auto start_time = chrono::high_resolution_clock::now();

        // Define U as the union of the subsets
        SetCover groupSol = solution;
        vector<int> group = solution.g.groups[ng];
        groupSol.U.clear();
        for(int ss : group) {
            if(CHECK) printf("%d ", ss);
            groupSol.U.add(scp.bF[ss]);
        }
        if(CHECK) printf("\n");
        // if unique sets were found, intersect with the original U
        for(int i=0; i<scp.nWX; i++) groupSol.U.S[i] &= solution.U.S[i];

        // upgrade rowMap
        int i=0;
        while(i < groupSol.rowMap.size()) {
            if(!groupSol.U.check(groupSol.rowMap[i].row)) {
                groupSol.rowMap.erase(groupSol.rowMap.begin() + i);
            }
            else i++;
        }

        if(CHECK) {
            printf("|U| = %d\n", groupSol.U.size());
            printf("|rowMap| = %ld\n", groupSol.rowMap.size());
        }

        // Initial solution
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

        printf("Group %d size: %ld\n", (ng + 1), sol.solution.size());

        groupSolutions[ng] = sol.solution;

        auto end_time = chrono::high_resolution_clock::now();

        printf("Time: %f\n", chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0);
    }

    // copy every group sol to the original solution
    for(vector<int> groupSol : groupSolutions)
        solution.solution.insert(solution.solution.end(), groupSol.begin(), groupSol.end());
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
    i=newSol.uniqueSets.size();
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

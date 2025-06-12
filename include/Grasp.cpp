#include <Grasp.h>

Grasp::Grasp(SCP &scp): scp(scp) {}

SetCover Grasp::search() {
    // Preprocess
    SetCover solution(scp);
    
    if(GROUP_SEG) searchPerGroup(solution);
    else {
        auto start_time = chrono::high_resolution_clock::now();
        // Initial solution
        vector<RowCovering> rowMap = solution.rowMap;
        random_device rd;
        mt19937 gen(rd());
        randSuccintSC(solution, false, gen);
        
        if(PRINT) printf("Initial Sol. Cardinality: %d\n", solution.size());
        
        bool improve = true;
        for(int iter=0; iter< MAX_ITER; iter++) {
            if(PRINT) {
                printf("--------------------------------------------\n");
                printf("IT: %d\n", (iter + 1));
            }
            updateSolution(solution, rowMap, improve, gen);
        }
        
        auto end_time = chrono::high_resolution_clock::now();
        if(1) printf("Time: %f\n", chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0);
    }

    // Add unique sets (grade 1) to the solution
    solution.solution.insert(solution.solution.end(), solution.uniqueSets.begin(), solution.uniqueSets.end());

    return solution;
}

void Grasp::searchPerGroup(SetCover& solution) {
    // Universe Segmentation
    auto start_time = chrono::high_resolution_clock::now();
    g = Group(scp.n, scp);
    g.findGroups(solution.rowMap);
    auto a = chrono::high_resolution_clock::now();
    g.distributeSubsets(solution.excludedSets, solution.rowMap);
    auto b = chrono::high_resolution_clock::now();
    if(1) printf("Time distribute: %f\n", chrono::duration_cast<chrono::microseconds>(b - a).count()/1000000.0);
    if(PRINT) printf("Groups: %d\n", g.sizeGroups());
    if(CHECK) g.printGroups();
    auto end_time = chrono::high_resolution_clock::now();
    if(1) printf("Time Segmentation: %f\n", chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0);

    int numGroups = g.sizeGroups();
    vector<vector<int>> groupSolutions(numGroups);

    start_time = chrono::high_resolution_clock::now();
    #pragma omp parallel for
    for(int ng=0; ng < numGroups; ng++) {
        auto start_time = chrono::high_resolution_clock::now();
        SetCover groupSol;
        groupSol.scp = &scp;
        // if the group is empty or with one subset, skip it
        if(g.groups[ng].size() < 2) {
            groupSolutions[ng] = g.groups[ng];
            auto end_time = chrono::high_resolution_clock::now();
            // if(1) printf("Time Preprocces Group %d: %f\n", (ng+1), chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0);
        }
        else {
            // Copy data of each group
            groupSol.X = g.U[ng];
            groupSol.U = g.U[ng];
            groupSol.rowMap = g.groupMap[ng];
            vector<RowCovering> rowMap = groupSol.rowMap;

            auto end_time = chrono::high_resolution_clock::now();

            
            // if(1) printf("Time Preprocces Group %d: %f\n", (ng+1), chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0);
            
            // Initial solution
            start_time = chrono::high_resolution_clock::now();
            random_device rd;
            mt19937 gen(rd() + omp_get_thread_num());
            randSuccintSC(groupSol, false, gen);
            end_time = chrono::high_resolution_clock::now();
            // printf("Group %d (%d, %ld) size: %ld\n", (ng + 1), groupSol.X.size(), g.groups[ng].size(), groupSol.solution.size());
            // if(1) printf("Time randInit Group %d: %f\n", (ng+1), chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0);
            
            if(PRINT) printf("Initial Sol. Cardinality: %d\n", groupSol.size());
                
            start_time = chrono::high_resolution_clock::now();
            if(groupSol.size() != 1) {
                bool improve = true;
                for(int iter=0; iter< MAX_ITER; iter++){
                    if(PRINT) {
                        printf("--------------------------------------------\n");
                        printf("Group %d IT: %d\n", (ng + 1), (iter + 1));
                    }
                    updateSolution(groupSol, rowMap, improve, gen);
                    
                }
            }

            groupSolutions[ng] = groupSol.solution;
            end_time = chrono::high_resolution_clock::now();
            // if(1) printf("Time Group %d: %f\n", (ng+1), chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0);
        }
        
        // if(1) printf("Group %d (%d, %ld) size: %ld\n", (ng + 1), groupSol.X.size(), g.groups[ng].size(), groupSol.solution.size());

    }
    end_time = chrono::high_resolution_clock::now();
    if(1) printf("Time grasp: %f\n", chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0);

    // copy every group sol to the original solution
    for(vector<int>& groupSol : groupSolutions) {
        solution.solution.insert(solution.solution.end(), groupSol.begin(), groupSol.end());
    }
    
    // Delete redundant subsets (MST)
    if(SEG_TYPE) solution.redundantSets();
}

void Grasp::updateSolution(SetCover& solution, const vector<RowCovering>& rowMap, bool& improve, mt19937& gen) {
    // Perturbation
    int i, col;
    vector<int> old_solution = solution.solution;
    int maxRemove = (int)ceil(solution.size() * MAX_RM);
    uniform_int_distribution<int> distRemove(1, maxRemove);
    int nRemove = distRemove(gen);

    if(PRINT) printf("%d subsets deleted\n", nRemove);

    if(CHECK) {
        printf("Solution = { ");
        for(i=0; i<solution.size(); i++) printf("S%d ", solution.solution[i]);
        printf("}\nDeleted = { ");
    }

    for(i=0; i<nRemove; i++) {
        uniform_int_distribution<int> distCol(0, solution.size() - 1);
        col = distCol(gen);

        if(CHECK) printf("%d ", solution.solution[col]);
            
        solution.erase(col);
    }

    if(CHECK) printf("}\n");

    // Update RowMap & U
    Set unionSC(solution.U.nW);
    for(int ss : solution.solution) {
        unionSC.add(solution.scp->bF[ss]);
    }
    for(const RowCovering& row : rowMap)  {
        if(!unionSC.check(row.row)) {
            solution.rowMap.push_back(row);
            solution.U.push_back(row.row);
        }
    }

    // New solution
    randSuccintSC(solution, improve, gen);

    // Delete redundant subsets
    solution.redundantSets();

    vector<int> new_solution = solution.solution;

    // Evaluate and Upgrade solution
    if(new_solution.size() > old_solution.size()) {
        solution.solution = old_solution;
        improve = false;
    } else improve = true;


    if(PRINT) {
        printf("\nSol.Cardinality: %d\nBest Cardinality: %d\n", new_solution.size(), solution.solution.size());
    }
}

void Grasp::randSuccintSC(SetCover &C, const bool& improve, mt19937& gen) {
    if(CHECK) printf("------------------------------------\n");
    int function;
    set<int> subsets;
    double coverage, bestCoverage;
    int grade, p, bestSet;

    double rand_subset, total = 0;
    vector<pair<int, double>> subsets_coverage;
    uniform_int_distribution<int> distFunction(0, 3);
    uniform_int_distribution<int> distImprove(0, 24);

    while( !C.rowMap.empty() ) { // Iterate until rowMap is empty
        function = distFunction(gen);
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
            coverage = C.U.intersectionLength(C.scp->bF[ss]);
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
        
        if(!improve && distImprove(gen) == 0) {
            if(CHECK) printf("random set\n");
            rand_subset = generate_canonical<double, 10>(gen);

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
            for(int e : C.scp->F[bestSet]) 
                printf("%d ", e);
            printf("\n");
        }

        subsets.clear();
        subsets_coverage.clear();
        total = 0;
    }
}

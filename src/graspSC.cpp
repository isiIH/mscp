#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include "../include/BasicCDS.h"
#include <cmath>
#include <set>
#include <map>
#include <unordered_map>
#include <numeric>
#include <cassert>

#define PRINT 1
#define CHECK 0

#include "../include/Group.h"
#include "../include/SCP.h"

using namespace std;
using namespace cds;


// Parámetros
#define RCL 0.7
#define MAX_ITER 100

// Structure with all globals parameters program
typedef struct {
    ulong sizeF, sizeNF;

    int last_visited;

    vector<int> greedy_sol;
    vector<int> aprox_sol;

    int function;
    bool improve;
    vector<float> worst_columns;
    vector<int> rep_colums;

} ParProg;

ParProg* par;
SCP prob;

void preprocess();

void greedy();

vector<int> graspSC();
vector<int> randSuccintSC(ulong* U, vector<int> init_sol);

int main(int argc, char** argv) {

    if(argc !=3){
		cout << "./opt <filename> <seed>" << endl;
		exit(EXIT_FAILURE);
	}

    // srand(atoi(argv[2]));
    srand(time(0));

    par = new ParProg();

    prob.readFile(argv[1]);
    auto start_time = chrono::high_resolution_clock::now();
    prob.analyzeF();
    auto end_time = chrono::high_resolution_clock::now();
    auto dur_analyze = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();

    if(PRINT) cout  << "X: " << prob.n << " | F: " << prob.m << endl;

    par->sizeF = prob.m*sizeof(ulong)*prob.n;
    par->sizeNF = prob.m*sizeof(ulong)*prob.nWX;

	if(PRINT) {
        cout << "nWX = " << prob.nWX << endl;
        cout << " size for F[] = " << par->sizeF/(1024.0*1024.0) << " MiB" << endl;
        cout << " size for nF[] = " << par->sizeNF/(1024.0*1024.0) << " MiB" << endl;
    }

    if(CHECK) {
        for(vector<int> set : prob.F) {
            for(int val : set) {
                cout << val << " ";
            }
            cout << endl;
        }

        prob.printSubsets(prob.bF);
        prob.printSubset(prob.X);
    }

    //GREEDY
    start_time = chrono::high_resolution_clock::now();
    greedy();
    end_time = chrono::high_resolution_clock::now();
    auto dur_greedyExh = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();
    dur_greedyExh += dur_analyze;

    //GRASP
    double dur_apr;
    vector<int> sol;
    int best_card = 9999999;
    for(int i=0; i<1; i++) {
        start_time = chrono::high_resolution_clock::now();
        sol = graspSC();
        end_time = chrono::high_resolution_clock::now();
        auto time = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();
        if(sol.size() < best_card) {
            par->aprox_sol = sol;
            best_card = sol.size();
            dur_apr = time;
        }
    }
    dur_apr += dur_analyze;

    if(CHECK) {
        cout << "SOL: { ";
        for(int ss : par->aprox_sol) {
            cout << ss << " ";
        }
        cout << "}" << endl;
    }
    if(PRINT) {
        cout << "------------------------" << endl;
        cout << "Greedy Cardinality: " << par->greedy_sol.size() << endl;
        cout << "Time [s]: " << dur_greedyExh/1000000.0 << endl;
        cout << "GraspSC Cardinality: " << par->aprox_sol.size() << endl;
        cout << "Time [s]: " << dur_apr/1000000.0 << endl;
    }

    assert(prob.isCovered(par->aprox_sol) && "Solución inválida");

    cout << argv[1] << " " << prob.n << " " << prob.m << " " << dur_greedyExh/1000000.0 << " " << par->greedy_sol.size() << " " << dur_apr/1000000.0 << " " << par->aprox_sol.size() << " " << endl;

    return 0;
}

void greedy() {
    int i;
    ulong* U = new ulong[prob.nWX];
    for(i=0; i<prob.nWX; i++) U[i] = prob.X[i];
    vector<int> C;
    int maxLengthSS = 0;
    int lengthSS;
    int posSet;

    map<int, ulong*> subsets;
    for (i=0; i<prob.bF.size(); i++) subsets[i] = prob.bF[i];

    while( prob.countSet(U, prob.nWX) > 0 ) {

        for(pair<int, ulong*> ss_pos : subsets){
            lengthSS = prob.intersectionLength(U, ss_pos.second);
            if(lengthSS > maxLengthSS) {
                maxLengthSS = lengthSS;
                posSet = ss_pos.first;
            }
        }

        for(i=0; i<prob.nWX; i++) U[i] = U[i] & ~subsets[posSet][i];
        C.push_back(posSet);
        subsets.erase(posSet);

        maxLengthSS = 0;
    }

    par->greedy_sol = C;
}

vector<int> graspSC() {
    // Lista de elementos ordenados por grado
    preprocess();

    int i;
    ulong* U = new ulong[prob.nWX];
    for(i=0; i<prob.nWX; i++) U[i] = prob.X[i];
    vector<int> sol, new_sol;
    par->worst_columns.assign(prob.bF.size(), 1);
    par->rep_colums.assign(prob.bF.size(), 0);
    ulong* unionSC;
    int col;
    int nRemove;
    vector<int> setsRemoved;
    par->improve = false;


    //Solución inicial
    sol = randSuccintSC(U, prob.unique_elements);

    if(PRINT) cout << "Initial Sol. Cardinality: " << sol.size() << endl;

    for(int iter=0; iter< MAX_ITER; iter++){
        //Perturbación
        new_sol = sol;
        nRemove = rand() % (int)ceil((new_sol.size()-prob.unique_elements.size()) * RCL) + 1;

        if(PRINT) {
            cout << "--------------------------------------------" << endl;
            cout << "IT: " << (iter+1) << endl;
            cout << nRemove << " subsets deleted" << endl;
        }

        if(CHECK) {
            cout << "{ ";
            for(int i=0; i<new_sol.size(); i++) cout << "S" << new_sol[i] << " ";
            cout << "}" << endl;
        }

        for(int i=0; i<nRemove; i++) {
            col = rand()%(new_sol.size()-prob.unique_elements.size()) + prob.unique_elements.size();
            setsRemoved.push_back(new_sol[col]);
            
            new_sol.erase(new_sol.begin() + col);
        }

        // Update U
        unionSC = prob.unionSets(new_sol);
        for(int ss : setsRemoved)  {
            for(int e : prob.F[ss]) {
                if(!checkBit(unionSC, (e-1)))
                    setBit64(U, (e-1));
            }
        }

        setsRemoved.clear();

        // sort(prob.mp.begin(), prob.mp.end(), [&](item a, item b){return a.rep < b.rep;});
        
        // Nueva solución
        new_sol = randSuccintSC(U, new_sol);

        // Eliminar subsets redundantes (que no agregan elementos nuevos)
        i=prob.unique_elements.size();
        while(i < new_sol.size()){
            vector<int> sol = new_sol;
            sol.erase(sol.begin() + i);
            if(prob.isCovered(sol)) {
                if(CHECK) cout << "Redundant subset erased: " << new_sol[i] << endl;
                new_sol.erase(new_sol.begin() + i);
            }
            else i++;
        }

        if(new_sol.size() < sol.size()) {
            sol = new_sol;
            par->improve = true;

            //Penalizar columnas repetidas en la solución anterior
            // for(int ss : new_sol)  {
            //     if(find(sol.begin(), sol.end(), ss) != sol.end()) {
            //         prob.rep_colums[ss]++;
            //         prob.worst_columns[ss] = 1.1;
            //         if(CHECK) cout << "subset " << ss << " repeated" << endl;
            //     } else {
            //         prob.worst_columns[ss] = 0.8;
            //         prob.rep_colums[ss] = 0;
            //     }
            // }
        } else par->improve = false;

        if(PRINT) {
            cout << endl;
            cout << "Sol. Cardinality: " << new_sol.size() << endl;
            // printSubsets(new_sol);
            cout << "Best Cardinality: " << sol.size() << endl;
        }
    }

    return sol;
    
}

vector<int> randSuccintSC(ulong* U, vector<int> init_sol) {
    // prob.function = rand() % 4;
    par->function = 0;
    vector<int> C = init_sol;
    int posSet;
    set<int> subsets;
    double coverage;
    double bestCoverage = numeric_limits<double>::max();;
    int grade;
    int p;
    // ulong* Ux = new ulong[prob.nWX];

    double total = 0;
    vector<pair<int, int>> subsets_coverage;
    double rand_subset;

    if(PRINT) {
        switch(par->function) {
            case 0: cout << "Using function (1/rowsCovered)" << endl; break;
            case 1: cout << "Using function (1/sqrt(rowsCovered))" << endl; break;
            case 2: cout << "Using function (1/log(1 + rowsCovered))" << endl; break;
            case 3: cout << "Using function (1/rowsCovered²)" << endl; break;
            default: break;
        }
    }

    while( prob.countSet(U, prob.nWX) > 0 ) {
        p = par->last_visited;
        
        while(p < prob.mp.size() && !checkBit(U, (prob.mp[p].value-1))) p++;
        grade = prob.mp[p].rep;
        // for(int i=0; i<prob.nWX; i++) Ux[i] = 0;
        while(p < prob.mp.size() && prob.mp[p].rep == grade) {
            if(checkBit(U, (prob.mp[p].value-1) ))
                for(int ss : prob.mp[p].subSets) subsets.insert(ss);
            // setBit64(Ux, prob.elem_pos[prob.mp[p].value]);
            p++;
        }
        if(CHECK) {
            for(int ss : subsets) cout << ss << " ";
            cout << endl;
        }

        for(int ss : subsets) {
            coverage = prob.intersectionLength(U, prob.bF[ss]);
            switch(par->function) {
                case 0: coverage = 1/coverage; break;
                case 1: coverage = 1/sqrt(coverage); break;
                case 2: coverage = 1/log(1 + coverage); break;
                case 3: coverage = 1/(coverage * coverage); break;
                default: break;
            }

            // coverage *= prob.worst_columns[ss];

            if(coverage < bestCoverage) {
                bestCoverage = coverage;
                posSet = ss;
            }
            if(!par->improve) {
                total += coverage;
                subsets_coverage.push_back(make_pair(posSet, total));
            }
        }
        if(!par->improve && rand() % 25 == 0) {
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

        for(int i=0; i<prob.nWX; i++) U[i] = U[i] & ~prob.bF[posSet][i];

        C.push_back(posSet);

        // for(int e : prob.F[posSet]) {
        //     prob.mp.erase(remove_if(prob.mp.begin(), prob.mp.end(), [e](const item& mp) {return mp.value == e;}), prob.mp.end());
        // }

        if(CHECK) {
            cout << "Best Coverage: " << bestCoverage << endl;
            cout << "Pos. Subset: " << posSet << endl;
            cout << "|U|: " << prob.countSet(U, prob.nWX) << endl;
            // prob.printSubset(U);
        }
        bestCoverage = numeric_limits<double>::max();;
        subsets.clear();
    }

    return C;
}

void preprocess() {
    if(PRINT) {
        cout << "------------------------" << endl;
        cout << "Executing PreSetCover..." << endl;
        cout << "------------------------" << endl;
    }

    // Universe Segmentation
    Graph graph;
    for(int i=0; i<prob.bF.size()-1; i++) {
            for(int j=i+1; j<prob.bF.size(); j++) {
                if(prob.intersectionLength(prob.bF[i], prob.bF[j]) != 0)
                    graph.add_edge(i, j);
            }
    }

    // graph.print();
    Group g(graph, prob.bF.size());
    g.create_groups();
    // g.print();
    cout << "Original Groups: " << g.groups() << endl;

    // Check for a segmentation ignoring a subset
    // for(int k=0; k<prob.bF.size(); k++) {
    //     Group gr(graph, prob.bF.size(), k);
    //     gr.create_groups();
    //     if(gr.groups() > 1)
    //         cout << "Ignoring subset " << k << ": " << gr.groups() << " groups" << endl;
    //     // gr.print();
    // }

    // Column domination
    vector<bool> visited(prob.bF.size(), false);
    int minSubset;
    for(int u=0; u<prob.bF.size(); u++) {
        if(!visited[u]) {
            for(int v : graph.adj_list[u]) {
                if(!visited[v] && prob.intersectionLength(prob.bF[u], prob.bF[v]) == prob.countSet(prob.bF[u], prob.nWX)) {
                    // prob.printSubset(prob.bF[u]);
                    // prob.printSubset(prob.bF[v]);
                    visited[u] = true;
                    setBit64(prob.excluded_subsets, u);
                    break;
                }
            }
        }
    }

    // Create map structure
    int pos = 0;
    prob.mp = vector<item>(prob.n);
    for(pair<int, vector<int>> values : prob.inSet){
        // prob.elem_pos[values.first] = pos;
        prob.mp[pos].value = values.first;
        for(int ss : values.second)
            values.second.erase(remove_if(values.second.begin(), values.second.end(), [](const int ss) {return checkBit(prob.excluded_subsets, ss);}), values.second.end());
        prob.mp[pos].subSets = values.second;
        prob.mp[pos].rep = prob.mp[pos].subSets.size();
        pos++;
    }
    sort(prob.mp.begin(), prob.mp.end(), [&](item a, item b){return a.rep < b.rep;});

    if(CHECK) {
        for(item mp_item : prob.mp) {
            cout << "(" << mp_item.value << ") |" << mp_item.rep << "| => ";
            for (int index : mp_item.subSets) {
                cout << index << " ";
            }
            cout << endl;
        }
    }

    // Add uniques elements - Row Reduction
    int setIndex;
    ulong* S;
    int p = 0;
    while(prob.mp[p].rep == 1) {
        setIndex = prob.mp[p].subSets[0];

        // Update the covered elements
        for(int e : prob.F[setIndex]) {
            // prob.mp.erase(remove_if(prob.mp.begin(), prob.mp.end(), [e](const item& mp) {return mp.value == e;}), prob.mp.end());
            cleanBit64(prob.X, (e-1));
        }

        prob.unique_elements.push_back(setIndex);
        // setBit64(prob.excluded_subsets, setIndex);

        par->last_visited = p;
        p++;
    }

    if(PRINT) {
        cout << "Added " << prob.unique_elements.size() << " subsets" << endl; 
        cout << "Excluded " << prob.countSet(prob.excluded_subsets, prob.nWF) << " subsets" << endl;
        cout << "|X| = " << prob.countSet(prob.X, prob.nWX) << endl;
        cout << "|F| = " << prob.bF.size() << endl;
    }
}
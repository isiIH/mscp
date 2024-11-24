#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <chrono>
#include "../include/BasicCDS.h"
#include <cmath>
#include <set>
#include <map>
#include <unordered_map>
#include <numeric>
#include <cassert>

using namespace std;
using namespace cds;

#define PRINT 1
#define CHECK 0

// Parámetros
#define RCL 0.7
#define TOLERANCE 50
#define ITER_FACTOR 0.1
#define MAX_ITER 300

typedef struct{
	int value;
	int rep;
	vector<int> subSets;
} item;

// Structure with all globals parameters program
typedef struct {
	ulong* X;
	vector<vector<int>> F;
    vector<ulong*> bF;

    // set<int> chi;
    // map<int,int> elem_pos;
    unordered_map<int, vector<int>> inSet;
    vector<item> mp;

    vector<int> unique_elements;
    vector<int> greedy_sol;
    vector<int> aprox_sol;

    ulong sizeF, sizeNF;
	ulong n, m, nWX;
    int nt;

    int function;
    vector<float> costFunc = {1,1,1,1};
    bool improve;
    vector<float> worst_columns;
    vector<int> rep_colums;

} ParProg;

ParProg* par;

void readFile(string filename);
void readFileScp(string filename);
void readFilePartition(string filename);
void analyzeF();
void createMap();
void preprocess();

void greedy();

int chooseFunction();
double jaccard(const ulong* A, const ulong* B);
vector<int> graspSC();
vector<int> randSuccintSC(ulong* U, vector<int> init_sol);

bool isCovered(vector<int> S);
ulong* unionSets(const vector<int> &S);
int countSet(const ulong* S);
int intersectionLength(const ulong* A, const ulong* B);

void printSubset(const ulong *S);
void printSubsets(const vector<ulong*> &C);

int main(int argc, char** argv) {

    if(argc !=3){
		cout << "./opt <filename> <seed>" << endl;
		exit(EXIT_FAILURE);
	}

    // srand(atoi(argv[2]));
    srand(time(0));

    par = new ParProg();

    readFile(argv[1]);
    auto start_time = chrono::high_resolution_clock::now();
    analyzeF();
    auto end_time = chrono::high_resolution_clock::now();
    auto dur_analyze = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();

    if(PRINT) cout  << "X: " << par->n << " | F: " << par->m << endl;

    par->sizeF = par->m*sizeof(ulong)*par->n;
    par->sizeNF = par->m*sizeof(ulong)*par->nWX;

	if(PRINT) {
        cout << "nWX = " << par->nWX << endl;
        cout << " size for F[] = " << par->sizeF/(1024.0*1024.0) << " MiB" << endl;
        cout << " size for nF[] = " << par->sizeNF/(1024.0*1024.0) << " MiB" << endl;
    }

    if(CHECK) {
        for(vector<int> set : par->F) {
            for(int val : set) {
                cout << val << " ";
            }
            cout << endl;
        }

        printSubsets(par->bF);
        printSubset(par->X);
    }

    //GREEDY
    start_time = chrono::high_resolution_clock::now();
    greedy();
    end_time = chrono::high_resolution_clock::now();
    auto dur_greedyExh = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();
    dur_greedyExh += dur_analyze;

    //GRASP
    int dur_apr;
    vector<int> sol;
    int best_card = 9999999;
    for(int i=0; i<1; i++) {
        start_time = chrono::high_resolution_clock::now();
        sol = graspSC();
        end_time = chrono::high_resolution_clock::now();
        int time = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();
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

    assert(isCovered(par->aprox_sol) && "Solución inválida");

    cout << argv[1] << " " << par->n << " " << par->m << " " << dur_greedyExh/1000000.0 << " " << par->greedy_sol.size() << " " << dur_apr/1000000.0 << " " << par->aprox_sol.size() << " " << endl;
    // cout << par->costFunc[0] << " " << par->costFunc[1] << " " << par->costFunc[2] << " " << par->costFunc[3] << endl;

    return 0;
}

void readFile(string filename) {
    if (filename.substr(0,3) == "scp") readFileScp(filename);
    else readFilePartition(filename);
}

void readFileScp(string filename) {
    if(PRINT) cout << "Reading file " << filename << "..." << endl;
    string nametxt = "test/" + filename;
    ifstream file(nametxt.c_str());
    if(file.fail()){
        cout << "File not found!" << endl;
        exit(EXIT_FAILURE);
    }
    string line,item;
    int i;

    //m & n
	getline(file>>std::ws,line);
    istringstream ss(line);
    ss >> (par->n) >> (par->m);

    //Costs
    i = 0;
    while(i < par->m)
    {
        getline(file>>std::ws,line);
        istringstream iss(line);
        while (getline(iss>>std::ws, item, ' ')){i++;}
    }

    //Sets
    int numCover;
    int j;
    par->F.resize(par->m);
    for(i=0; i<par->n; i++) {
        getline(file>>std::ws,line);
        numCover = stoi(line);

        j = 0;
        while(j < numCover){
            getline(file>>std::ws,line);
            istringstream iss(line);
            while (getline(iss, item, ' ')) {
                par->F[stoi(item)-1].push_back(i+1);
                j++;
            }
        }
    }
    file.close();
}

void readFilePartition(string filename) {
    if(PRINT) cout << "Reading file " << filename << "..." << endl;
    string nametxt = "test/" + filename;
    ifstream file(nametxt.c_str());
    if(file.fail()){
        cout << "File not found!" << endl;
        exit(EXIT_FAILURE);
    }
    string line,item;

    //m & n
	getline(file>>std::ws,line);
    istringstream ss(line);
    ss >> (par->n) >> (par->m);

    //Sets
    vector<int> sub;
    for (int i = 0; i < par->m; i++) {
        getline(file>>std::ws,line);
        istringstream ss(line);
        getline(ss>>std::ws, item, ' ');
        getline(ss>>std::ws, item, ' ');

        while (getline(ss>>std::ws, item, ' ')) {
            sub.push_back(stoi(item));
        }
        (par->F).push_back(sub);
        sub.clear();
    }
    file.close();
}

void analyzeF() {
    // for( int i=0; i<par->F.size(); i++ ) {
    //     for( int e : par->F[i] ) {
    //         par->chi.insert(e);
    //         par->inSet[e].push_back(i);
    //     }
    // }

    // par->n = par->chi.size();

    // par->nWX = (par->n)/(sizeof(ulong)*8);
    // if ((par->n)%(sizeof(ulong)*8)>0) par->nWX++;
    // par->X = new ulong[par->nWX];
    // fill(par->X, par->X + par->nWX, 0);
    // for(int i=0; i<par->n; i++) setBit64(par->X, i);

    // ulong *bset;
    // for( int i=0; i<par->F.size(); i++ ) {
    //     bset = new ulong[par->nWX];
    //     fill(bset, bset + par->nWX, 0);

    //     for( int e : par->F[i] ) {
    //         setBit64(bset, par->elem_pos[e]);
    //     }

    //     par->bF.push_back(bset);
    // }

    // if(CHECK) {
    //     cout << "Universe elements = " << endl;
    //     for( pair<int, int> values : par->elem_pos ) if(getBit64(par->X, values.second)) cout << values.first << " ";
    //     cout << endl;
    //     cout << "X = " << countSet(par->X) << endl;
    //     cout << "n = " << par->n << endl;
    //     cout << "F = " << par->bF.size() << endl;
    //     cout << "m = " << par->m << endl;
    // }
    par->nWX = (par->n)/(sizeof(ulong)*8);
    if ((par->n)%(sizeof(ulong)*8)>0) par->nWX++;
    
    par->X = new ulong[par->nWX];
    fill(par->X, par->X + par->nWX, 0);

    // par->mp = vector<item>(par->n);
    ulong *bset;
    for(int i=0; i<par->F.size(); i++){
        bset = new ulong[par->nWX];
        fill(bset, bset + par->nWX, 0);

        for(int e : par->F[i]) {
            setBit64(par->X, (e-1));
            par->inSet[e].push_back(i);
            // par->mp[(e-1)].value = e;
            // par->mp[(e-1)].subSets.push_back(i);

            setBit64(bset, (e-1));
        }

        par->bF.push_back(bset);
    }

    // for(int i=0; i<par->mp.size(); i++) par->mp[i].rep = par->mp[i].subSets.size();

    // sort(par->mp.begin(), par->mp.end(), [&](item a, item b){return a.rep < b.rep;});

    if(CHECK) {
        cout << "X = " << countSet(par->X) << endl;
        cout << "F = " << par->bF.size() << endl;
    }
}

void greedy() {
    int i;
    ulong* U = new ulong[par->nWX];
    for(i=0; i<par->nWX; i++) U[i] = par->X[i];
    vector<int> C = par->unique_elements;
    int maxLengthSS = 0;
    int lengthSS;
    int posSet;

    map<int, ulong*> subsets;
    for (i=0; i<par->bF.size(); i++) subsets[i] = par->bF[i];

    while( countSet(U) > 0 ) {

        for(pair<int, ulong*> ss_pos : subsets){
            lengthSS = intersectionLength(U, ss_pos.second);
            if(lengthSS > maxLengthSS) {
                maxLengthSS = lengthSS;
                posSet = ss_pos.first;
            }
        }

        for(i=0; i<par->nWX; i++) U[i] = U[i] & ~subsets[posSet][i];
        C.push_back(posSet);
        subsets.erase(posSet);

        maxLengthSS = 0;
    }

    par->greedy_sol = C;
}

double jaccard(const ulong* A, const ulong* B) {
    double cont = 0.0;
    for(int i = 0; i < par->nWX; i++) {
        int inter_set = __builtin_popcountl(A[i] & B[i]);
        int union_set = __builtin_popcountl(A[i] | B[i]);
        
        if (union_set > 0) {
            cont += (double)(inter_set) / union_set;
        }
    }
    return cont;
}

vector<int> graspSC() {
    // Lista de elementos ordenados por grado
    createMap();

    int i;
    ulong* U = new ulong[par->nWX];
    for(i=0; i<par->nWX; i++) U[i] = par->X[i];
    vector<int> sol, new_sol;
    par->worst_columns.assign(par->bF.size(), 1);
    par->rep_colums.assign(par->bF.size(), 0);
    ulong* unionSC;
    int col;
    int nRemove;
    vector<int> setsRemoved;
    par->improve = false;
    int iter = 0;
    int tol = 0;


    //Solución inicial
    sol = randSuccintSC(U, par->unique_elements);

    if(PRINT) cout << "Initial Sol. Cardinality: " << sol.size() << endl;

    while((par->improve || tol < (int)(TOLERANCE + sol.size() * ITER_FACTOR)) && iter < MAX_ITER){
        //Perturbación
        new_sol = sol;
        nRemove = rand() % (int)ceil((new_sol.size()-par->unique_elements.size()) * RCL) + 1;

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
            col = rand()%(new_sol.size()-par->unique_elements.size()) + par->unique_elements.size();
            setsRemoved.push_back(new_sol[col]);
            new_sol.erase(new_sol.begin() + col);
        }

        // Actualizar U y map
        unionSC = unionSets(new_sol);
        for(int ss : setsRemoved)  {
            for(int e : par->F[ss]) {
                if(checkBit(U, (e-1)) == 0 && checkBit(unionSC, (e-1)) == 0) {
                    item it_map;
                    it_map.value = e;
                    it_map.subSets = par->inSet[e];
                    it_map.rep = it_map.subSets.size();
                    par->mp.push_back(it_map);
                    setBit64(U, (e-1));
                }

            }
        }

        setsRemoved.clear();

        sort(par->mp.begin(), par->mp.end(), [&](item a, item b){return a.rep < b.rep;});

        if(CHECK) {

            for(item mp_item : par->mp) {
                cout << "(" << mp_item.value << ") |" << mp_item.rep << "| => ";
                for (int index : mp_item.subSets) {
                    cout << index << " ";
                }
                cout << endl;
            }
        }
        
        // Nueva solución
        new_sol = randSuccintSC(U, new_sol);

        // Eliminar subsets redundantes (que no agregan elementos nuevos)
        i=par->unique_elements.size();
        while(i < new_sol.size()){
            vector<int> sol = new_sol;
            sol.erase(sol.begin() + i);
            if(isCovered(sol)) {
                if(CHECK) cout << "Redundant subset erased: " << new_sol[i] << endl;
                new_sol.erase(new_sol.begin() + i);
            }
            else i++;
        }

        if(new_sol.size() < sol.size()) {
            sol = new_sol;
            par->improve = true;
            tol = 0;
            // par->costFunc[par->function] = min(2.0, par->costFunc[par->function] + 0.2);

            //Penalizar columnas repetidas en la solución anterior
            // for(int ss : new_sol)  {
            //     if(find(sol.begin(), sol.end(), ss) == sol.end()) {
            //         par->rep_colums[ss]++;
            //         par->worst_columns[ss] = 1.1;
            //         if(CHECK) cout << "subset " << ss << " repeated" << endl;
            //     } else {
            //         par->worst_columns[ss] = 1;
            //         par->rep_colums[ss] = 0;
            //     }
            // }
        } else {
            tol++;
            par->improve = false;
        }


        if(PRINT) {
            cout << endl;
            cout << "Sol. Cardinality: " << new_sol.size() << endl;
            // printSubsets(new_sol);
            cout << "Best Cardinality: " << sol.size() << endl;
        }
        iter++;
    }

    return sol;
    
}

vector<int> randSuccintSC(ulong* U, vector<int> init_sol) {
    par->function = rand() % 4;
    // par->function = chooseFunction();
    // par->function = 0;
    vector<int> C = init_sol;
    int posSet;
    set<int> subsets;
    vector<pair<int, int>> subsets_coverage;
    double coverage;
    double bestCoverage = 0.0;
    int grade;
    int p;
    // ulong* Ux = new ulong[par->nWX];

    double total;
    double rand_subset;

    if(PRINT) {
        switch(par->function) {
            case 0: cout << "Using function (1/rowsCovered)" << endl; break;
            case 1: cout << "Using function (1/sqrt(rowsCovered))" << endl; break;
            case 2: cout << "Using function (1/log(1 + rowsCovered))" << endl; break;
            case 3: cout << "Using function (1/rowsCovered²)" << endl; break;
        }
    }

    while( countSet(U) > 0 ) {
        total = 0;
        grade = par->mp[0].rep;
        p = 0;
        // for(int i=0; i<par->nWX; i++) Ux[i] = 0;
        while(p < par->mp.size() && par->mp[p].rep == grade) {
            for(int ss : par->mp[p].subSets) subsets.insert(ss);
            // setBit64(Ux, par->elem_pos[par->mp[p].value]);
            p++;
        }
        if(CHECK) {
            for(int ss : subsets) cout << ss << " ";
            cout << endl;
        }
        for(int ss : subsets) {
            coverage = intersectionLength(U, par->bF[ss]);
            switch(par->function) {
                case 1: coverage = sqrt(coverage); break;
                case 2: coverage = log(1 + coverage); break;
                case 3: coverage = (coverage * coverage); break;
                default: break;
            }

            coverage /= par->worst_columns[ss];

            if(coverage > bestCoverage) {
                bestCoverage = coverage;
                posSet = ss;
            }
            // if(!par->improve) {
            //     total += coverage;
            //     subsets_coverage.push_back(make_pair(posSet, total));
            // }
        }
        // if(!par->improve && rand() % 8 == 0) {
        //     if(CHECK) cout << "random set" << endl;
        //     rand_subset = ((double) rand()) / RAND_MAX;
        //     for(int i=0; i<subsets_coverage.size(); i++) {
        //         if(rand_subset <= subsets_coverage[i].second / total) {
        //             posSet = subsets_coverage[i].first;
        //             break;
        //         }
        //     }

        //     // posSet = subsets_coverage[rand() % (subsets_coverage.size())];
        // }

        for(int i=0; i<par->nWX; i++) U[i] = U[i] & ~par->bF[posSet][i];

        C.push_back(posSet);

        for(int e : par->F[posSet]) {
            par->mp.erase(remove_if(par->mp.begin(), par->mp.end(), [e](const item& mp) {return mp.value == e;}), par->mp.end());
        }

        if(CHECK) {
            cout << "Best Coverage: " << bestCoverage << endl;
            cout << "Pos. Subset: " << posSet << endl;
            cout << "|U|: " << countSet(U) << endl;
            printSubset(U);
        }
        bestCoverage = 0;
        subsets.clear();
        subsets_coverage.clear();
    }

    return C;
}

void createMap() {
    int pos = 0;
    par->mp = vector<item>(par->n);
    for(pair<int, vector<int>> values : par->inSet){
        // par->elem_pos[values.first] = pos;
        par->mp[pos].value = values.first;
        par->mp[pos].subSets = values.second;
        par->mp[pos].rep = values.second.size();
        pos++;
    }
    sort(par->mp.begin(), par->mp.end(), [&](item a, item b){return a.rep < b.rep;});

    preprocess();

    if(CHECK) {
        for(item mp_item : par->mp) {
            cout << "(" << mp_item.value << ") |" << mp_item.rep << "| => ";
            for (int index : mp_item.subSets) {
                cout << index << " ";
            }
            cout << endl;
        }
    }
}

void preprocess() {
    if(PRINT) {
        cout << "------------------------" << endl;
        cout << "Executing PreSetCover..." << endl;
        cout << "------------------------" << endl;
    }
    // Add uniques elements
    int setIndex;
    ulong* S;
    while(par->mp[0].rep == 1) {
        setIndex = par->mp[0].subSets[0];

        // Eliminar subsets del map que no se usen
        for(int e : par->F[setIndex]) {
            par->mp.erase(remove_if(par->mp.begin(), par->mp.end(), [e](const item& mp) {return mp.value == e;}), par->mp.end());
            cleanBit64(par->X, (e-1));
        }

        par->unique_elements.push_back(setIndex);
    }

    if(PRINT) {
        cout << "Added " << par->unique_elements.size() << " subsets " << endl; 
        cout << "|X| = " << countSet(par->X) << endl;
        cout << "|F| = " << par->bF.size() << endl;
    }
}

int chooseFunction() {
    float total = 0;
    for(int i=0; i<4; i++) total += par->costFunc[i];
    float rand_func = ((float) rand()) / RAND_MAX;
    float cum_score = 0;

    for(int i=0; i<4; i++) {
        cum_score += par->costFunc[i] / total;
        if(rand_func <= cum_score) return i;
    }
    return 0;
}

bool isCovered(vector<int> S) {
    ulong* coveredElements = unionSets(S);
    for (int i = 0; i < par->nWX; i++) if ((coveredElements[i] & par->X[i]) != par->X[i]) {
        delete[] coveredElements;
        return false;
    }
    delete[] coveredElements;
    return true;
}

ulong* unionSets(const vector<int> &S) {
    ulong* C = new ulong[par->nWX];
    fill(C, C + par->nWX, 0);
    for(const int s_idx : S) for(int i=0; i<par->nWX; i++) C[i] |= par->bF[s_idx][i];
    return C;
}

int intersectionLength(const ulong* A, const ulong* B) {
    int cont = 0;
    for(int i=0; i<par->nWX; i++) cont += __builtin_popcountl(A[i] & B[i]);
    return cont;
}

int countSet(const ulong* S){
    int cont = 0;
    for(int i=0; i<par->nWX; i++) {
        cont += __builtin_popcountl(S[i]);
    }
    return cont;
}

void printSubset(const ulong *S) {
    for (int i=0; i<par->nWX; i++){
        printBitsUlong(S[i]);
        cout << " - ";
    }
    cout << endl;
}

void printSubsets(const vector<ulong*> &C) {
    for(ulong* S : C) {
        printSubset(S);
    }
}

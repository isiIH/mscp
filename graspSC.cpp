#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <chrono>
#include "include/BasicCDS.h"
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

#define IT 100
#define RCL 1

typedef struct{
	int value;
	int rep;
	vector<int> subSets;
} item;

typedef struct {
    int posSet;
    int size;
    int elems;
} cvg;

// Structure with all globals parameters program
typedef struct {
	ulong* X;
	vector<vector<int>> F;
    vector<ulong*> bF;

    set<int> chi;
    map<int,int> elem_pos;
    map<int, ulong*> subsets_pos;
    unordered_map<int, vector<int>> inSet;
    vector<item> mp;
    vector<cvg> subset_cvg;
    vector<int> unique_elements;
    vector<int> greedy_sol;
    vector<int> aprox_sol;
    ulong sizeF, sizeNF;
	ulong n, m, nWX;
    int nt;
} ParProg;

ParProg* par;

void readFile(string filename);
void readFileScp(string filename);
void readFilePartition(string filename);
void analyzeF();
void preprocess();

void exhaustive_sol();
void greedy();

double jaccard(const ulong* A, const ulong* B);
int coverageSubset(vector<int> sets, const int pos);
void graspSC();
vector<int> randGreedySC(ulong* U, vector<int> init_sol);
vector<int> randSuccintSC(ulong* U, vector<int> init_sol);

bool isCovered(vector<int> S);
ulong* unionSets(const vector<int> &S);
int countSet(const ulong* S);
int intersectionLength(const ulong* A, const ulong* B);

void printSubset(const ulong *S);
void printSubsets(const vector<ulong*> &C);

int main(int argc, char** argv) {

    if(argc !=2){
		cout << "./opt <filename>" << endl;
		exit(EXIT_FAILURE);
	}

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
    }

    //PREPROCESS
    start_time = chrono::high_resolution_clock::now();
    preprocess();
    end_time = chrono::high_resolution_clock::now();
    auto dur_preprocess = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();

    //GREEDY
    start_time = chrono::high_resolution_clock::now();
    greedy();
    end_time = chrono::high_resolution_clock::now();
    auto dur_greedyExh = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();
    dur_greedyExh += dur_preprocess + dur_analyze;

    //NEW EXHAUSTIVE ALGORITHM
    start_time = chrono::high_resolution_clock::now();
    graspSC();
    // par->aprox_sol = randSuccintSC(par->unique_elements);
    end_time = chrono::high_resolution_clock::now();
    auto dur_apr = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();
    dur_apr += dur_preprocess + dur_analyze;

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

    return 0;
}

void readFile(string filename) {
    if (filename.substr(0,3) == "scp") readFileScp(filename);
    else readFilePartition(filename);
}

void readFileScp(string filename) {
    if(PRINT) cout << "Reading file " << filename << "..." << endl;
    string nametxt = "test_scp/" + filename;
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
    string nametxt = "test_partition/" + filename;
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
    for( int i=0; i<par->F.size(); i++ ) {
        for( int e : par->F[i] ) {
            par->chi.insert(e);
            par->inSet[e].push_back(i);
        }
    }

    par->n = par->chi.size();

    par->nWX = (par->n)/(sizeof(ulong)*8);
    if ((par->n)%(sizeof(ulong)*8)>0) par->nWX++;
    par->X = new ulong[par->nWX];
    fill(par->X, par->X + par->nWX, 0);

    int pos = 0;
    par->mp = vector<item>(par->n);
    for(pair<int, vector<int>> values : par->inSet){
        setBit64(par->X, pos);
        par->elem_pos[values.first] = pos;
        par->mp[pos].value = values.first;
        par->mp[pos].subSets = values.second;
        par->mp[pos].rep = values.second.size();
        pos++;
    }

    ulong *bset;
    for( int i=0; i<par->F.size(); i++ ) {
        bset = new ulong[par->nWX];
        fill(bset, bset + par->nWX, 0);

        for( int e : par->F[i] ) {
            setBit64(bset, par->elem_pos[e]);
        }

        par->subsets_pos[i] = bset;

        par->bF.push_back(bset);
    }

    if(CHECK) {
        for(item mp_item : par->mp) {
            cout << "(" << mp_item.value << ") |" << mp_item.rep << "| => ";
            for (int index : mp_item.subSets) {
                cout << index << " ";
            }
            cout << endl;
        }
    }

    sort(par->mp.begin(), par->mp.end(), [&](item a, item b){return a.rep < b.rep;});

    if(CHECK) {
        cout << "Universe elements = " << endl;
        for( pair<int, int> values : par->elem_pos ) if(getBit64(par->X, values.second)) cout << values.first << " ";
        cout << endl;
        cout << "X = " << countSet(par->X) << endl;
        cout << "n = " << par->n << endl;
        cout << "F = " << par->bF.size() << endl;
        cout << "m = " << par->m << endl;

        for(item mp_item : par->mp) {
            cout << " - " << mp_item.value << " - " << endl;
            cout << mp_item.rep << " subsets." << endl;
            // for (int setIndex : mp_item.subSets) printSubset(par->bF[setIndex]);
        }
    }
    // par->nWX = (par->n)/(sizeof(ulong)*8);
    // if ((par->n)%(sizeof(ulong)*8)>0) par->nWX++;
    
    // par->X = new ulong[par->nWX];
    // fill(par->X, par->X + par->nWX, 0);

    // par->mp = vector<item>(par->n);
    // ulong *bset;
    // for(int i=0; i<par->F.size(); i++){
    //     bset = new ulong[par->nWX];
    //     fill(bset, bset + par->nWX, 0);

    //     for(int e : par->F[i]) {
    //         setBit64(par->X, (e-1));
    //         par->mp[(e-1)].value = e;
    //         par->mp[(e-1)].subSets.push_back(i);

    //         setBit64(bset, (e-1));
    //     }

    //     par->bF.push_back(bset);
    // }

    // for(int i=0; i<par->mp.size(); i++) par->mp[i].rep = par->mp[i].subSets.size();

    // if(CHECK) {
    //     for(item mp_item : par->mp) {
    //         cout << "(" << mp_item.value << ") |" << mp_item.rep << "| => ";
    //         for (int index : mp_item.subSets) {
    //             cout << index << " ";
    //         }
    //         cout << endl;
    //     }
    // }

    // sort(par->mp.begin(), par->mp.end(), [&](item a, item b){return a.rep < b.rep;});

    // if(CHECK) {
    //     cout << "X = " << countSet(par->X) << endl;
    //     cout << "F = " << par->bF.size() << endl;

    //     for(item mp_item : par->mp) {
    //         cout << " - " << mp_item.value << " - " << endl;
    //         cout << mp_item.rep << " subsets." << endl;
    //         // for (int setIndex : mp_item.subSets) printSubset(par->bF[setIndex]);
    //     }
    // }
}

void greedy() {
    int i;
    ulong* U = new ulong[par->nWX];
    for(i=0; i<par->nWX; i++) U[i] = par->X[i];
    vector<int> C = par->unique_elements;
    int maxLengthSS = 0;
    int lengthSS;
    int posSet;

    while( countSet(U) > 0 ) {

        for(pair<int, ulong*> ss_pos : par->subsets_pos){
            lengthSS = intersectionLength(U, ss_pos.second);
            if(lengthSS > maxLengthSS) {
                maxLengthSS = lengthSS;
                posSet = ss_pos.first;
            }
        }

        for(i=0; i<par->nWX; i++) U[i] = U[i] & ~par->subsets_pos[posSet][i];
        C.push_back(posSet);
        par->subsets_pos.erase(posSet);

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

void graspSC() {
    vector<int> new_sol;
    ulong* U = new ulong[par->nWX];
    for(int i=0; i<par->nWX; i++) U[i] = par->X[i];
    int ss;
    ulong* unionSC;
    int cvg;
    int low_cvg = par->n+1;
    int col;
    int nRemove;
    vector<int> setsRemoved;

    //Solución inicial
    // par->aprox_sol = randSuccintSC(U, par->unique_elements);
    par->aprox_sol = randGreedySC(U, par->unique_elements);

    par->mp.clear();

    if(PRINT) cout << "Initial Sol. Cardinality: " << par->aprox_sol.size() << endl;

    for(int iter=0; iter<IT; iter++){
        //Perturbación
        new_sol = par->aprox_sol;
        nRemove = rand() % (int)ceil((new_sol.size() * RCL)) + 1;

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
            // Cálculo de número de elementos cubiertos de cada set si se eliminara del SC
            for(int s_idx = 0; s_idx < new_sol.size(); s_idx++) {
                cvg = coverageSubset(new_sol, s_idx);

                if(CHECK) {
                    cout << "S" << new_sol[s_idx] << ": " << cvg << endl;
                }

                if(cvg < low_cvg) {
                    low_cvg = cvg;
                    col = s_idx;
                }
            }

            if(CHECK) {
                cout << "Removing S" << new_sol[col] << endl;
            }
            setsRemoved.push_back(new_sol[col]);
            new_sol.erase(new_sol.begin() + col);


            low_cvg = par->n+1;
        }

        // Actualizar U y map
        unionSC = unionSets(new_sol);
        for(int ss : setsRemoved)  {
            for(int e : par->F[ss]) {
                if(checkBit(U, par->elem_pos[e]) == 0 && checkBit(unionSC, par->elem_pos[e]) == 0) {
                    item it_map;
                    it_map.value = e;
                    it_map.subSets = par->inSet[e];
                    it_map.rep = it_map.subSets.size();
                    par->mp.push_back(it_map);
                    setBit64(U, par->elem_pos[e]);
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
        // new_sol = randSuccintSC(U, new_sol);
        new_sol = randSuccintSC(U, new_sol);

        if(new_sol.size() < par->aprox_sol.size()) {
            par->aprox_sol = new_sol;
        }

        if(PRINT) {
            cout << "Sol. Cardinality: " << new_sol.size() << endl;
            // printSubsets(new_sol);
            cout << "Best Cardinality: " << par->aprox_sol.size() << endl;
        }
    }

    // for(ulong* ss: par->unique_elements) par->aprox_sol.push_back(ss);
    
}

vector<int> randGreedySC(ulong* U, vector<int> init_sol) {
    int i;
    vector<int> C = init_sol;
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

    return C;
}

vector<int> randSuccintSC(ulong* U, vector<int> init_sol) {
    vector<int> C = init_sol;
    int posSet;
    vector<int> subsets;
    double coverage;
    double bestCoverage = 0;
    int grade;
    int p;

    while( countSet(U) > 0 ) {
        grade = par->mp[0].rep;
        p = 0;
        while(p < par->mp.size() && par->mp[p].rep == grade) {
            for(int ss : par->mp[p].subSets) subsets.push_back(ss);
            p++;
        }
        if(1) {
            for(int ss : subsets) {
                // cout << ss << endl;
                // cout << jaccard(unionSC, subsets[i]) << endl;
                // printSubset(subsets[i]);

                // coverage = jaccard(unionSC, par->bF[subsets[i]]);
                coverage = intersectionLength(U, par->bF[ss]);

                if(coverage > bestCoverage) {
                    bestCoverage = coverage;
                    posSet = ss;
                }
            }
        } else {
            posSet = subsets[rand() % (subsets.size())];
        }

        for(int i=0; i<par->nWX; i++) U[i] = U[i] & ~par->bF[posSet][i];

        C.push_back(posSet);

        for(int e : par->F[posSet]) {
            par->mp.erase(remove_if(par->mp.begin(), par->mp.end(), [e](const item& mp) {return mp.value == e;}), par->mp.end());
        }

        if(CHECK) {
            cout << "Best Coverage: " << bestCoverage << endl;
            cout << "Pos. Subset: " << posSet << endl;
            cout << "|U|: " << countSet(U) << endl;
        }
        bestCoverage = 0;
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
    // Add uniques elements
    int setIndex;
    ulong* S;
    while(par->mp[0].rep == 1) {
        setIndex = par->mp[0].subSets[0];

        // Eliminar subsets del map que no se usen
        for(int e : par->F[setIndex]) {
            par->mp.erase(remove_if(par->mp.begin(), par->mp.end(), [e](const item& mp) {return mp.value == e;}), par->mp.end());
            cleanBit64(par->X,par->elem_pos[e]);
        }

        par->unique_elements.push_back(setIndex);
        // par->bF.erase(find(par->bF.begin(), par->bF.end(), S));
        // par->m--;
    }

    // Inicializar vector de subconjuntos con numero de elementos con cierto grado
    // par->subset_cvg = vector<cvg>(par->m);
    // for(int i=0; i<par->m; i++) {
    //     par->subset_cvg[i].posSet = i;
    //     par->subset_cvg[i].size = countSet(par->bF[i]);
    //     par->subset_cvg[i].elems = 0;
    // }
    // for( item mp_item : par->mp ) {
    //     if(mp_item.rep >= GRADE) {
    //         for(int ss : mp_item.subSets) {
    //             par->subset_cvg[ss].elems++;
    //         }
    //     }
    // } 

    // sort(par->subset_cvg.begin(), par->subset_cvg.end(), [&](cvg a, cvg b){return (a.elems > b.elems) || (a.elems == b.elems && a.size > b.size);});
    // while(par->subset_cvg[par->subset_cvg.size()-1].elems == 0) par->subset_cvg.pop_back();

    if(PRINT) {
        cout << "Added " << par->unique_elements.size() << " subsets " << endl; 
        cout << "|X| = " << countSet(par->X) << endl;
        cout << "|F| = " << par->bF.size() << endl;
        for(item mp_item : par->mp) {
            cout << " - " << mp_item.value << " - " << endl;
            cout << mp_item.rep << " subsets." << endl;
            // for (int setIndex : mp_item.subSets) printSubset(par->bF[setIndex]);
        }
        // for(cvg mp_ss : par->subset_cvg) {
        //     cout << "S" << (mp_ss.posSet+1) << ": |" << mp_ss.elems << "| elems with grade >= " << GRADE << endl;
        // }
    }
}

int coverageSubset(vector<int> sets, int pos) {
    ulong* S = par->bF[sets[pos]];
    sets.erase(sets.begin() + pos);
    ulong* sets_union = unionSets(sets);
    int cont = 0;
    for(int i=0; i<par->nWX; i++) cont += __builtin_popcountl(S[i] & ~sets_union[i]);
    return cont;
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

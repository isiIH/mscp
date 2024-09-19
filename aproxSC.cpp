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

using namespace std;
using namespace cds;

#define PRINT 1
#define CHECK 1

#define GRADE 2

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

    vector<item> mp;
    vector<cvg> subset_cvg;
    vector<ulong*> unique_elements;
    vector<ulong*> greedy_sol;
    vector<ulong*> aprox_sol;
    ulong sizeF, sizeNF;
	ulong n, m, nWX;
    int nt;
} ParProg;

ParProg* par;

void readFile(string filename);
void readFileScp(string filename);
void readFilePartition(string filename);
void analizeF();
void preprocess();

void exhaustive_sol();
void greedy();

double jaccard(const ulong* A, const ulong* B);

bool isCovered(vector<ulong*> chosenSets);
ulong* unionF(const vector<ulong*> &F);
int countSet(const ulong* S);
int intersectionLength(const ulong* A, const ulong* B);

void printSubset(const ulong *S);
void printSubsets(const vector<ulong*> &C);

int main(int argc, char** argv) {

    if(argc !=2){
		cout << "./opt <filename>" << endl;
		exit(EXIT_FAILURE);
	}

    par = new ParProg();

    readFile(argv[1]);
    auto start_time = chrono::high_resolution_clock::now();
    analizeF();
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
    
    end_time = chrono::high_resolution_clock::now();
    auto dur_apr = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();
    dur_apr += dur_preprocess + dur_analyze;

    if(CHECK) {
        cout << "SOL: " << endl;
        printSubsets(par->aprox_sol);
    }
    if(PRINT) {
        cout << "------------------------" << endl;
        cout << "Greedy Cardinality: " << par->greedy_sol.size() << endl;
        cout << "Time [s]: " << dur_greedyExh/1000000.0 << endl;
        cout << "AproxSC Cardinality: " << par->aprox_sol.size() << endl;
        cout << "Time [s]: " << dur_apr/1000000.0 << endl;

        // printSubsets(par->greedy_sol);
    }

    cout << argv[1] << " " << par->n << " " << par->m << " " << dur_apr/1000000.0 << " " << par->aprox_sol.size() << " " << par->greedy_sol.size() << endl;

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

void analizeF() {
    par->nWX = (par->n)/(sizeof(ulong)*8);
    if ((par->n)%(sizeof(ulong)*8)>0) par->nWX++;
    
    par->X = new ulong[par->nWX];
    fill(par->X, par->X + par->nWX, 0);

    par->mp = vector<item>(par->n);
    ulong *bset;
    for(int i=0; i<par->F.size(); i++){
        bset = new ulong[par->nWX];
        fill(bset, bset + par->nWX, 0);

        for(int e : par->F[i]) {
            setBit64(par->X, (e-1));
            par->mp[(e-1)].value = e;
            par->mp[(e-1)].subSets.push_back(i);

            setBit64(bset, (e-1));
        }

        par->bF.push_back(bset);
    }

    for(int i=0; i<par->mp.size(); i++) par->mp[i].rep = par->mp[i].subSets.size();

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
}

void greedy() {
    int i;
    ulong* U = new ulong[par->nWX];
    for(i=0; i<par->nWX; i++) U[i] = par->X[i];
    vector<ulong*> subsets = par->bF;
    vector<ulong*> C = par->unique_elements;
    int posSet;
    int maxLengthSS = 0;
    int lengthSS;

    while( countSet(U) > 0 ) {

        for( i=0; i<subsets.size(); i++ ){
            lengthSS = intersectionLength(U, subsets[i]);
            if(lengthSS > maxLengthSS) {
                maxLengthSS = lengthSS;
                posSet = i;
            }
        }

        for(i=0; i<par->nWX; i++) U[i] = U[i] & ~subsets[posSet][i];
        C.push_back(subsets[posSet]);
        subsets.erase(subsets.begin()+posSet);

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
        S = par->bF[setIndex];

        // Eliminar subsets del map que no se usen
        for(int e : par->F[setIndex]) {
            par->mp.erase(remove_if(par->mp.begin(), par->mp.end(), [e](const item& mp) {return mp.value == e;}), par->mp.end());
            cleanBit64(par->X,(e-1));
        }

        par->unique_elements.push_back(S);
        // par->bF.erase(find(par->bF.begin(), par->bF.end(), S));
        // par->m--;
    }

    // Inicializar vector de subconjuntos con numero de elementos con cierto grado
    par->subset_cvg = vector<cvg>(par->m);
    for(int i=0; i<par->m; i++) {
        par->subset_cvg[i].posSet = i;
        par->subset_cvg[i].size = countSet(par->bF[i]);
        par->subset_cvg[i].elems = 0;
    }
    for( item mp_item : par->mp ) {
        if(mp_item.rep >= GRADE) {
            for(int ss : mp_item.subSets) {
                par->subset_cvg[ss].elems++;
            }
        }
    } 

    sort(par->subset_cvg.begin(), par->subset_cvg.end(), [&](cvg a, cvg b){return (a.elems > b.elems) || (a.elems == b.elems && a.size > b.size);});
    while(par->subset_cvg[par->subset_cvg.size()-1].elems == 0) par->subset_cvg.pop_back();

    par->n = countSet(par->X);
    if(PRINT) {
        cout << "Added " << par->unique_elements.size() << " subsets " << endl; 
        cout << "|X| = " << par->n << endl;
        cout << "|F| = " << par->m << endl;
        for(item mp_item : par->mp) {
            cout << " - " << mp_item.value << " - " << endl;
            cout << mp_item.rep << " subsets." << endl;
            // for (int setIndex : mp_item.subSets) printSubset(par->bF[setIndex]);
        }
        for(cvg mp_ss : par->subset_cvg) {
            cout << "S" << (mp_ss.posSet+1) << ": |" << mp_ss.elems << "| elems with grade >= " << GRADE << endl;
        }
    }
}

bool isCovered(vector<ulong*> chosenSets) {
    ulong* coveredElements = unionF(chosenSets);
    for (int i = 0; i < par->nWX; i++) if ((coveredElements[i] & par->X[i]) != par->X[i]) {
        delete[] coveredElements;
        return false;
    }
    delete[] coveredElements;
    return true;
}

ulong* unionF(const vector<ulong*> &F) {
    ulong* C = new ulong[par->nWX];
    fill(C, C + par->nWX, 0);
    for(const ulong* subset : F) for(int i=0; i<par->nWX; i++) C[i] |= subset[i];
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
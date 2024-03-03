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

using namespace std;
using namespace cds;

#define PRINT 1
#define CHECK 0

typedef struct{
	int value;
	int rep;
	vector<int> subSets;
} item;

// Structure with all globals parameters program
typedef struct {
    int search;
	ulong* X;
	vector<vector<int>> F;
    set<int> chi;
    map<int,int> elem_pos;
    vector<ulong*> bF;
    vector<item> mp;
    vector<ulong*> unique_elements;
    vector<ulong*> greedy_sol;
    vector<ulong*> exh_sol;
    ulong sizeF, sizeNF;
	ulong n, m, nWX;
} ParProg;

ParProg* par;

void readFile(string filename);
void readFileScp(string filename);
void readFilePartition(string filename);
void analizeF();
void preprocess();

bool kSol(int i, int k, vector<ulong*> chosenSets);
void binarySearch(int mi, int ma, vector<ulong*> &chosenSets);

void exhaustive_sol();
void greedy();

bool isCovered(vector<ulong*> chosenSets);
ulong* unionF(const vector<ulong*> &F);
int countSet(const ulong* S);
int intersectionLength(const ulong* A, const ulong* B);

void printSubset(const ulong *S);
void printSubsets(const vector<ulong*> &C);

int main(int argc, char** argv) {

    if(argc !=3){
		cout << "./opt <filename> <search>" << endl;
		exit(EXIT_FAILURE);
	}

    par = new ParProg();
    par->search = atoi(argv[2]);
    if(par->search < 0 || par->search > 2){
        cout << "Invalid Search Type!\n0: Secuential Search\n1: Binary Search\n2: Exponential Search" << endl;
        exit(EXIT_FAILURE);
    }

    readFile(argv[1]);
    auto start_time = chrono::high_resolution_clock::now();
    analizeF();
    auto end_time = chrono::high_resolution_clock::now();
    auto dur_analyze = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();


    cout  << "X: " << par->n
        << " | F: " << par->m << endl;

    par->sizeF = par->m*sizeof(ulong)*par->n;
    par->sizeNF = par->m*sizeof(ulong)*par->nWX;

	cout << "nWX = " << par->nWX << endl;
	cout << " size for F[] = " << par->sizeF/(1024.0*1024.0) << " MiB" << endl;
	cout << " size for nF[] = " << par->sizeNF/(1024.0*1024.0) << " MiB" << endl;

    if(CHECK) {
        for(vector<int> set : par->F) {
            for(int val : set) {
                cout << val << " - ";
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
    exhaustive_sol();
    end_time = chrono::high_resolution_clock::now();
    auto dur_opt = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();
    dur_opt += dur_preprocess + dur_analyze;

    if(CHECK) {
        cout << "SOL: " << endl;
        printSubsets(par->exh_sol);
    }
    cout << "Greedy Cardinality: " << par->greedy_sol.size() << endl;
    cout << "Time [μs]: " << dur_greedyExh << endl;
	cout << "Optimal Cardinality: " << par->exh_sol.size() << endl;
    cout << "Time [μs]: " << dur_opt << endl;

    return 0;
}

void readFile(string filename) {
    if (filename.substr(0,3) == "scp") readFileScp(filename);
    else readFilePartition(filename);
}

void readFileScp(string filename) {
    cout << "Reading file " << filename << "..." << endl;
    string nametxt = "scp/" + filename;
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
    unordered_map<int, vector<int>> inSet;
    for( int i=0; i<par->F.size(); i++ ) {
        for( int e : par->F[i] ) {
            par->chi.insert(e);
            inSet[e].push_back(i);
        }
    }

    par->n = par->chi.size();

    par->nWX = (par->n)/(sizeof(ulong)*8);
    if ((par->n)%(sizeof(ulong)*8)>0) par->nWX++;
    par->X = new ulong[par->nWX];
    fill(par->X, par->X + par->nWX, 0);

    int pos = 0;
    par->mp = vector<item>(par->n);
    for(pair<int, vector<int>> values : inSet){
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
            inSet[e].push_back(i);
        }

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

void exhaustive_sol() {
    cout << "-----------------------------------------------------------" << endl;
    cout << "Executing new exhaustive algorithm";
    switch (par->search) {
        case 0: cout << " with sequential search..." << endl;
        break;
        case 1: cout << " with binary search..." << endl;
        break;
        case 2: cout << " with exponential search..." << endl;
    }
    cout << "-----------------------------------------------------------" << endl;

    // Calcular mínimo número de subconjuntos
    vector<ulong*> Fsort = par->bF;
    sort(Fsort.begin(), Fsort.end(), [&](ulong* a, ulong* b){return countSet(a) > countSet(b);});
    int minSS = 0;
    int k = 0;
    while(minSS < par->n) {
        minSS += countSet(Fsort[k]);
        k++;
    }

    //Iterar desde K hasta encontrar el óptimo
    vector<ulong*> chosenSets;
    cout << "Search Range = [" << k << " - " << par->greedy_sol.size() - par->unique_elements.size() << "]" << endl;

    //Búsqueda secuencial
    if(par->search == 0){
        bool found = false;
        while(!found) {
            if(PRINT) cout << "K = " << k << endl;
            auto start = chrono::high_resolution_clock::now();
            found = kSol(0, k, chosenSets);
            auto end = chrono::high_resolution_clock::now();
            cout << "Time [μs]: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
            //No se encuentra una solución de tamaño k, aumentamos en 1
            k++;
        }
    //Búsqueda binaria
    } else if(par->search == 1) {
        int ma = par->greedy_sol.size() - par->unique_elements.size();
        binarySearch(k, ma, chosenSets);
    //Búsqueda exponencial
    } else if(par->search == 2) {
        int exp = 1;
        int greedySize = par->greedy_sol.size() - par->unique_elements.size();
        bool found = false;
        while(k <= greedySize && !found) {
            cout << "K = " << k << endl;
            auto start = chrono::high_resolution_clock::now();
            found = kSol(0, k, chosenSets);
            auto end = chrono::high_resolution_clock::now();
            cout << "Time [μs]: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
            
            if(!found){
                k += exp;
                exp *= 2;
            }
        }

        //Realizar búsqueda binaria en un rango más pequeño
        int mi = k - exp/2 + 1;
        int ma = min(k-1, greedySize);
        cout << "Search range for binary search: [" << mi << " - " << ma << "]" << endl;
        binarySearch(mi, ma, chosenSets);
    }

    par->exh_sol.insert(par->exh_sol.end(), par->unique_elements.begin(), par->unique_elements.end());
}

bool kSol(int index, int k, vector<ulong*> chosenSets) {
    //|chosenSets| = k 
    if(k == 0) {
        //Verificar si se cubre el universo
        if(isCovered(chosenSets)) {
            cout << "Solution found with K = " << chosenSets.size() << endl;
            // for(ulong* S : chosenSets) par->exh_sol.push_back(S);
            par->exh_sol = chosenSets;
            return true;
        }
        return false;
    }

    for (int j = index; j<par->m-k+1; j++) {
        //Incluir el subconjunto
        chosenSets.push_back(par->bF[j]);

        if ( kSol(j+1, k-1, chosenSets) ) return true;
        
        //No incluir el subconjunto
        chosenSets.pop_back();
    }

    return false;
}

void binarySearch(int mi, int ma, vector<ulong*> &chosenSets) {
    int m = mi;
    bool found;
    while(mi <= ma) {
        if(PRINT) cout << "K = " << m << endl;
        auto start = chrono::high_resolution_clock::now();
        found = kSol(0, m, chosenSets);
        auto end = chrono::high_resolution_clock::now();
        cout << "Time [μs]: " << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
        if (found) ma = m - 1;
        else mi = m + 1;
        m = mi + (ma - mi)/2;
    }
}

void preprocess() {
    cout << "------------------------" << endl;
    cout << "Executing PreSetCover..." << endl;
    cout << "------------------------" << endl;
    // Add uniques elements
    int setIndex;
    ulong* S;
    while(par->mp[0].rep == 1) {
        setIndex = par->mp[0].subSets[0];
        S = par->bF[setIndex];

        // Eliminar subsets del map que no se usen
        for(int e : par->F[setIndex]) {
            par->mp.erase(remove_if(par->mp.begin(), par->mp.end(), [e](const item& mp) {return mp.value == e;}), par->mp.end());
            cleanBit64(par->X,par->elem_pos[e]);
        }
        par->unique_elements.push_back(S);
    }

    for(ulong* S : par->unique_elements) {
        par->bF.erase(find(par->bF.begin(), par->bF.end(), S));
        par->m--;
    }

    cout << "Added " << par->unique_elements.size() << " subsets " << endl; 
    par->n = countSet(par->X);
    cout << "|X| = " << par->n << endl;
    cout << "|F| = " << par->m << endl;
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
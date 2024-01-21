#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <chrono>
#include "include/BasicCDS.h"
#include <cmath>

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
    vector<ulong*> bF;
    vector<item> mp;
    vector<ulong*> unique_elements;
    vector<ulong*> greedy_sol;
    vector<ulong*> exh_sol;
    ulong sizeF, sizeNF;
	ulong n, m, nWX, nWF;
} ParProg;

ParProg* par;

void readFile(string filename);
void preprocess();
void createMap();
bool kSol(int i, int k, vector<ulong*> chosenSets);
void binarySearch(int m, int mi, int ma, vector<ulong*> &chosenSets);

void exhaustive_sol();
void greedy();

ulong* unionF(vector<ulong*> &F);
int countSet(ulong* S);
int intersectionLength(ulong* A, ulong* B);

void printSubset(ulong *S, int size = par->nWX);
void printSubsets(vector<ulong*> C);

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

    cout  << "X: " << par->n
        << " | F: " << par->m << endl;

    par->sizeF = par->m*sizeof(ulong)*par->n;
    par->sizeNF = par->m*sizeof(ulong)*par->nWX;

	cout << "nWX = " << par->nWX << endl;
	cout << " size for F[] = " << par->sizeF/(1024.0*1024.0) << " MiB (using ulong per cell)" << endl;
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

    auto start_time = chrono::high_resolution_clock::now();
    preprocess();
    auto end_time = chrono::high_resolution_clock::now();
	cout << "Duración en milisegundos: " << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000.0 << endl;

    //SOL. CLASSIC GREEDY ALGORITHM
    start_time = chrono::high_resolution_clock::now();
    greedy();
    end_time = chrono::high_resolution_clock::now();

    if(CHECK) {
        cout << "SOL: " << endl;
        printSubsets(par->greedy_sol);
    }
    cout << "Solution Cardinality: " << par->greedy_sol.size() << endl;

    cout << "Duración en milisegundos: " << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000.0 << endl;


    //NEW EXHAUSTIVE ALGORITHM
    start_time = chrono::high_resolution_clock::now();
    exhaustive_sol();
    end_time = chrono::high_resolution_clock::now();

    if(CHECK) {
        cout << "SOL: " << endl;
        printSubsets(par->exh_sol);
    }
    cout << "Solution Cardinality: " << par->exh_sol.size() << endl;

    cout << "Duración en segundos: " << chrono::duration_cast<chrono::seconds>(end_time - start_time).count() << endl;
	cout << "Duración en microsegundos: " << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count() << endl;

    return 0;
}

void readFile(string filename) {
    cout << "Reading file " << filename << "..." << endl;
    string nametxt = "test/" + filename;
    ifstream file(nametxt.c_str());
    if(file.fail()){
        cout << "File not found!" << endl;
        exit(EXIT_FAILURE);
    }
    string line,item;

    //m & n
	getline(file,line);
    istringstream ss(line);
    ss >> (par->m) >> (par->n);

    par->nWX = (par->n)/(sizeof(ulong)*8);
    if ((par->n)%(sizeof(ulong)*8)>0) par->nWX++;
    par->X = new ulong[par->nWX];

    //Costs
    int i = 0;
    while(i < par->n) {
        getline(file>>std::ws,line);
        istringstream iss(line);
        while (getline(iss, item, ' ')){i++;}
    }

    //Sets
    int numCover;
    i = 0;
    int j;
    ulong *bset;
    vector<int> set;
    while(i < par->m) {
        getline(file>>std::ws,line);
        numCover = stoi(line);

        bset = new ulong[par->nWX];
        for (j=0; j < par->nWX; j++) bset[j] = 0;

        j = 0;
        while(j < numCover){
            getline(file>>std::ws,line);
            istringstream iss(line);
            while (getline(iss, item, ' ')) {
                set.push_back(stoi(item));
                setBit64(bset, stoi(item)-1);
                setBit64(par->X, stoi(item)-1);
                j++;
            }
        }
        (par->F).push_back(set);
        set.clear();
        (par->bF).push_back(bset);
        i++;
    }
}

void greedy() {
    cout << "-------------------------------------" << endl;
    cout << "Executing classic greedy algorithm..." << endl;
    cout << "-------------------------------------" << endl;
    ulong* U = new ulong[par->nWX];
    int i;
    for(i=0; i<par->nWX; i++) U[i] = par->X[i];
    vector<ulong*> C = par->unique_elements;
    ulong* maxS;
    int maxLengthSS = 0;
    int lengthSS;

    while( countSet(U) > 0 ) {

        for( ulong* S : par->bF ){
            if( find(C.begin(), C.end(), S) == C.end()) {
                lengthSS = intersectionLength(U, S);
                if(lengthSS > maxLengthSS) {
                    maxLengthSS = lengthSS;
                    maxS = S;
                }
            }
        }

        for(i=0; i<par->nWX; i++) U[i] = U[i] & ~maxS[i];
        C.push_back(maxS);

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

    //Búsqueda secuencial
    if(par->search == 0){
        while(true) {
            if(PRINT) cout << "K = " << k << endl;
            if ( kSol(0, k, chosenSets) ) break;
            //No se encuentra una solución de tamaño k, aumentamos en 1
            k++;
        }
    //Búsqueda binaria
    } else if(par->search == 1) {
        int ma = par->greedy_sol.size() - par->unique_elements.size();
        binarySearch(k, k, ma, chosenSets);
    //Búsqueda exponencial
    } else if(par->search == 2) {
        int exp = 1;
        int greedySize = par->greedy_sol.size() - par->unique_elements.size();
        while(k <= greedySize && !kSol(0, k, chosenSets)) {
            cout << "K = " << k << endl;
            k += exp;
            exp *= 2;
        }

        //Realizar búsqueda binaria en un rango más pequeño
        int mi = k - exp/2 + 1;
        int ma = min(k-1, greedySize);
        binarySearch(mi, mi, ma, chosenSets);
    }

    par->exh_sol.insert(par->exh_sol.end(), par->unique_elements.begin(), par->unique_elements.end());
}

bool kSol(int index, int k, vector<ulong*> chosenSets) {
    //|chosenSets| = k 
    if(k == 0) {
        //Verificar si se cubre el universo
        ulong* coveredElements = unionF(chosenSets);
        ulong* xCover = new ulong[par->nWX]; 
        for(int i=0; i<par->nWX; i++) xCover[i] = coveredElements[i] & par->X[i];
        if(countSet(xCover) == countSet(par->X)) {
            cout << "Solution found with K = " << chosenSets.size() << endl;
            // for(ulong* S : chosenSets) par->exh_sol.push_back(S);
            par->exh_sol = chosenSets;
            delete[] coveredElements;
            delete[] xCover;
            return true;
        } 
        delete[] coveredElements;
        delete[] xCover;
        return false;
    }

    for (int j = index; j<par->bF.size()-k+1; j++) {
        //Incluir el subconjunto
        chosenSets.push_back(par->bF[j]);

        if ( kSol(j+1, k-1, chosenSets) ) return true;
        
        //No incluir el subconjunto
        chosenSets.pop_back();
    }

    return false;
}

void preprocess() {
    cout << "------------------------" << endl;
    cout << "Executing PreSetCover..." << endl;
    cout << "------------------------" << endl;
    // Create a map structure for each element of the universe
    createMap();
    if(CHECK) {
        for(item mp_item : par->mp) {
            cout << " - " << mp_item.value << " - " << endl;
            cout << mp_item.rep << " subsets." << endl;
            for (int setIndex : mp_item.subSets) printSubset(par->bF[setIndex]);
        }
    }

    // Add uniques elements
    int setIndex;
    ulong* S;
    while(par->mp[0].rep == 1) {
        setIndex = par->mp[0].subSets[0];
        S = par->bF[setIndex];

        // Eliminar subsets del map que no se usen
        for(int e : par->F[setIndex]) {
            for (int j=0; j<par->mp.size(); j++) {
                if (par->mp[j].value == e) {
                    par->mp.erase(par->mp.begin()+j);
                }
            }
            cleanBit64(par->X,e-1);
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

void binarySearch(int m, int mi, int ma, vector<ulong*> &chosenSets) {
    while(mi <= ma) {
        if(PRINT) cout << "K = " << m << endl;
        if (kSol(0, m, chosenSets)) ma = m - 1;
        else mi = m + 1;
        m = mi + (ma - mi)/2;
    }
}

void createMap() {
    par->mp = vector<item>(par->n);
    for(int i=0; i<par->n; i++) {
        par->mp[i].value = i+1;
        for( int j = 0; j<par->bF.size(); j++ ) if(getBit64(par->bF[j], i)) par->mp[i].subSets.push_back(j);
        par->mp[i].rep = par->mp[i].subSets.size();
    }
    sort(par->mp.begin(), par->mp.end(), [&](item a, item b){return a.rep < b.rep;});
}

ulong* unionF(vector<ulong*> &F) {
    int i,j;
    ulong* C = new ulong[par->nWX];
    for(j=0; j<par->nWX; j++) C[j] = 0;

    for(i=0; i<F.size(); i++) {
        for(j=0; j<par->nWX; j++) C[j] = C[j] | F[i][j];
    }

    return C;
}

int intersectionLength(ulong* A, ulong* B) {
    ulong* interSS = new ulong[par->nWX];
    for(int i=0; i<par->nWX; i++) interSS[i] = A[i] & B[i];
    return countSet(interSS);
}

int countSet(ulong* S){
    int cont = 0;
    for(int i=0; i<par->nWX; i++) {
        cont += __builtin_popcountl(S[i]);
    }

    return cont;
}

void printSubset(ulong *S, int size) {
    for (int i=0; i<size; i++){
        printBitsUlong(S[i]);
        cout << " - ";
    }
    cout << endl;
}

void printSubsets(vector<ulong*> C) {
    for(ulong* S : C) {
        printSubset(S);
    }
}
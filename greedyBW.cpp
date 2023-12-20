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

#define PRINT 0
#define CHECK 0

typedef struct{
	int value;
	int rep;
	vector<ulong*> subSets;
} item;

// Structure with all globals parameters program
typedef struct {
	ulong* X;
	vector<vector<int>> F;
    vector<ulong*> bF;
    vector<item> mp;
    vector<ulong*> greedy_sol;
    vector<ulong*> greedy2_sol;
    ulong sizeF, sizeNF;
	ulong n, m, nWX, nWF;
} ParProg;

ParProg* par;

void readFile(string filename);

void greedy();

void greedy2();
vector<ulong*> exhaustiveSC(ulong* X, ulong *U, vector<ulong*> &F);
void createMap();
vector<ulong*> sets(vector<ulong*> &F, const int e);
vector<ulong*> setsOfLength(vector<item> &mp, int &i, const int l);

int intersectionLength(ulong* A, ulong* B);
ulong* unionF(vector<ulong*> &F);
int countSet(ulong* S);

void printSubset(ulong* C);
void printSubsets(vector<ulong*> C);

int main(int argc, char** argv) {

    if(argc !=2){
		cout << "./greedy <filename>" << endl;
		exit(EXIT_FAILURE);
	}

    par = new ParProg();

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

    createMap();
    if(CHECK) {
        for(item mp_item : par->mp) {
            cout << " - " << mp_item.value << " - " << endl;
            cout << mp_item.rep << " subsets." << endl;
            printSubsets(mp_item.subSets);
        }
    }

    //SOL. CLASSIC GREEDY ALGORITHM
    cout << "-------------------------" << endl;
    cout << "Executing classic greedy algorithm..." << endl;

    auto start_time = chrono::high_resolution_clock::now();
    greedy();
    auto end_time = chrono::high_resolution_clock::now();

    if(PRINT) {
        cout << "SOL: " << endl;
        printSubsets(par->greedy_sol);
    }
    cout << "Solution Cardinality: " << par->greedy_sol.size() << endl;

    cout << "Duración en segundos: " << chrono::duration_cast<chrono::seconds>(end_time - start_time).count() << endl;
	cout << "Duración en microsegundos: " << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count() << endl;

    //SOL. NEW GREEDY ALGORITHM
    cout << "-------------------------" << endl;

    start_time = chrono::high_resolution_clock::now();
    greedy2();
    end_time = chrono::high_resolution_clock::now();

    if(PRINT) {
        cout << "SOL: " << endl;
        printSubsets(par->greedy2_sol);
    }
    cout << "Solution Cardinality: " << par->greedy2_sol.size() << endl;

    cout << "Duración en segundos: " << chrono::duration_cast<chrono::seconds>(end_time - start_time).count() << endl;
	cout << "Duración en microsegundos: " << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count() << endl;

    return 0;
}

void readFile(string filename) {
    cout << "Reading file " << filename << "..." << endl;
    string nametxt = "test/" + filename;
    ifstream file(nametxt.c_str());
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
    while(i < par->n)
    {
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
    ulong* U = new ulong[par->nWX];
    int i;
    for(i=0; i<par->nWX; i++) U[i] = par->X[i];
    vector<ulong*> C;
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

void greedy2(){
    cout << "Executing new greedy algorithm..." << endl;
    ulong* U = new ulong[par->nWX];
    int i;
    for(i=0; i<par->nWX; i++) U[i] = par->X[i];
    vector<ulong*> C;
    int j;
    int k = 1;
    vector<ulong*> subF;
    ulong* subU;
    vector<ulong*> solExhaustive;

    while(countSet(U) > 0) {
        subF = setsOfLength(par->mp,i,k);

        if(!subF.empty()){
            subU = unionF(subF);
            if(countSet(subU) > countSet(U)) {
                for(j=0; j<par->nWX; j++) subU[j] = subU[j] & U[j];
            }
            if(CHECK) {
                cout << "----------------------" << endl;
                cout << "Rep: " << k << endl;
                cout << "U: " << countSet(U) << endl;
                cout << "SubU: " << countSet(subU) << endl;
                cout << "SubF: " << subF.size() << endl;
            }

            solExhaustive = exhaustiveSC(U, subU, subF);
            if(CHECK) {
                cout << "Sol. Exhaustiva: " << solExhaustive.size() << endl;
            }

            for (ulong* S : solExhaustive){
                for(j=0; j<par->nWX; j++) U[j] = U[j] & ~S[j];
                for(item &it : par->mp){
                    if(find(it.subSets.begin(), it.subSets.end(), S) != it.subSets.end()) {
                        it.subSets.clear();
                    }
                }
                C.push_back(S);
            }
        }
        k++;
    }

    par->greedy2_sol = C;

}

vector<ulong*> eshaustiveSC(ulong* X, vector<ulong*> &F) {
    
}

// vector<ulong*> exhaustiveSC(ulong* X, ulong *U, vector<ulong*> &F) {
//     vector<ulong*> C;
//     int n = countSet(U);
//     int m = F.size();
//     int bestLen = m+1;
//     int MaxENotCover = 0;
//     vector<ulong*> bestC;
//     int k;

//     ulong* unionC;
//     ulong* xCover = new ulong[par->nWX];
//     ulong* eNotCover = new ulong[par->nWX];

//     for(ulong i=1;i<(1<<m);i++){
//         for(ulong j=0;j<m;j++){
//             if(i&(1<<j)){
//                 C.push_back(F[j]);
//             }
//         }

//         unionC = unionF(C);
//         for(k=0; k<par->nWX; k++) {
//             xCover[k] = unionC[k] | U[k];
//             eNotCover[k] = unionC[k] | X[k];
//         }

//         if(countSet(xCover) == n && C.size() <= bestLen && countSet(eNotCover) > MaxENotCover){
//             bestLen = C.size();
//             bestC = C;
//             MaxENotCover = countSet(eNotCover);
//         }

//         C.clear();
//     }

//     return bestC;
// }

void createMap() {
    par->mp = vector<item>(par->n);
    for(int i=0; i<par->n; i++) {
        par->mp[i].value = i+1;
        par->mp[i].subSets = sets(par->bF,i);
        par->mp[i].rep = par->mp[i].subSets.size();
    }
    sort(par->mp.begin(), par->mp.end(), [&](item a, item b){return a.rep < b.rep;});
}

vector<ulong*> sets(vector<ulong*> &F, const int e) {
    vector<ulong*> C;
    for( ulong* S : F ) if(getBit64(S, e)) C.push_back(S);

    return C;
}

vector<ulong*> setsOfLength(vector<item> &mp, int &i, const int l) {
    vector<ulong*> C;
    while (mp[i].rep == l) {
        for (ulong* S : mp[i].subSets){
            if(find(C.begin(), C.end(), S) == C.end()) {
                C.push_back(S);
            }
        }
        i++;
    }

    return C;
}

ulong* unionF(vector<ulong*> &F) {
    int i,j;
    ulong* C = new ulong[par->nWX];

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
        uint cnt = 0;
        ulong mask = 0x8000000000000000;

        for(cnt=0;cnt<W64;++cnt){
            if((S[i] & mask) != 0) cont++;
            mask >>= 1;
        }
    }

    return cont;
}

void printSubset(ulong *S) {
    for (int i=0; i<par->nWX; i++){
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

// /*
// PROYECTO #1: Adaptar solución greedy exhaustiva a bitwise operator

// probar con n = 64 primero, luego agregar más celdas

// PROYECTO #2: Óptimo que parte de la mímina cantidad posible de subconjuntos
// min <= Opt <= Greedy
// hacer secuencial, empezando desde el mínimo (no es necesario hacer greedy)
// aplicar bitwise
// */


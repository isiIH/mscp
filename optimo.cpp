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

// struct slong{
//     ulong* x;
//     int size;

//     slong(int size) {
//         this->size = size;
//         this->x = new ulong[size];
//     }

//     void set(int i, ulong n){
//         x[i] = n;
//     }

//     ulong operator[] (int i) {
//         return x[i];
//     }

//     slong operator& (slong n) // slong & slong -> slong
//     {
//         slong a(this->size);
//         for(int i=0; i<size; i++) a.set(i,x[i] & n[i]);
//         return a;
//     }
    
//     template<typename T> slong operator<< (T n) // slong << 500 -> slong
//     {
//         int i = (ulong)n / 64;
//         int j = (ulong)n % 64;

//         ulong num = x[0];
//         for(int k=0; k<i; k++) x.set(k,0);
//         x.set(i,num<<j);
//         return x;
//     }

//     template<typename T> slong operator+ (T n) // slong << 500 -> slong
//     {
//         if(__builtin_popcount(x[0]) == sizeof(ulong)*8) {}
//         x.set(0,x[0]++);
//     }
// };

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
    vector<ulong*> exh_sol;
    ulong sizeF, sizeNF;
	ulong n, m, nWX, nWF;
} ParProg;

ParProg* par;

void readFile(string filename);

void exhaustive_sol();

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

    //SOL. CLASSIC GREEDY ALGORITHM
    cout << "-------------------------" << endl;
    cout << "Executing new exhaustive algorithm..." << endl;


    auto start_time = chrono::high_resolution_clock::now();
    exhaustive_sol();
    auto end_time = chrono::high_resolution_clock::now();

    if(PRINT) {
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
    string line,item;

    //m & n
	getline(file,line);
    istringstream ss(line);
    ss >> (par->m) >> (par->n);


    par->nWX = (par->n)/(sizeof(ulong)*8);
    if ((par->n)%(sizeof(ulong)*8)>0) par->nWX++;
    par->X = new ulong[par->nWX];

    par->nWF = (par->m)/(sizeof(ulong)*8);
    if ((par->m)%(sizeof(ulong)*8)>0) par->nWF++;

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

void exhaustive_sol() {
    // Calcular mínimo número de subconjuntos
    vector<ulong*> Fsort = par->bF;
    sort(Fsort.begin(), Fsort.end(), [&](ulong* a, ulong* b){return countSet(a) > countSet(b);});
    // printSubsets(Fsort);
    int minSS = 0;
    int k = 0;
    while(minSS < par->n) {
        minSS += countSet(Fsort[k]);
        k++;
    }

    cout << "K = " << k << endl;
    cout << "MinSS = " << minSS << endl;

    ulong* eCover = new ulong[par->nWX];
    ulong* unionC;

    ulong i;
    ulong j;
    int l;

    //Iterar desde K hasta encontrar el óptimo
    while(true) {
        for(i=1;i<((ulong)1<<par->m);i++){
            if (__builtin_popcount(i) == k) {
                vector<ulong*> C;
                for(j=0;j<par->m;j++){
                    if(i&((ulong)1<<j)){
                        // cout << j << endl;
                        C.push_back(Fsort[j]);
                    }
                }

                // cout << "----------------" << endl;
                // printSubsets(C);
                // cout << "----------------" << endl;

                // cout << C.size() << endl;

                unionC = unionF(C);
                for(l=0; l<par->nWX; l++) eCover[l] = unionC[l] & par->X[l];
                
                if(countSet(eCover) == par->n){
                    par->exh_sol = C;
                    return;
                }

                C.clear();
            }
        }

        //No se encuentra una solución de tamaño k, aumentamos en 1
        k++;
    }
}

ulong* unionF(vector<ulong*> &F) {
    int i,j;
    ulong* C = new ulong[par->nWX];

    for(i=0; i<F.size(); i++) {
        for(j=0; j<par->nWX; j++) C[j] = C[j] | F[i][j];
    }

    return C;
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
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

struct slong{
    ulong* x;
    int size;

    slong(int size) {
        this->size = size;
        this->x = new ulong[size];
        for (int i=0; i<size; i++) this->x[i] = 0;
    }

    void set(int i, ulong n){
        x[i] = n;
    }

    ulong operator[] (int i) {
        return x[i];
    }

    void clear() {
        for (int i=0; i<size; i++) this->x[i] = 0; 
    }

    void nextComb() {
        bool carry = true;
        int i = 0;
        while(carry) {
            x[i]++;
            if(x[i] != 0) carry = false;
            i++;
        }
    }

    slong operator& (slong n) // slong & slong -> slong
    {
        slong a(this->size);
        for(int i=0; i<size; i++) a.set(i,x[i] & n[i]);
        return a;
    }
    
    // template<typename T> slong operator<< (T n) // slong << 500 -> slong
    // {
    //     int i = (ulong)n / 64;
    //     int j = (ulong)n % 64;

    //     ulong num = x[0];
    //     for(int k=0; k<i; k++) x.set(k,0);
    //     x.set(i,num<<j);
    //     return x;
    // }



    // template<typename T> slong operator+ (T n) // slong + n -> slong
    // {
    //     ulong i = 0;
    //     bool carry = true;
    //     while(carry) {
    //         if(getBit64(x, i)) {
    //             cleanBit64(x, i);
    //             i++;
    //         } else {
    //             carry = false;
    //             setBit64(x,i);
    //         }
    //     }
    // }
};

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
    vector<item> mp;
    vector<ulong*> exh_sol;
    ulong sizeF, sizeNF;
	ulong n, m, nWX, nWF;
} ParProg;

ParProg* par;

void readFile(string filename);
void preprocess();
void createMap();

void exhaustive_sol();

ulong* unionF(vector<ulong*> &F);
int countSet(ulong* S, int size=par->nWX);
vector<ulong*> getComb(ulong* comb, vector<ulong*> F);

void printSubset(ulong *S, int size = par->nWX);
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

    auto start_time = chrono::high_resolution_clock::now();
    preprocess();
    auto end_time = chrono::high_resolution_clock::now();
	cout << "Duración en milisegundos: " << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000.0 << endl;


    //SOL. CLASSIC GREEDY ALGORITHM
    cout << "-------------------------" << endl;
    cout << "Executing new exhaustive algorithm..." << endl;

    start_time = chrono::high_resolution_clock::now();
    exhaustive_sol();
    end_time = chrono::high_resolution_clock::now();

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

    // cout << "K = " << k << endl;
    // cout << "MinSS = " << minSS << endl;

    ulong* eCover = new ulong[par->nWX];
    ulong* unionC;

    ulong i;
    int j;

    //Iterar desde K hasta encontrar el óptimo

    slong setsComb = slong(par->nWF);
    vector<ulong*> C;
    while(true) {
        cout << "K = " << k << endl;
        while(countSet(setsComb.x, par->nWF) < par->m) {
            setsComb.nextComb();
            if( countSet(setsComb.x, par->nWF) == k ) {
                // printSubset(setsComb.x, par->nWF);
                C = getComb(setsComb.x,Fsort);
                // printSubsets(C);
                // cout << "----------------------" << endl;

                unionC = unionF(C);
                for(j=0; j<par->nWX; j++) eCover[j] = unionC[j] & par->X[j];
                
                if(countSet(eCover) == par->n){
                    cout << "Solution found with K = " << k << endl;
                    for(ulong* S : C) par->exh_sol.push_back(S);
                    return;
                }
                C.clear();
            }
        }

        //No se encuentra una solución de tamaño k, aumentamos en 1
        k++;
        setsComb.clear();
    }
}

void preprocess() {
    cout << "-------------------------" << endl;
    cout << "Executing PreSetCover..." << endl;
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
        par->exh_sol.push_back(S);
    }

    for(ulong* S : par->exh_sol) {
        par->bF.erase(find(par->bF.begin(), par->bF.end(), S));
        par->m--;
    }

    cout << "Added " << par->exh_sol.size() << " subsets " << endl; 
    par->n = countSet(par->X);
    cout << "|X| = " << par->n << endl;
    cout << "|F| = " << par->m << endl;

    par->nWF = (par->m)/(sizeof(ulong)*8);
    if ((par->m)%(sizeof(ulong)*8)>0) par->nWF++;

    cout << "nWF = " << par->nWF << endl;
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

    for(i=0; i<F.size(); i++) {
        for(j=0; j<par->nWX; j++) C[j] = C[j] | F[i][j];
    }

    return C;
}

int countSet(ulong* S, int size){
    int cont = 0;
    for(int i=0; i<size; i++) {
        uint cnt = 0;
        ulong mask = 0x8000000000000000;

        for(cnt=0;cnt<W64;++cnt){
            if((S[i] & mask) != 0) cont++;
            mask >>= 1;
        }
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

vector<ulong*> getComb(ulong* comb, vector<ulong*> F) {
    vector<ulong*> C;

    for(int i=0; i<par->nWF; i++) {
        uint cnt = 0;
        ulong mask = 1;

        for(cnt=0;cnt<W64;++cnt){
            if((comb[i] & mask) != 0) C.push_back(F[cnt]);
            mask <<= 1;
        }
    }

    return C;
}

void printSubsets(vector<ulong*> C) {
    for(ulong* S : C) {
        printSubset(S);
    }
}
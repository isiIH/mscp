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
#define MAX_F_SIZE 32

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
    vector<ulong*> unique_elements;
    vector<ulong*> greedy_sol;
    vector<ulong*> greedy2_sol;
    ulong sizeF, sizeNF;
	ulong n, m, nWX;
} ParProg;

ParProg* par;

void readFile(string filename);
void preprocess();

vector<ulong*> greedy(const ulong* X, const vector<ulong*> &F);
void greedyExh();
void optimalSol(int i, const ulong* X, const vector<ulong*> &F, vector<ulong*> chosenSets, vector<ulong*> &minSetCover, int &minSetSize);
void createMap();
ulong* setsOfLength(const vector<item> &mp, const int l, vector<ulong*> &C);

bool isCovered(const vector<ulong*> subsets, const ulong* X);
int intersectionLength(const ulong* A, const ulong* B);
ulong* unionF(const vector<ulong*> &F);
int countSet(const ulong* S);

void printSubset(const ulong* C);
void printSubsets(const vector<ulong*> &C);

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

    auto start_time = chrono::high_resolution_clock::now();
    // preprocess();
    auto end_time = chrono::high_resolution_clock::now();
	cout << "Duración en microsegundos: " << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count() << endl;

    //SOL. CLASSIC GREEDY ALGORITHM
    start_time = chrono::high_resolution_clock::now();
    par->greedy_sol = greedy(par->X, par->bF);
    end_time = chrono::high_resolution_clock::now();
    cout << isCovered(par->greedy_sol, par->X) << endl;

    if(CHECK) {
        cout << "SOL: " << endl;
        printSubsets(par->greedy_sol);
    }
    cout << "Solution Cardinality: " << par->greedy_sol.size() << endl;

    cout << "Duración en microsegundos: " << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count() << endl;

    //SOL. NEW GREEDY ALGORITHM
    createMap();
    start_time = chrono::high_resolution_clock::now();
    greedyExh();
    end_time = chrono::high_resolution_clock::now();

    if(CHECK) {
        cout << "SOL: " << endl;
        printSubsets(par->greedy2_sol);
    }
    cout << "Solution Cardinality: " << par->greedy2_sol.size() << endl;

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

vector<ulong*> greedy(const ulong* X, const vector<ulong*> &F) {
    cout << "-------------------------------------" << endl;
    cout << "Executing classic greedy algorithm..." << endl;
    cout << "-------------------------------------" << endl;
    int i;
    ulong* U = new ulong[par->nWX];
    for(i=0; i<par->nWX; i++) U[i] = X[i];
    vector<ulong*> subsets = F;
    vector<ulong*> C;
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

    return C;
}

void greedyExh(){
    cout << "---------------------------------" << endl;
    cout << "Executing new greedy algorithm..." << endl;
    cout << "---------------------------------" << endl;

    ulong* U = new ulong[par->nWX];
    for(int i=0; i<par->nWX; i++) U[i] = par->X[i];
    int j;
    int k = par->mp[0].rep;
    vector<ulong*> subF;
    ulong* subU;
    int minSetSize;
    vector<ulong*> chosenSets;
    vector<ulong*> minSetCover;

    while(countSet(U) > 0) {
        subU = setsOfLength(par->mp,k, subF);

        if(!subF.empty()){

            if(PRINT) {
                cout << "Rep: " << k << endl;
                cout << "U: " << countSet(U) << endl;
                cout << "SubU: " << countSet(subU) << endl;
                cout << "SubF: " << subF.size() << endl;
            }

            if(subF.size() <= MAX_F_SIZE) {
                minSetSize = par->m+1;
                optimalSol(0, subU, subF, chosenSets, minSetCover, minSetSize);
                chosenSets.clear();
            } else {
                minSetCover = greedy(subU, subF);
            }

            if(PRINT) {
                cout << "Sol. Exhaustiva: " << minSetCover.size() << endl;
                cout << "----------------------" << endl;
            }

            for (ulong* S : minSetCover){
                for(j=0; j<par->nWX; j++) U[j] = U[j] & ~S[j];

                par->mp.erase(remove_if(par->mp.begin(), par->mp.end(), [S](const item& it) {
                    return find_if(it.subSets.begin(), it.subSets.end(), [=](const int& set) { return par->bF[set] == S; }) != it.subSets.end();
                }), par->mp.end());

                par->greedy2_sol.push_back(S);
            }

            subF.clear();
            minSetCover.clear();
        }
        k++;
    }

}

void optimalSol(int i, const ulong* X, const vector<ulong*> &F, vector<ulong*> chosenSets, vector<ulong*> &minSetCover, int &minSetSize) {
    //Si ya no hay más conjuntos que agregar o si no hay una mejor solución por esta rama
    if(i == F.size() || chosenSets.size() >= minSetSize-1) return;

    for (int j = i; j<F.size(); j++) {
        //Incluir el subconjunto
        chosenSets.push_back(F[j]);

        //Verificar si se cubre el universo
        if(isCovered(chosenSets, X) & chosenSets.size() < minSetSize) {
            if(PRINT) cout << "NEW MIN = " << chosenSets.size() << endl;
            minSetSize = chosenSets.size();
            minSetCover = chosenSets;
            return;
        }

        optimalSol(j+1, X, F, chosenSets, minSetCover, minSetSize);
        
        //No incluir el subconjunto
        chosenSets.pop_back();
    }
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
            par->mp.erase(remove_if(par->mp.begin(), par->mp.end(), [e](const item& mp) {return mp.value == e;}), par->mp.end());
            cleanBit64(par->X,e-1);
        }
        
        par->greedy_sol.push_back(S);
        par->greedy2_sol.push_back(S);
    }

    cout << "Added " << par->greedy_sol.size() << " subsets " << endl; 
    cout << "|X| = " << countSet(par->X) << endl;

    if(CHECK) {
        cout << "----------------------" << endl;
        for(item mp_item : par->mp) {
            cout << " - " << mp_item.value << " - " << endl;
            cout << mp_item.rep << " subsets." << endl;
            for (int setIndex : mp_item.subSets) printSubset(par->bF[setIndex]);
        }
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

ulong* setsOfLength(const vector<item> &mp, const int l, vector<ulong*> &C) {
    ulong* subU= new ulong[par->nWX];
    fill(subU, subU + par->nWX, 0);
    int i=0;
    while (i < mp.size() && mp[i].rep == l) {
        for (int index : mp[i].subSets){
            ulong* S = par->bF[index];
            if(find(C.begin(), C.end(), S) == C.end()) {
                C.push_back(S);
            }
        }
        setBit64(subU, mp[i].value-1);
        i++;
    }
    
    sort(C.begin(), C.end(), [&](ulong* a, ulong* b){return countSet(a) > countSet(b);});

    return subU;
}

bool isCovered(const vector<ulong*> subsets, const ulong* X) {
    ulong* coveredElements = unionF(subsets);
    for (int i = 0; i < par->nWX; i++) if ((coveredElements[i] & X[i]) != X[i]) {
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
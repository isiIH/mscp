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
#define MAX_F_SIZE 32

typedef struct{
	int value;
	int rep;
	vector<int> subSets;
} item;

// Structure with all globals parameters program
typedef struct {
	vector<set<int>> F;
    set<int> chi;
    map<int,int> elem_pos;
	ulong* X;
    vector<ulong*> bF;
    vector<item> mp;
    vector<ulong*> greedy_sol;
    vector<ulong*> greedy_exh_sol;
    ulong sizeF, sizeNF;
	ulong n, m, nWX;
} ParProg;

ParProg* par;

void readFile(string filename);
void readFileScp(string filename);
void readFilePartition(string filename);
void analizeF();
void preprocess();

vector<ulong*> greedy(const ulong* X, const vector<ulong*> &F);
void greedyExh();
vector<ulong*> findGreedySets(ulong* elems, const int k);
void optimalSol(int i, const ulong* elems, const vector<ulong*> &F, vector<ulong*> chosenSets, vector<ulong*> &minSetCover, int &minSetSize, int &maxCover);
ulong* setsOfLength(const int l, vector<ulong*> &C);

bool isCovered(const ulong* sumF, const ulong* X);
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
    auto start_time = chrono::high_resolution_clock::now();
    analizeF();
    auto end_time = chrono::high_resolution_clock::now();
    auto dur_analyze = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();

    if(PRINT) {
        cout  << "X: " << par->n
        << " | F: " << par->m << endl;
    }

    par->sizeF = par->m*sizeof(ulong)*par->n;
    par->sizeNF = par->m*sizeof(ulong)*par->nWX;

	cout << "nWX = " << par->nWX << endl;
	cout << " size for F[] = " << par->sizeF/(1024.0*1024.0) << " MiB" << endl;
	cout << " size for nF[] = " << par->sizeNF/(1024.0*1024.0) << " MiB" << endl;

    if(CHECK) {
        for(set<int> sub : par->F) {
            for(int val : sub) {
                cout << val << " - ";
            }
            cout << endl;
        }

        printSubsets(par->bF);
    }

    //SOL. CLASSIC GREEDY ALGORITHM
    start_time = chrono::high_resolution_clock::now();
    par->greedy_sol = greedy(par->X, par->bF);
    end_time = chrono::high_resolution_clock::now();
    auto dur_greedy = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();

    if(CHECK) {
        cout << "SOL: " << endl;
        printSubsets(par->greedy_sol);
    }

    //PREPROCESS
    start_time = chrono::high_resolution_clock::now();
    preprocess();
    end_time = chrono::high_resolution_clock::now();
    auto dur_preprocess = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();

    //SOL. NEW GREEDY ALGORITHM
    start_time = chrono::high_resolution_clock::now();
    greedyExh();
    end_time = chrono::high_resolution_clock::now();
    auto dur_greedyExh = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();
    dur_greedyExh += dur_preprocess + dur_analyze;

    if(CHECK) {
        cout << "SOL: " << endl;
        printSubsets(par->greedy_exh_sol);
    }
    cout << "Greedy Cardinality: " << par->greedy_sol.size() << endl;
    cout << "Time [μs]: " << dur_greedy << endl;
    cout << "Greedy Exhaustive Cardinality: " << par->greedy_exh_sol.size() << endl;
    cout << "Time [μs]: " << dur_greedyExh << endl;

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
                par->F[stoi(item)-1].insert(i+1);
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
    set<int> sub;
    for (int i = 0; i < par->m; i++) {
        getline(file>>std::ws,line);
        istringstream ss(line);
        getline(ss>>std::ws, item, ' ');
        getline(ss>>std::ws, item, ' ');

        while (getline(ss>>std::ws, item, ' ')) {
            sub.insert(stoi(item));
        }
        (par->F).push_back(sub);
        sub.clear();
    }
    file.close();
}

void analizeF() {
    unordered_map<int, vector<int>> inSet;
    unordered_map<int,int> cont;
    for( int i=0; i<par->F.size(); i++ ) {
        for( int e : par->F[i] ) {
            par->chi.insert(e);
            cont[e]++;
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
    for(pair<int, int> values : cont){
        setBit64(par->X, pos);
        par->elem_pos[values.first] = pos;
        par->mp[pos].value = values.first;
        par->mp[pos].subSets = inSet[values.first];
        par->mp[pos].rep = values.second;
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

vector<ulong*> greedy(const ulong* X, const vector<ulong*> &F) {
    if(PRINT) {
        cout << "-------------------------------------" << endl;
        cout << "Executing classic greedy algorithm..." << endl;
        cout << "-------------------------------------" << endl;
    }
    int i;
    ulong* U = new ulong[par->nWX];
    for(i=0; i<par->nWX; i++) U[i] = X[i];
    vector<ulong*> subsets = F;
    vector<ulong*> C;
    int posSet;
    int maxLengthSS = 0;
    int lengthSS;

    while( countSet(U) > 0 ) {
        if(CHECK) cout << "POSIBLE SETS:" << endl;
        for( i=0; i<subsets.size(); i++ ){
            lengthSS = intersectionLength(U, subsets[i]);
            if(lengthSS > maxLengthSS) {
                maxLengthSS = lengthSS;
                posSet = i;
            }
            if(CHECK) {
                for( pair<int, int> values : par->elem_pos ) if(getBit64(subsets[i], values.second)) cout << values.first << " ";
                cout << " | " << lengthSS << endl;
            }
        }

        if(CHECK) {
            cout << "U = " << countSet(U) << ", { ";
            for( pair<int, int> values : par->elem_pos ) if(getBit64(U, values.second)) cout << values.first << " ";
            cout << "}" << endl << "Best Set: { ";
            for( pair<int, int> values : par->elem_pos ) if(getBit64(subsets[posSet], values.second)) cout << values.first << " ";
            cout << "} > " << maxLengthSS << endl;
        }

        for(i=0; i<par->nWX; i++) U[i] = U[i] & ~subsets[posSet][i];
        C.push_back(subsets[posSet]);
        subsets.erase(subsets.begin()+posSet);

        maxLengthSS = 0;
    }

    return C;
}

vector<ulong*> findGreedySets(ulong* elems, const int k) {
    int i;
    vector<ulong*> subsets;
    vector<ulong*> C;
    int posSet;
    int maxLengthSS = 0;
    int lengthSS;

    while( countSet(elems) > 0 ) {
        if(CHECK) cout << "POSIBLE SETS:" << endl;
        setsOfLength(k, subsets);
        for( i=0; i<subsets.size(); i++ ){
            lengthSS = intersectionLength(par->X, subsets[i]);
            if(lengthSS > maxLengthSS) {
                maxLengthSS = lengthSS;
                posSet = i;
            }
            if(CHECK) {
                for( pair<int, int> values : par->elem_pos ) if(getBit64(subsets[i], values.second)) cout << values.first << " ";
                cout << " | " << lengthSS << endl;
            }
        }

        if(CHECK) {
            cout << "X = " << countSet(par->X) << ", { ";
            for( pair<int, int> values : par->elem_pos ) if(getBit64(par->X, values.second)) cout << values.first << " ";
            cout << "}" << endl << "Best Set: { ";
            for( pair<int, int> values : par->elem_pos ) if(getBit64(subsets[posSet], values.second)) cout << values.first << " ";
            cout << "}" <<endl;
        }

        for(i=0; i<par->nWX; i++) {
            par->X[i] = par->X[i] & ~subsets[posSet][i];
            elems[i] = elems[i] & ~subsets[posSet][i];
        }
        C.push_back(subsets[posSet]);

        maxLengthSS = 0;
    }

    return C;
}

void greedyExh(){
    if(PRINT) {
        cout << "---------------------------------" << endl;
        cout << "Executing new greedy algorithm..." << endl;
        cout << "---------------------------------" << endl;
    }

    int j;
    int k = par->mp[0].rep;
    vector<ulong*> subF;
    ulong* subU;
    ulong* sumF;
    int minSetSize;
    int maxCover;
    vector<ulong*> chosenSets;
    vector<ulong*> minSetCover;

    while(countSet(par->X) > 0) {
        subU = setsOfLength(k, subF);

        if(!subF.empty()){

            if(CHECK) {
                cout << "Universe = { ";
                    for( pair<int, int> values : par->elem_pos ) if(getBit64(par->X, values.second)) cout << values.first << " ";
                    cout << "}" << endl;
                cout << "Sets = " << endl;
                for( ulong* set : subF ){
                    cout << "{ ";
                    for( pair<int, int> values : par->elem_pos ) if(getBit64(set, values.second)) cout << values.first << " ";
                    cout << "}" << endl;
                }
            }

            if(CHECK) {
                cout << "Rep: " << k << endl;
                cout << "U: " << countSet(par->X) << endl;
                cout << "SubU: " << countSet(subU) << endl;
                cout << "SubF: " << subF.size() << endl;
            }

            if(subF.size() <= MAX_F_SIZE) {
                sort(subF.begin(), subF.end(), [&](ulong* a, ulong* b){return countSet(a) > countSet(b);});
                minSetSize = par->m+1;
                maxCover = 0;
                optimalSol(0, subU, subF, chosenSets, minSetCover, minSetSize, maxCover);
                chosenSets.clear();
            } else {
                minSetCover = findGreedySets(subU, k);
            }

            if(CHECK) {
                cout << "Sol. Exhaustiva: " << minSetCover.size() << endl;
                for( ulong* set : minSetCover ){
                    cout << "{ ";
                    for( pair<int, int> values : par->elem_pos ) if(getBit64(set, values.second)) cout << values.first << " ";
                    cout << "}" << endl;
                }
                cout << "----------------------" << endl;
            }

            for (ulong* S : minSetCover){
                for(j=0; j<par->nWX; j++) par->X[j] = par->X[j] & ~S[j];

                par->mp.erase(remove_if(par->mp.begin(), par->mp.end(), [S](const item& it) {
                    return find_if(it.subSets.begin(), it.subSets.end(), [=](const int& set) { return par->bF[set] == S; }) != it.subSets.end();
                }), par->mp.end());

                par->greedy_exh_sol.push_back(S);
            }

            subF.clear();
            minSetCover.clear();
        }
        k++;
    }

}

void optimalSol(int i, const ulong* elems, const vector<ulong*> &F, vector<ulong*> chosenSets, vector<ulong*> &minSetCover, int &minSetSize, int &maxCover) {
    //Si ya no hay más conjuntos que agregar o si no hay una mejor solución por esta rama
    if(i == F.size() || chosenSets.size() >= minSetSize) return;

    for (int j = i; j<F.size(); j++) {
        //Incluir el subconjunto
        chosenSets.push_back(F[j]);

        //Verificar si se cubre el universo
        ulong* sumF = unionF(chosenSets);
        int coveredElements = intersectionLength(par->X, sumF);
        if(isCovered(sumF, elems) && (chosenSets.size() < minSetSize || coveredElements > maxCover)) {
            if(CHECK) cout << "NEW MIN = " << chosenSets.size() << " maxCover = " << coveredElements << endl;
            minSetSize = chosenSets.size();
            maxCover = coveredElements;
            minSetCover = chosenSets;
            delete[] sumF;
        } else {
            delete[] sumF;
            optimalSol(j+1, elems, F, chosenSets, minSetCover, minSetSize, maxCover);
        }
        
        //No incluir el subconjunto
        chosenSets.pop_back();
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

        par->greedy_exh_sol.push_back(S);
    }

    if(PRINT) {
        cout << "Added " << par->greedy_exh_sol.size() << " subsets " << endl; 
        cout << "|X| = " << countSet(par->X) << endl;
    }

    if(CHECK) {
        cout << "----------------------" << endl;
        for(item mp_item : par->mp) {
            cout << " - " << mp_item.value << " - " << endl;
            cout << mp_item.rep << " subsets." << endl;
            for (int setIndex : mp_item.subSets) cout << setIndex << " ";
            cout << endl;
            // for (int setIndex : mp_item.subSets) printSubset(par->bF[setIndex]);
        }
    }
}

ulong* setsOfLength(const int l, vector<ulong*> &C) {
    C.clear();
    ulong* subU= new ulong[par->nWX];
    fill(subU, subU + par->nWX, 0);
    int i=0;
    while (i < par->mp.size() && par->mp[i].rep == l) {
        if(checkBit(par->X, par->elem_pos[par->mp[i].value])) {
            for (int index : par->mp[i].subSets){
                ulong* S = par->bF[index];
                if(find(C.begin(), C.end(), S) == C.end()) {
                    C.push_back(S);
                }
            }
            setBit64(subU, par->elem_pos[par->mp[i].value]);
        }
        i++;
    }

    // sort(C.begin(), C.end(), [&](ulong* a, ulong* b){return countSet(a) > countSet(b);});

    return subU;
}

bool isCovered(const ulong* sumF, const ulong* X) {
    for (int i = 0; i < par->nWX; i++) if ((sumF[i] & X[i]) != X[i]) return false;
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
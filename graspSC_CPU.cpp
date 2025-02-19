#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "../include/BasicCDS.h"
#include <cmath>
#include <set>
#include <map>
#include <unordered_map>
#include <numeric>
#include <cassert>
#include <omp.h>
#include <limits.h>
#include <execution>

using namespace std;
using namespace cds;

#define PRINT 1
#define CHECK 0

// Parámetros
#define RCL 0.7
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

    vector<item> mp;

    vector<int> unique_elements;
    vector<int> greedy_sol;
    vector<int> aprox_sol;

    ulong sizeF, sizeNF;
	ulong n, m, nWX;
    int nt;

    int function;
    bool improve;
    vector<float> worst_columns;
    vector<int> rep_colums;

} ParProg;

ParProg* par;

void readFile(string filename);
void readFileScp(string filename);
void readFilePartition(string filename);
void analyzeF();
void preprocess();

void greedy();

double jaccard(const ulong* A, const ulong* B);
vector<int> graspSC();
vector<int> randSuccintSC(ulong* U, vector<int> init_sol, bool r);

bool isCovered(vector<int> S);
ulong* unionSets(const vector<int> &S);
int countSet(const ulong* S);
int intersectionLength(const ulong* A, const ulong* B);

void printSubset(const ulong *S);
void printSubsets(const vector<ulong*> &C);

int main(int argc, char** argv) {

    if(argc !=4){
		cout << "./opt <filename> <nt> <seed>" << endl;
		exit(EXIT_FAILURE);
	}

    par = new ParProg();

    par->nt = atoi(argv[2]);
    omp_set_num_threads(par->nt);

    // srand(atoi(argv[3]));
    srand(time(0));

    readFile(argv[1]);
    double start_time = omp_get_wtime();
    analyzeF();
    double end_time = omp_get_wtime();
    double dur_analyze = end_time - start_time;

    cout << dur_analyze << endl;
    return 0;

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
    start_time = omp_get_wtime();
    greedy();
    end_time = omp_get_wtime();
    auto dur_greedyExh = end_time - start_time;
    dur_greedyExh += dur_analyze;

    //GRASP
    double dur_apr;
    vector<int> sol;
    int best_card = INT_MAX;
    for(int i=0; i<1; i++) {
        start_time = omp_get_wtime();
        sol = graspSC();
        end_time = omp_get_wtime();
        double time = end_time - start_time;
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
        cout << "Time [s]: " << dur_greedyExh << endl;
        cout << "GraspSC Cardinality: " << par->aprox_sol.size() << endl;
        cout << "Time [s]: " << dur_apr << endl;
    }

    assert(isCovered(par->aprox_sol) && "Solución inválida");

    cout << argv[1] << " " << par->n << " " << par->m << " " << dur_greedyExh << " " << par->greedy_sol.size() << " " << dur_apr << " " << par->aprox_sol.size() << " " << endl;

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
    par->nWX = (par->n)/(sizeof(ulong)*8);
    if ((par->n)%(sizeof(ulong)*8)>0) par->nWX++;
    
    par->X = new ulong[par->nWX];
    fill(par->X, par->X + par->nWX, 0);

    par->mp = vector<item>(par->n);

    for(int i=0; i<par->m; i++){
        ulong *bset = new ulong[par->nWX];
        fill(bset, bset + par->nWX, 0);

        for(int e : par->F[i]) {
            if(!checkBit(par->X, (e-1))) setBit64(par->X, (e-1));
            par->mp[(e-1)].subSets.push_back(i);
            setBit64(bset, (e-1));
        }

        par->bF.push_back(bset);
    }

    #pragma omp parallel for shared(par)
    for(int i=0; i<par->n; i++) {
        par->mp[i].value = i+1;
        par->mp[i].rep = par->mp[i].subSets.size();
    }

    sort(std::execution::par_unseq, par->mp.begin(), par->mp.end(), [&](item a, item b){return a.rep < b.rep;});

    if(CHECK) {
        cout << "X = " << countSet(par->X) << endl;
        cout << "F = " << par->bF.size() << endl;
    }
}

void greedy() {
    int i;
    ulong* U = new ulong[par->nWX];
    for(i=0; i<par->nWX; i++) U[i] = par->X[i];
    vector<int> C;
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
    preprocess();

    vector<int> best_sol;
    int best_size = INT_MAX;
    par->improve = false;
    vector<vector<int>> local_solutions(par->nt);

    #pragma omp parallel default(none) shared(par, local_solutions, best_sol, best_size)
    {
        ulong* U = new ulong[par->nWX];
        for(int i=0; i<par->nWX; i++) U[i] = par->X[i];
        int th = omp_get_thread_num();

        //Solución inicial
        local_solutions[th] = randSuccintSC(U, par->unique_elements, true);
        printf("Init sol. th%d: %ld\n", th, local_solutions[th].size());

        #pragma omp barrier

        #pragma omp master
        {
            best_sol = *min_element(local_solutions.begin(), local_solutions.end(), 
                    [](const vector<int>& a, const vector<int>& b) {
                        return a.size() < b.size();
                    });
            best_size = best_sol.size();
            printf("BEST INITIAL SIZE: %d\n", best_size);
        }

        vector<int> new_sol, setsRemoved;
        ulong* unionSC;
        int numRemove, col;

        for(int iter=0; iter<MAX_ITER; iter++) {
            //Perturbación
            new_sol = local_solutions[th];
            numRemove = rand() % (int)ceil((new_sol.size()-par->unique_elements.size()) * RCL) + 1;

            for(int i=0; i<numRemove; i++) {
                col = rand()%(new_sol.size()-par->unique_elements.size()) + par->unique_elements.size();
                setsRemoved.push_back(new_sol[col]);
                new_sol.erase(new_sol.begin() + col);
            }

            // Actualizar U y map
            unionSC = unionSets(new_sol);
            for(int ss : setsRemoved)  {
                for(int e : par->F[ss]) {
                    if(!checkBit(unionSC, (e-1)))
                        setBit64(U, (e-1));
                }
            }

            setsRemoved.clear();

            // Nueva solución
            new_sol = randSuccintSC(U, new_sol, false);

            // Eliminar subsets redundantes (que no agregan elementos nuevos)
            int i=par->unique_elements.size();
            vector<int> sol_i;
            while(i < new_sol.size()){
                sol_i = new_sol;
                sol_i.erase(sol_i.begin() + i);
                if(isCovered(sol_i)) {
                    new_sol.erase(new_sol.begin() + i);
                }
                else i++;
            }

            #pragma omp reduction(||:par->improve)
            if(new_sol.size() < local_solutions[th].size()) {
                local_solutions[th] = new_sol;
                par->improve = true;
            } else par->improve = false;

            // if(th == 0) {
            //     printf("ITER %d\nSol. Card: %ld\nBest Local Card: %ld\n", iter, new_sol.size(), local_solutions[th].size());
            // }

            #pragma omp barrier

            #pragma omp master 
            {
                best_sol = *min_element(local_solutions.begin(), local_solutions.end(), 
                    [](const vector<int>& a, const vector<int>& b) {
                        return a.size() < b.size();
                    });
                best_size = best_sol.size();
                printf("It:%d - Best Sol. %d\n", iter, best_size);
            }

            #pragma omp barrier

            // Mejorar la mitad de las soluciones si hay mejora (random)
            if(par->improve && rand()%2 == 0) local_solutions[omp_get_thread_num()] = best_sol;
        }
    }

    return best_sol;
    
}

vector<int> randSuccintSC(ulong* U, vector<int> init_sol, bool r) {
    // par->function = rand() % 4;
    par->function = 0;
    vector<int> C = init_sol;
    int posSet;
    set<int> subsets;
    double coverage;
    double bestCoverage = numeric_limits<double>::max();
    int grade;
    int p;
    // ulong* Ux = new ulong[par->nWX];

    double total = 0;
    vector<pair<int, int>> subsets_coverage;
    double rand_subset;

    while( countSet(U) > 0 ) {
        p = 0;
        while(p < par->mp.size() && !checkBit(U, (par->mp[p].value-1))) p++;
        grade = par->mp[p].rep; 
        // for(int i=0; i<par->nWX; i++) Ux[i] = 0;
        while(p < par->mp.size() && par->mp[p].rep == grade) {
            if(checkBit(U, (par->mp[p].value-1)))
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
                case 0: coverage = 1/coverage; break;
                case 1: coverage = 1/sqrt(coverage); break;
                case 2: coverage = 1/log(1 + coverage); break;
                case 3: coverage = 1/(coverage * coverage); break;
                default: break;
            }

            // coverage *= par->worst_columns[ss];

            if(coverage < bestCoverage) {
                bestCoverage = coverage;
                posSet = ss;
            }
            if(r && !par->improve) {
                total += coverage;
                subsets_coverage.push_back(make_pair(posSet, total));
            }
        }
        if(r && !par->improve && rand() % 100 == 0) {
            if(CHECK) cout << "random set" << endl;
            rand_subset = ((double) rand()) / RAND_MAX;
            for(int i=0; i<subsets_coverage.size(); i++) {
                if(rand_subset <= subsets_coverage[i].second / total) {
                    posSet = subsets_coverage[i].first;
                    break;
                }
            }
            subsets_coverage.clear();
            total = 0;
        }

        for(int i=0; i<par->nWX; i++) U[i] = U[i] & ~par->bF[posSet][i];

        C.push_back(posSet);

        if(CHECK) {
            cout << "Best Coverage: " << bestCoverage << endl;
            cout << "Pos. Subset: " << posSet << endl;
            cout << "|U|: " << countSet(U) << endl;
            printSubset(U);
        }
        bestCoverage = numeric_limits<double>::max();
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
    // Add uniques elements - Row Reduction
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
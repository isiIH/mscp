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
#include <cuda_runtime.h>

using namespace std;
using namespace cds;

#define PRINT 1
#define CHECK 0

#define BS 1024

typedef struct{
	int value;
	int rep;
	vector<int> subSets;
} item;

// Structure with all globals parameters program
typedef struct {
    int search;
    int k = 0;
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
    int nt;
} ParProg;

ParProg* par;

void readFile(string filename);
void readFileScp(string filename);
void readFilePartition(string filename);
void analizeF();
void preprocess();

void linearSearch();
void binarySearch(int l, int r);
void exponentialSearch();
void reverseSearch();

void exhaustive_sol();
void greedy();

__device__ int coefBin(int n, int k);
__global__ void generateCombinationsKernel(int m, int k, ulong* d_combinations, ulong* X, ulong* sol, int nWX, ulong combCount, bool* found);
__device__ bool isCovered(ulong* chosenSets, int k, ulong* X, ulong nWX, ulong* sharedMem);
bool launchKernel(int k);

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
    if(par->search < 0 || par->search > 3){
        cout << "Invalid Search Type!\n0: Secuential Search\n1: Binary Search\n2: Exponential Search\n3: Reverse Search" << endl;
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
    cout << "------------------------" << endl;
    cout << "Greedy Cardinality: " << par->greedy_sol.size() << endl;
    cout << "Time [s]: " << dur_greedyExh/1000000.0 << endl;
	cout << "Optimal Cardinality: " << par->exh_sol.size() << endl;
    cout << "Time [s]: " << dur_opt/1000000.0 << endl;

    return 0;
}

void readFile(string filename) {
    if (filename.substr(0,3) == "scp") readFileScp(filename);
    else readFilePartition(filename);
}

void readFileScp(string filename) {
    cout << "Reading file " << filename << "..." << endl;
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

void checkCudaError(cudaError_t result, const char* msg) {
    if (result != cudaSuccess) {
        cerr << msg << ": " << cudaGetErrorString(result) << endl;
        exit(EXIT_FAILURE);
    }
}

void exhaustive_sol() {
    cout << "-----------------------------------------------------------" << endl;
    cout << "Executing new exhaustive algorithm";
    switch (par->search) {
        case 0: cout << " with sequential search..." << endl; break;
        case 1: cout << " with binary search..." << endl; break;
        case 2: cout << " with exponential search..." << endl; break;
        case 3: cout << " with reverse search..." << endl; break;
    }
    cout << "-----------------------------------------------------------" << endl;

    // Calcular mínimo número de subconjuntos
    vector<ulong*> Fsort = par->bF;
    sort(Fsort.begin(), Fsort.end(), [&](ulong* a, ulong* b){return countSet(a) > countSet(b);});
    int minSS = 0;
    while(minSS < par->n) {
        minSS += countSet(Fsort[par->k]);
        par->k++;
    }

    //Iterar desde K hasta encontrar el óptimo
    int max_k = par->greedy_sol.size() - par->unique_elements.size();
    cout << "Search Range = [" << par->k << " - " << max_k << "]" << endl;

    switch (par->search) {
        case 0: //Búsqueda secuencial
            linearSearch(); break;
        case 1: //Búsqueda binaria
            binarySearch(par->k, max_k); break;
        case 2: //Búsqueda exponencial
            exponentialSearch(); break;
        case 3: //Búsqueda secuencial inversa (greedy--)
            reverseSearch(); break;
    }
    par->exh_sol.insert(par->exh_sol.end(), par->unique_elements.begin(), par->unique_elements.end());
}

bool launchKernel(int k) {
    //Calcular num combinaciones
    ulong combCount = 1;
    for(int i=0; i<k; i++) {
        combCount *= (par->m-i);
        combCount /= (i+1);
    }
    cout << "Comb. count: " << combCount << endl;

    dim3 threadsPerBlock(BS);
    dim3 numBlocks((combCount + BS -1) / BS);
    size_t sharedMemSize = BS * k * par->nWX * sizeof(ulong) + BS * par->nWX * sizeof(ulong);

    cout << "tpb: " << BS << endl;
    cout << "num_blocks: " << (combCount + BS -1) / BS << endl;

    int uSize = par->nWX * sizeof(ulong);
    int fSize = par->m * uSize;
    int solSize = k * uSize;

    cout << "uSize: " << uSize << endl;
    cout << "fSize: " << fSize << endl;
    cout << "solSize: " << solSize << endl;

    // Family of subsets
    ulong* d_comb;
    checkCudaError(cudaMalloc(&d_comb, fSize), "Failed to allocate device memory for d_comb");
    for (int i = 0; i < par->m; i++) {
        checkCudaError(cudaMemcpy(&d_comb[i * par->nWX], par->bF[i], uSize, cudaMemcpyHostToDevice), "Failed to copy data from host to device");
    }

    // Universe
    ulong* d_X;
    checkCudaError(cudaMalloc(&d_X, uSize), "Failed to allocate device memory for d_X");
    checkCudaError(cudaMemcpy(d_X, par->X, uSize, cudaMemcpyHostToDevice), "Failed to copy data from host to device");

    // Solution mscp
    ulong* h_sol = new ulong[k*par->nWX];
    ulong* d_sol;
    checkCudaError(cudaMalloc(&d_sol, solSize), "Failed to allocate device memory for d_sol");

    //found
    bool* d_found;
    checkCudaError(cudaMalloc(&d_found, sizeof(bool)), "Failed to allocate device memory for d_found");
    checkCudaError(cudaMemset(d_found, 0, sizeof(bool)), "Failed to set device memory for d_found");

    generateCombinationsKernel<<<numBlocks, threadsPerBlock, sharedMemSize>>>(par->m, k, d_comb, d_X, d_sol, par->nWX, combCount, d_found);
    checkCudaError(cudaDeviceSynchronize(), "Kernel execution failed");

    checkCudaError(cudaMemcpy(h_sol, d_sol, solSize, cudaMemcpyDeviceToHost), "Failed to copy data from device to host");
    bool h_found;
    checkCudaError(cudaMemcpy(&h_found, d_found, sizeof(bool), cudaMemcpyDeviceToHost), "Failed to copy found flag from device to host");

    if(h_found) {
        for(int i=0; i<k; i++) {
            ulong* ss = new ulong[par->nWX];
            for(int j=0; j<par->nWX; j++) {
                if(CHECK) printSubset(h_sol);
                ss[j] = h_sol[i*par->nWX+j];
            }
            par->exh_sol.push_back(ss);
        }
    }

    cudaFree(d_comb);
    cudaFree(d_sol);
    cudaFree(d_X);
    cudaFree(d_found);

    return h_found;
}

__device__ int coefBin(int n, int k){

    if((2*k) > n) k = n - k;
    int combCount = 1;
    for (int i = 0; i < k; i++) {
        combCount *= (n - i);
        combCount /= (i + 1);
    }
    return combCount;
}

__global__ void generateCombinationsKernel(int m, int k, ulong* d_combinations, ulong* X, ulong* sol, int nWX, ulong combCount, bool* found) {
    extern __shared__ ulong shared_mem[];
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    // tid fuera de rango
    if (tid >= combCount || *found) return;

    int combId = tid;
    int offset = 0;
    ulong* uComb = shared_mem + threadIdx.x * k * nWX; //comb única generada por tid

    // Calcular combinación única por tid
    for (int i = 0; i < k; i++) {
        for (int j = offset; j < m; j++) {
            int countRest = coefBin((m-j-1), (k-i-1)); //num combinaciones restantes luego de agregar 1
            if (combId < countRest) { //valor de j se agrega a la combinación
                for(int l=0; l<nWX; l++) uComb[i*nWX+l] = d_combinations[j*nWX+l];
                offset = j + 1;
                break;
            } else { //j++, disminución de num combinaciones
                combId -= countRest;
            }
        }
    }

    //Verificar si la comb cubre el universo
    if(!*found && isCovered(uComb, k, X, nWX, shared_mem + (blockDim.x * k * nWX) + (threadIdx.x * nWX))) {
        printf("tid = %d found a solution!\n", tid);

        //Copiar comb a sol
        if (atomicCAS((int*)found, 0, 1) == 0) {
            printf("tid = %d copying the uComb!\n", tid);
            // Copy comb to sol
            for (int i = 0; i < k; i++) {
                for (int j = 0; j < nWX; j++) sol[i * nWX + j] = uComb[i * nWX + j];
            }
        }

    }
}

__device__ bool isCovered(ulong* chosenSets, int k, ulong* X, ulong nWX, ulong* sharedMem) {
    // Calcular unión de chosenSets
    ulong* C = sharedMem;
    for (int i = 0; i < nWX; i++) C[i] = 0;
    for(int i=0; i<k; i++) for(int j=0; j<nWX; j++) C[j] |= chosenSets[i*nWX+j];

    // Calcular número de conjuntos cubiertos
    bool isCov = true;
    for (int i = 0; i < nWX; i++) if ((C[i] & X[i]) != X[i]) {
        isCov = false;
        break;
    }
    return isCov;
}

void linearSearch() {
    bool found = false;
    while(!found) {
        if(PRINT) cout << "K = " << par->k << endl;
        auto start = chrono::high_resolution_clock::now();
        found = launchKernel(par->k);
        auto end = chrono::high_resolution_clock::now();

        auto time = std::chrono::duration_cast<chrono::microseconds>(end - start).count();
        if(PRINT) cout << "Time: " << time / 1000000.0 << "[s]" << endl;
        //No se encuentra una solución de tamaño k, aumentamos en 1
        par->k++;
    }
}

void binarySearch(int l, int r) {
    int m = l;
    bool found;
    while(l <= r) {
        if(PRINT) cout << "K = " << m << endl;
        auto start = chrono::high_resolution_clock::now();
        found = launchKernel(m);
        auto end = chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<chrono::microseconds>(end - start).count();
        if(PRINT) cout << "Time: " << time << "[s]" << endl;
        if (found) r = m - 1;
        else l = m + 1;
        m = l + (r - l)/2;
    }
}

void exponentialSearch() {
    int exp = 1;
    int greedySize = par->greedy_sol.size() - par->unique_elements.size();
    bool found = false;

    while(par->k <= greedySize && !found) {
        cout << "K = " << par->k << endl;
        auto start = chrono::high_resolution_clock::now();
        found = launchKernel(par->k);
        auto end = chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<chrono::microseconds>(end - start).count();
        if(PRINT) cout << "Time: " << time / 1000000.0 << "[s]" << endl;
        
        if(!found){
            par->k += exp;
            exp *= 2;
        }
    }

    //Realizar búsqueda binaria en un rango más pequeño
    int l = par->k - exp/2 + 1;
    int r = min(par->k-1, greedySize);
    cout << "Search range for binary search: [" << l << " - " << r << "]" << endl;
    binarySearch(l, r);
}

void reverseSearch() {
    bool found = true;
    int k = par->greedy_sol.size() - par->unique_elements.size();
    while(k >= par->k && found) {
        if(PRINT) cout << "K = " << k << endl;
        auto start = chrono::high_resolution_clock::now();
        found = launchKernel(k);
        auto end = chrono::high_resolution_clock::now();
        auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
        if(PRINT) cout << "Time: " << time / 1000000.0 << "[s]" << endl;
        //No se encuentra una solución de tamaño k, disminuimos en 1
        k--;
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
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <chrono>

// Si -> ulong
// U = U || Si

using namespace std;

#define PRINT 0
#define CHECK 0

typedef struct{
	int value;
	int rep;
	vector<vector<int>> subSets;
} item;

void readFile(string filename, vector<int> &X, vector<vector<int>> &F);

vector<vector<int>> greedy(vector<int> &X, vector<vector<int>> &F);

vector<vector<int>> greedy2(vector<int> &X, vector<vector<int>> &F);
vector<vector<int>> exhaustiveSC(vector<int> &X, vector<int> &U, vector<vector<int>> &F);
vector<item> createMap(vector<int> &X, vector<vector<int>> &F);
vector<vector<int>> sets(vector<vector<int>> &F, const int e);
vector<vector<int>> setsOfLength(vector<item> &mp, int &i, const int l);

vector<int> unionF(vector<vector<int>> &F);
int intersectionLength(vector<int> &A, vector<int> &B);
vector<int> intersection(vector<int> &A, vector<int> &B);
vector<int> subtract(vector<int> &A, vector<int> &B);
vector<int> unionSS(vector<int> &A, vector<int> &B);
int countSet(vector<int>& S);

void printSubset(vector<int> &S);
void printSubsets(vector<vector<int>> &F);

int main(int argc, char** argv) {

    if(argc !=2){
		cout << "./greedy <filename>" << endl;
		exit(EXIT_FAILURE);
	}

    vector<int> X;
    vector<vector<int>> F;

    string filename = argv[1];
    readFile(filename, X, F);

    cout  << "X: " << countSet(X)
        << " | F: " << F.size() << endl;

    if(PRINT) {
        cout << endl <<  "SETS: " << endl;
        printSubsets(F);
        cout << endl;
    }

    //SOL. CLASSIC GREEDY ALGORITHM
    cout << "-------------------------" << endl;

    auto start_time = chrono::high_resolution_clock::now();
    vector<vector<int>> sol = greedy(X,F);
    auto end_time = chrono::high_resolution_clock::now();

    if(PRINT) {
        cout << "SOL: " << endl;
        printSubsets(sol);
    }
    cout << "Solution Cardinality: " << sol.size() << endl;

    cout << "Duración en segundos: " << chrono::duration_cast<chrono::seconds>(end_time - start_time).count() << endl;
	cout << "Duración en microsegundos: " << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count() << endl;

    //SOL. NEW GREEDY ALGORITHM
    cout << "-------------------------" << endl;

    start_time = chrono::high_resolution_clock::now();
    sol = greedy2(X,F);
    end_time = chrono::high_resolution_clock::now();

    if(PRINT) {
        cout << "SOL: " << endl;
        printSubsets(sol);
    }
    cout << "Solution Cardinality: " << sol.size() << endl;

    cout << "Duración en segundos: " << chrono::duration_cast<chrono::seconds>(end_time - start_time).count() << endl;
	cout << "Duración en microsegundos: " << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count() << endl;

    return 0;
}

void readFile(string filename, vector<int> &X, vector<vector<int>> &F) {
    cout << "Reading file " << filename << "..." << endl;
    string nametxt = "test/" + filename;
    ifstream file(nametxt.c_str());
    int m, n;
    string line,item;

    //m & n
	getline(file,line);
    istringstream ss(line);
    ss >> m >> n;

    X.insert(X.begin(),n,1);

    //Costs
    int i = 0;
    while(i < n)
    {
        getline(file>>std::ws,line);
        istringstream iss(line);
        while (getline(iss, item, ' ')){i++;}
    }

    //Sets
    int numCover;
    i = 0;
    int j;
    while(i < m)
    {
        getline(file>>std::ws,line);
        numCover = stoi(line);
        vector<int> set(n,0);

        j = 0;
        while(j < numCover){
            getline(file>>std::ws,line);
            istringstream iss(line);
            while (getline(iss, item, ' ')) {
                    set[stoi(item)-1] = 1;
                    j++;
            }
        }

        F.push_back(set);
        i++;
    }
}

vector<vector<int>> greedy(vector<int> &X, vector<vector<int>> &F){
    cout << "Executing classic greedy algorithm..." << endl;
    vector<int> U = X;
    vector<vector<int>> C;
    vector<int> maxS;
    int maxLengthSS = 0;
    int lengthSS;

    while( countSet(U) > 0 ) {

        for( vector<int> S : F ){
            if( find(C.begin(), C.end(), S) == C.end()) {
                lengthSS = intersectionLength(U, S);
                // cout << "S:" << endl;
                // printSubset(S);
                // cout << lengthSS << endl;
                if(lengthSS > maxLengthSS) {
                    maxLengthSS = lengthSS;
                    maxS = S;
                }
            }
        }

        U = subtract(U, maxS);
        // cout << countSet(U) << endl;
        C.push_back(maxS);

        // cout << "maxS: ";
        // printSubset(maxS);
        // cout << "U: ";
        // printSubset(U);
        // cout << "C: ";
        // printSubsets(C);
        // cout << endl;

        maxLengthSS = 0;
        maxS.clear();
    }

    return C;
}

vector<vector<int>> greedy2(vector<int> &X, vector<vector<int>> &F){
    cout << "Executing new greedy algorithm..." << endl;
    vector<int> U = X;
    vector<vector<int>> C;
    vector<item> mp = createMap(X,F);
    sort(mp.begin(), mp.end(), [&](item a, item b){return a.rep < b.rep;});
    // for(item ite : mp) {
    //     cout << "------------" << endl;
    //     cout << ite.value << endl;
    //     printSubsets(ite.subSets);
    //     cout << ite.rep << endl;
    // }

    int i = 0;
    int k = 1;
    vector<vector<int>> subF;
    vector<int> subU;
    vector<vector<int>> solExhaustive;

    while(countSet(U) > 0) {
        subF = setsOfLength(mp,i,k);

        if(!subF.empty()){
            subU = unionF(subF);
            if(countSet(subU) > countSet(U)) {
                subU = intersection(subU, U);
            }
            if(CHECK) {
                cout << "----------------------" << endl;
                cout << "Rep: " << k << endl;
                cout << "U: " << countSet(U) << endl;
                cout << "SubU: " << countSet(subU) << endl;
                cout << "SubF: " << subF.size() << endl;
            }
            // vector<vector<int>> solExhaustive = greedy(subU, subF);
            solExhaustive = exhaustiveSC(U, subU, subF);
            if(CHECK) {
                cout << "Sol. Exhaustiva: " << solExhaustive.size() << endl;
            }

            for (vector<int> S : solExhaustive){
                U = subtract(U, S);
                for(item &it : mp){
                    if(find(it.subSets.begin(), it.subSets.end(), S) != it.subSets.end()) {
                        it.subSets.clear();
                    }
                }
                C.push_back(S);
            }
        }
        k++;
    }


    return C;

}

vector<vector<int>> exhaustiveSC(vector<int> &X, vector<int> &U, vector<vector<int>> &F) {
    int n = countSet(U);
    int m = F.size();
    vector<vector<int>> C;
    int bestLen = m+1;
    int MaxENotCover = 0;
    vector<vector<int>> bestC;

    vector<int> unionC;
    int xCover;
    int eNotCover;

    for(long long int i=1;i<(1<<m);i++){
        for(long long int j=0;j<m;j++){
            if(i&(1<<j)){
                C.push_back(F[j]);
            }
        }

        // cout << "sets" << endl;
        // printSubsets(C);
        // cout << "--------------" << endl;

        unionC = unionF(C);
        xCover = intersectionLength(unionC, U);
        eNotCover = intersectionLength(unionC, X);
        // cout << xCover << endl;

        if(xCover == n && C.size() <= bestLen && eNotCover > MaxENotCover){
            bestLen = C.size();
            bestC = C;
            // cout << "BESTC" << endl;
            // printSubsets(bestC);
            // cout << "--------------" << endl;
            MaxENotCover = eNotCover;
        }

        C.clear();
    }

    return bestC;
}

vector<item> createMap(vector<int> &X, vector<vector<int>> &F) {
    vector<item> mp(X.size());
    for(int i=0; i<X.size(); i++) {
        mp[i].value = i+1;
        mp[i].subSets = sets(F,i);
        mp[i].rep = mp[i].subSets.size();
    }

    return mp;
}

vector<vector<int>> sets(vector<vector<int>> &F, const int e) {
    vector<vector<int>> C;
    for( vector<int> S : F ) {
        if(S[e] == 1) {
            C.push_back(S);
        }
    }

    return C;
}

vector<vector<int>> setsOfLength(vector<item> &mp, int &i, const int l) {
    vector<vector<int>> C;
    while (mp[i].rep == l) {
        for (vector<int> S : mp[i].subSets){
            if(find(C.begin(), C.end(), S) == C.end()) {
                C.push_back(S);
            }
        }
        i++;
    }

    return C;
}

// vector<vector<int>> exhaustiveSC(vector<int> &X, vector<vector<int>> &F) {
//     int n = countSet(X);
//     int m = F.size();
//     vector<vector<int>> C;
//     int bestLen = m+1;
//     vector<vector<int>> bestC;

//     for(long long int i=1;i<(1<<m);i++){
//         for(int j=0;j<m;j++){
//             if(i&(1<<j)){
//                 C.push_back(F[j]);
//             }
//         }

//         vector<int> unionC = unionF(C);
//         int eCover = intersectionLength(unionC, X);
        
//         if(eCover == n && C.size() <= bestLen){
//             bestLen = C.size();
//             bestC = C;
//         }

//         C.clear();
//     }

//     return bestC;
// }

vector<int> unionF(vector<vector<int>> &F) {
    vector<int> C(F[0].size(),0);
    for( vector<int> &S : F ) {
        for(int i=0; i<S.size(); i++){
            if(S[i] == 1){
                C[i] = 1;
            }
        }
    }

    return C;
}

int intersectionLength(vector<int> &A, vector<int> &B) {
    // vector<int> C(A.size(), 0);
    int c = 0;

    for(int i = 0; i<A.size(); i++){
        if(A[i] == 1 && B[i] == 1) {
            // C[i] = 1;
            c++;
        }
    }

    return c;
}

vector<int> intersection(vector<int> &A, vector<int> &B) {
    vector<int> C(A.size(), 0);

    for(int i = 0; i<A.size(); i++){
        if(A[i] == 1 && B[i] == 1) {
            C[i] = 1;
        }
    }

    return C;
}

vector<int> subtract(vector<int> &A, vector<int> &B) {
    vector<int> C(A.size(), 0);

    for(int i = 0; i<A.size(); i++){
        if(A[i] != B[i]) {
            C[i] = A[i];
        }
    }

    return C;
}

vector<int> unionSS(vector<int> &A, vector<int> &B) {
    vector<int> C(A.size(), 0);

    for(int i = 0; i<A.size(); i++){
        if(A[i] == 1 || B[i] == 1) {
            C[i] = 1;
        }
    }

    return C;
}

int countSet(vector<int>& S){
    return count(S.begin(), S.end(), 1);
}

void printSubset(vector<int> &S) {
    for(int i=1; i<=S.size(); i++){
        if(S[i-1] == 1){
            cout << i << " ";
        }
    }
    cout << endl;
}

void printSubsets(vector<vector<int>> &F) {
    for (int i=0; i<F.size(); i++) {
        cout << "Set #" << i+1 << ": ";
        printSubset(F[i]);
        cout << endl;
    }
}

/*
PROYECTO #1: Adaptar solución greedy exhaustiva a bitwise operator

probar con n = 64 primero, luego agregar más celdas

PROYECTO #2: Óptimo que parte de la mímina cantidad posible de subconjuntos
min <= Opt <= Greedy
hacer secuencial, empezando desde el mínimo (no es necesario hacer greedy)
aplicar bitwise
*/


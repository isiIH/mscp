#include <iostream>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>

using namespace std;

typedef struct{
	int value;
	int rep;
	vector<int> subSets;
} item;

class SCP {
public:
    ulong* X;
    vector<vector<int>> F;
    vector<ulong*> bF;

    unordered_map<int, vector<int>> inSet;
    vector<item> mp;

    vector<int> unique_elements;
    ulong* excluded_subsets;

    ulong n, m, nWX, nWF;

    SCP() {
    }

    void readFile(string filename) {
        if (filename.substr(0,3) == "scp") readFileScp(filename);
        else readFilePartition(filename);
    }

    void readFileScp(string filename) {
        cout << "Reading file " << filename << "..." << endl;
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
        ss >> n >> m;

        //Costs
        i = 0;
        while(i < m)
        {
            getline(file>>std::ws,line);
            istringstream iss(line);
            while (getline(iss>>std::ws, item, ' ')){i++;}
        }

        //Sets
        int numCover;
        int j;
        F.resize(m);
        for(i=0; i<n; i++) {
            getline(file>>std::ws,line);
            numCover = stoi(line);

            j = 0;
            while(j < numCover){
                getline(file>>std::ws,line);
                istringstream iss(line);
                while (getline(iss, item, ' ')) {
                    F[stoi(item)-1].push_back(i+1);
                    j++;
                }
            }
        }
        file.close();
    }

    void readFilePartition(string filename) {
        cout << "Reading file " << filename << "..." << endl;
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
        ss >> n >> m;

        //Sets
        vector<int> sub;
        for (int i = 0; i < m; i++) {
            getline(file>>std::ws,line);
            istringstream ss(line);
            getline(ss>>std::ws, item, ' ');
            getline(ss>>std::ws, item, ' ');

            while (getline(ss>>std::ws, item, ' ')) {
                sub.push_back(stoi(item));
            }
            F.push_back(sub);
            sub.clear();
        }
        file.close();
    }

    void analyzeF() {
        nWX = n/(sizeof(ulong)*8);
        if (n%(sizeof(ulong)*8)>0) nWX++;

        nWF = m/(sizeof(ulong)*8); 
        if(m%(sizeof(ulong)*8) > 0) nWF++;

        excluded_subsets = new ulong[nWF];
        for(int i=0; i<nWF; i++) excluded_subsets[i] = 0;
        
        X = new ulong[nWX];
        fill(X, X + nWX, 0);

        ulong *bset;
        for(int i=0; i<F.size(); i++){
            bset = new ulong[nWX];
            fill(bset, bset + nWX, 0);

            for(int e : F[i]) {
                setBit64(X, (e-1));
                inSet[e].push_back(i);
                setBit64(bset, (e-1));
            }

            bF.push_back(bset);
        }

        if(CHECK) {
            cout << "X = " << countSet(X, nWX) << endl;
            cout << "F = " << bF.size() << endl;
        }
    }

    bool isCovered(const vector<int> &S) {
        ulong* coveredElements = unionSets(S);
        for (int i = 0; i < nWX; i++) if ((coveredElements[i] & X[i]) != X[i]) {
            delete[] coveredElements;
            return false;
        }
        delete[] coveredElements;
        return true;
    }

    ulong* unionSets(const vector<int> &S) {
        ulong* C = new ulong[nWX];
        fill(C, C + nWX, 0);
        for(const int s_idx : S) for(int i=0; i<nWX; i++) C[i] |= bF[s_idx][i];
        return C;
    }

    int intersectionLength(const ulong* A, const ulong* B) {
        int cont = 0;
        for(int i=0; i<nWX; i++) cont += __builtin_popcountl(A[i] & B[i]);
        return cont;
    }

    int countSet(const ulong* S, const int word){
        int cont = 0;
        for(int i=0; i<word; i++) {
            cont += __builtin_popcountl(S[i]);
        }
        return cont;
    }

    void printSubset(const ulong *S) {
        for (int i=0; i<nWX; i++){
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

};
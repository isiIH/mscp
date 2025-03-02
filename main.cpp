#include <chrono>
#include <cassert>
#include <map>

using namespace std;

#include <config.h>
#include <Set.h>
#include <SCP.h>
#include <SetCover.h>
#include <Grasp.h>

vector<int> greedy(const SCP &scp);

int main(int argc, char** argv) {
    if(argc !=3){
		cout << "./opt <filename> <seed>" << endl;
		exit(EXIT_FAILURE);
	}

    // srand(atoi(argv[2]));
    srand(time(0));

    // Read file
    SCP scp = SCP(argv[1]);

    // Create succint array
    auto start_time = chrono::high_resolution_clock::now();
    scp.analyzeF();
    auto end_time = chrono::high_resolution_clock::now();
    auto dur_analyze = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();

    if(PRINT) cout  << "X: " << scp.n << " | F: " << scp.m << endl;

    // Calculate original F size and succint F representation
    ulong sizeF = scp.m*sizeof(ulong)*scp.n;
    ulong sizeNF = scp.m*sizeof(ulong)*scp.nWX;

	if(PRINT) {
        cout << "nWX = " << scp.nWX << endl;
        cout << " size for F[] = " << sizeF/(1024.0*1024.0) << " MiB" << endl;
        cout << " size for nF[] = " << sizeNF/(1024.0*1024.0) << " MiB" << endl;
    }

    if(CHECK) {
        for(vector<int> set : scp.F) {
            for(int val : set) {
                cout << val << " ";
            }
            cout << endl;
        }

        scp.printSubsets(scp.bF);
        scp.X.print();
    }

    // GREEDY-ALG
    start_time = chrono::high_resolution_clock::now();
    vector<int> greedySol = vector<int>(0); // greedy(scp);
    end_time = chrono::high_resolution_clock::now();
    auto dur_greedyExh = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();
    dur_greedyExh += dur_analyze;

    // GRASP-ALG
    Grasp alg(scp);
    double dur_apr;
    int best_card = numeric_limits<int>::max();;
    for(int i=0; i<1; i++) {
        start_time = chrono::high_resolution_clock::now();
        SetCover sol = alg.search();
        end_time = chrono::high_resolution_clock::now();
        auto time = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();
        if(sol.size() < best_card) {
            alg.bestSol = sol;
            best_card = sol.size();
            dur_apr = time;
        }
    }
    dur_apr += dur_analyze;

    if(CHECK) {
        cout << "SOL: { ";
        for(int ss : alg.bestSol.solution) {
            cout << ss << " ";
        }
        cout << "}" << endl;
    }
    if(PRINT) {
        cout << "------------------------" << endl;
        cout << "Greedy Cardinality: " << greedySol.size() << endl;
        cout << "Time [s]: " << dur_greedyExh/1000000.0 << endl;
        cout << "GraspSC Cardinality: " << alg.bestSol.size() << endl;
        cout << "Time [s]: " << dur_apr/1000000.0 << endl;
    }

    assert(alg.scp.X.size() == scp.n);
    assert(alg.bestSol.isCovered());

    cout << argv[1] << " " << scp.n << " " << scp.m << " " << alg.bestSol.uniqueSets.size() << " " << alg.bestSol.excludedSets.size() << " " << dur_greedyExh/1000000.0 << " " << greedySol.size() << " " << dur_apr/1000000.0 << " " << alg.bestSol.size() << " " << endl;

    return 0;
}

vector<int> greedy(const SCP &scp) {
    int i;
    Set U(scp.X);
    vector<int> C;
    int maxLengthSS = 0;
    int lengthSS;
    int posSet;

    map<int, Set> subsets;
    for (i=0; i<scp.bF.size(); i++) subsets[i] = scp.bF[i];

    while( U.size() > 0 ) {

        for(pair<int, Set> ss_pos : subsets){
            lengthSS = U.intersectionLength(ss_pos.second);
            if(lengthSS > maxLengthSS) {
                maxLengthSS = lengthSS;
                posSet = ss_pos.first;
            }
        }

        U.substract(subsets[posSet]);
        C.push_back(posSet);
        subsets.erase(posSet);

        maxLengthSS = 0;
    }

    return C;
}
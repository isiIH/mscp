#include <SCP.h>

SCP::SCP() {}

SCP::SCP(const string filename) {
    readFile(filename);
}

void SCP::readFile(const string filename) {
    if (filename.substr(0,3) == "scp") readFileScp(filename);
    else readFilePartition(filename);
}

void SCP::readFileScp(const string filename) {
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

void SCP::readFilePartition(const string filename) {
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

void SCP::analyzeF() {
    nWX = n/(sizeof(ulong)*8);
    if (n%(sizeof(ulong)*8)>0) nWX++;

    nWF = m/(sizeof(ulong)*8); 
    if(m%(sizeof(ulong)*8) > 0) nWF++;
    
    X = Set(nWX);

    for(int i=0; i<F.size(); i++){
        Set bset(nWX);

        for(int e : F[i]) {
            X.push_back(e-1);
            bset.push_back(e-1);
        }

        bF.push_back(bset);
    }

    if(CHECK) {
        cout << "X = " << X.size() << endl;
        cout << "F = " << bF.size() << endl;
    }
}

void SCP::printSubsets(const vector<Set> &C) {
    for(Set S : C) {
        S.print();
    }
}
#include <SCP.h>

SCP::SCP() {}

SCP::SCP(const string filename) {
    readFile(filename);
}

SCP::SCP(const SCP& other) :
    X(other.X),
    F(other.F),
    bF(other.bF),
    n(other.n),
    m(other.m),
    nWX(other.nWX),
    nWF(other.nWF) {}

SCP& SCP::operator=(const SCP& other) {
    if (this != &other) {
        X = other.X;
        F = other.F;
        bF = other.bF;
        n = other.n;
        m = other.m;
        nWX = other.nWX;
        nWF = other.nWF;
    }
    return *this;
}

void SCP::readFile(const string filename) {
    if (filename.substr(0,3) == "scp") readFileScp(filename);
    else readFilePartition(filename);
}

void SCP::readFileScp(const string filename) {
    // Format:
    // number of rows (m), number of columns (n)
    // the cost of each column c(j),j=1,...,n
    // for each row i (i=1,...,m): the number of columns which cover
    // row i followed by a list of the columns which cover row i
    cout << "Reading file " << filename << "..." << endl;
    string nametxt = "test/" + filename + ".txt";
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
//    The format of these data files is:
//    number of rows, number of columns (n)
//    for each column j (j=1,...,n) in turn:
//       column cost, number of rows covered by j, list of the rows covered by j
    cout << "Reading file " << filename << "..." << endl;
    string nametxt = "test/" + filename + ".txt";
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

// Analyze F to create the universe set X and the bitset representation bF of subsets
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
#include <Set.h>

Set::Set() : nW(0), S(nullptr) {}

Set::Set(const int nW) {
    initialize(nW);
    clear();
}

Set::Set(const Set &X) {
    initialize(X.nW);
    for(int i=0; i<nW; i++) S[i] = X.S[i];
}

Set::~Set() {
    delete[] S;
}

Set& Set::operator=(const Set &X) {
    if(this != &X) {
        Set temp(X);
        swap(S, temp.S);
        swap(nW, temp.nW);
    }

    return *this;
}

void Set::initialize(int nW) {
    this->nW = nW;
    S = new ulong[nW];
} 

void Set::push_back(int i) {
    setBit64(S, i);
}

void Set::erase(int i) {
    cleanBit64(S, i);
}

bool Set::check(int i) {
    return checkBit(S, i);
}

void Set::clear() {
    for(int i=0; i<nW; i++) S[i] = 0;
}

int Set::intersectionLength(const Set B) {
    int cont = 0;
    for(int i=0; i<nW; i++) cont += __builtin_popcountl(S[i] & B.S[i]);
    return cont;
}

void Set::substract(const Set B) {
    for(int i=0; i<nW; i++) S[i] = S[i] & ~B.S[i];
}

int Set::size() {
    int rowsCovered = 0;
    for(int i=0; i<nW; i++)
        rowsCovered += __builtin_popcountl(S[i]);
    return rowsCovered;
}

void Set::print() {
    for (int i=0; i<nW; i++){
        printBitsUlong(S[i]);
        cout << " - ";
    }
    cout << endl;
}
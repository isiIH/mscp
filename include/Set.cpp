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
    if (this == &X) {
        return *this;
    }

    if (S != nullptr) {
        delete[] S;
        S = nullptr;
    }

    nW = X.nW;
    if (nW > 0) {
        S = new unsigned long[nW];
        std::copy(X.S, X.S + nW, S);
    } else {
        S = nullptr;
    }

    return *this;
}

// Initialize the set with nW words
void Set::initialize(const int nW) {
    this->nW = nW;
    S = new ulong[nW];
} 

// Add element i to the set
void Set::push_back(const int i) {
    setBit64(S, i);
}

// Remove element i from the set
void Set::erase(const int i) {
    cleanBit64(S, i);
}

// Check if element i is in the set
bool Set::check(const int i) const {
    return checkBit(S, i);
}

// Remove all elements from the set
void Set::clear() {
    for(int i=0; i<nW; i++) S[i] = 0;
}

// Compute the length of the intersection between two sets
int Set::intersectionLength(const Set &B) const {
    int cont = 0;
    for(int i=0; i<nW; i++) cont += __builtin_popcountl(S[i] & B.S[i]);
    return cont;
}

// Compute the union between two sets
void Set::add(const Set& B) {
    for(int i=0; i<nW; i++) S[i] |= B.S[i];
}

// Compute the difference between two sets
void Set::substract(const Set &B) {
    for(int i=0; i<nW; i++) S[i] = S[i] & ~B.S[i];
}

// Compute the number of elements covered by the set
int Set::size() const {
    int rowsCovered = 0;
    for(int i=0; i<nW; i++)
        rowsCovered += __builtin_popcountl(S[i]);
    return rowsCovered;
}
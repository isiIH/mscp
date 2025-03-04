#ifndef SET_H
#define SET_H

#include <iostream>

#include <config.h>
#include <BasicCDS.h>

using namespace std;
using namespace cds;

class Set {
public:
    ulong* S;
    int nW;

    Set();
    Set(const int nW);
    Set(const Set &X);
    ~Set();
    Set& operator=(const Set &X);

    void initialize(const int nW);
    void push_back(const int i);
    void erase(const int i);
    bool check(const int i) const;
    void clear();
    int intersectionLength(const Set &B);
    void add(const Set& B);
    void substract(const Set &B);
    int size() const;
    void print() const;
};

#endif
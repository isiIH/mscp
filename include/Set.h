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

    void initialize(int nW);
    void push_back(int i);
    void erase(int i);
    bool check(int i);
    void clear();
    int intersectionLength(const Set B);
    void substract(const Set B);
    int size();
    void print();
};

#endif
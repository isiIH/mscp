#ifndef ROW_COVERING_H
#define ROW_COVERING_H

#include <vector>

#include <Set.h>
#include <SCP.h>

using namespace std;

class RowCovering {
public:
    int row; // element to be covered
    vector<int> col_covering; // Columns that cover row i
    int n_columns; // size of col_covering

    RowCovering();
    RowCovering(const Set &ignSets, const vector<Set> &bF, const int row);

    void createRowCovering(const Set &ignSets, const vector<Set> &bF);
    int countIntersection(const vector<int>& B) const;
};

#endif
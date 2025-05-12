#include <RowCovering.h>

RowCovering::RowCovering() {}

RowCovering::RowCovering(const Set &ignSets, const vector<Set> &bF, const int row) : row(row), n_columns(0) {
    createRowCovering(ignSets, bF);
}

void RowCovering::createRowCovering(const Set &ignSets, const vector<Set> &bF) {
    for(int i=0; i<bF.size(); i++) {
        if(!ignSets.check(i) && bF[i].check(row)) { // check if column i cover the row
            n_columns++;
            col_covering.push_back(i);
        }         
    }
}

int RowCovering::countIntersection(const vector<int>& B) const {
    int i = 0, j = 0, count = 0;
    int sizeA = col_covering.size(), sizeB = B.size();
    while (i < sizeA && j < sizeB) {
        if (col_covering[i] < B[j]) {
            i++;
        } else if (B[j] < col_covering[i]) {
            j++;
        } else { // col_covering[i] == B[j]
            count++; i++; j++;
        }
    }
    return count;
}
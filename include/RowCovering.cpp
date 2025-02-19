#include <RowCovering.h>

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
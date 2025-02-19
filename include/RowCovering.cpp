#include <RowCovering.h>

RowCovering::RowCovering(SCP &scp, int row) : row(row), n_columns(0) {
    createRowCovering(scp);
}

void RowCovering::createRowCovering(SCP &scp) {
    for(int i=0; i<scp.m; i++) {
        if(!scp.excludedSets.check(i) && scp.bF[i].check(row)) { // check if column i cover the row
            n_columns++;
            col_covering.push_back(i);
        }         
    }
}
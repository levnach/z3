/*++
Copyright (c) 2017 Microsoft Corporation

Module Name:

    <name>

Abstract:

    <abstract>
    Creates the Hermite Normal Form of a matrix in place.
Author:
    Nikolaj Bjorner (nbjorner)
    Lev Nachmanson (levnach)

Revision History:


--*/
#include "util/lp/numeric_pair.h"
#include "util/ext_gcd.h"
namespace lp {
template <typename M> // M is the matrix type
class hnf {
    M &          m_A;
    vector<mpq>  m_buffer;
    unsigned     m() const { return m_A.row_count(); }
    unsigned     n() const { return m_A.column_count(); }
    void process_row_column(unsigned i, unsigned j){ 
        if (is_zero(m_A[i][j]))
            return;
        mpq p,q,r;
        std::cout << "A[" << i << "][" << i << "] = " << m_A[i][i] << "," << std::endl;
        std::cout << "A[" << i << "][" << j << "] = " << m_A[i][j] << "," << std::endl;
        while (j < n()) {
            mpq aii = m_A[i][i];
            mpq aij = m_A[i][j];
            extended_gcd(aii, aij, r, p, q);
            buffer_p_col_i_plus_q_col_j(p, i, q, j);
            replace_column_j_by_lcom(-aii/r, i, aij/r, j);
            copy_buffer_to_col_i(i);
            std::cout << "r = " << r << ", p = " << p << ", q = " << q << std::endl;
            j++;
        }
        // continue step 3 here
    }
    
    void process_row(unsigned i) {
        for (unsigned j = i + 1; j < n(); j++)
            process_row_column(i, j);
    }
    
    void calculate() {
        std::cout << "working" << std::endl;
        for (unsigned i = 0; i < m_A.row_count(); i++) {
            process_row(i);
        }
    }
    
public:
    hnf(M & A) : m_A(A), m_buffer(A.column_count()) {
        calculate();
    }
};
}

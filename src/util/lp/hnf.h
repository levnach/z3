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
    M &          m_H;
    M            m_U;
    M            m_A_orig;  // debug only
    M            m_U_reverse;  // debug only
    vector<mpq>  m_buffer;
    unsigned     m_m;
    unsigned     m_n;
    

    void buffer_p_col_i_plus_q_col_j_H(const mpq & p, unsigned i, const mpq & q, unsigned j) {
        for (unsigned k = i; k < m_m; k++) {
            m_buffer[k] = p * m_H[k][i] + q * m_H[k][j];
        }
    }

    void buffer_p_col_i_plus_q_col_j_U(const mpq & p, unsigned i, const mpq & q, unsigned j) {
        for (unsigned k = 0; k < m_n; k++) {
            m_buffer[k] = p * m_U[k][i] + q * m_U[k][j];
        }
    }

    void pivot_column_i_to_column_j_H(mpq u, unsigned i, mpq v, unsigned j) {
        lp_assert(is_zero(u * m_H[i][i] + v * m_H[i][j]));
        m_H[i][j] = zero_of_type<mpq>();
        for (unsigned k = i + 1; k < m_m; k ++)
            m_H[k][j] = u * m_H[k][i] + v * m_H[k][j];
                  
    }

    void pivot_column_i_to_column_j_U(mpq u, unsigned i, mpq v, unsigned j) {
        for (unsigned k = 0; k < m_n; k ++)
            m_U[k][j] = u * m_U[k][i] + v * m_U[k][j];
                  
    }

    void copy_buffer_to_col_i_H(unsigned i) {
        for (unsigned k = i; k < m_m; k++)
            m_H[k][i] = m_buffer[k];
    }

    void copy_buffer_to_col_i_U(unsigned i) {
        for (unsigned k = 0; k < m_n; k++)
            m_U[k][i] = m_buffer[k];
    }

    // multiply by (a, b)
    //             (c, d)
    // from the left where i and j are the modified columns
    // the [i][i] = a, and [i][j] = b for the matrix we multiply by
   
    
    void multiply_U_reverse_from_left_by(unsigned i, unsigned j, const mpq & a, const mpq & b, const mpq & c, const mpq d) {
        // the new i-th row goes to the buffer
        for (unsigned k = 0; k < m_n; k++) {
            m_buffer[k] = a * m_U_reverse[i][k] + b * m_U_reverse[j][k];
        }

        // calculate the new j-th row in place
        for (unsigned k = 0; k < m_n; k++) {
            m_U_reverse[j][k] = c * m_U_reverse[i][k] + d * m_U_reverse[j][k];
        }

        // copy the buffer into i-th row
        for (unsigned k = 0; k < m_n; k++) {
            m_U_reverse[i][k] = m_buffer[k];
        }
    }
    
    void handle_column_ij_in_row_i(unsigned i, unsigned j) {
        lp_assert(is_correct());
        const mpq& aii = m_H[i][i];
        const mpq& aij = m_H[i][j];
        mpq p,q,r;
        extended_gcd(aii, aij, r, p, q);
        mpq aii_over_r = aii / r;
        mpq aij_over_r = aij / r;
        
        
        buffer_p_col_i_plus_q_col_j_H(p, i, q, j);
        pivot_column_i_to_column_j_H(- aij_over_r, i, aii_over_r, j);
        copy_buffer_to_col_i_H(i);

        
        buffer_p_col_i_plus_q_col_j_U(p, i, q, j);
        pivot_column_i_to_column_j_U(- aij_over_r, i, aii_over_r, j);
        copy_buffer_to_col_i_U(i);

        // U was multiplied from the right by (p, - aij_over_r)
        //                                    (q, aii_over_r  )
        // We need to multiply U_reverse by   (aii_over_r, aij_over_r)
        //                                    (-q        , p)
        // from the left
        
        multiply_U_reverse_from_left_by(i, j, aii_over_r, aij_over_r, -q, p);
        lp_assert(is_correct());
    }

    void switch_sign_for_column(unsigned i) {
        for (unsigned k = i; k < m_m; k++)
            m_H[k][i].neg();
        for (unsigned k = 0; k < m_n; k++)
            m_U[k][i].neg();

        // switch sign for the i-th row in the reverse m_U_reverse
        for (unsigned k = 0; k < m_n; k++)
            m_U_reverse[i][k].neg();
        
    }
    
    void process_row_column(unsigned i, unsigned j){
        std::cout << "process_row_column (" << i << ", " << j << ")" << std::endl;
        if (is_zero(m_H[i][j]))
            return;
        handle_column_ij_in_row_i(i, j);
        
        lp_assert(m_H == m_A_orig * m_U);
    }

    void replace_column_j_by_j_minus_u_col_i_H(unsigned i, unsigned j, const mpq & u) {
        lp_assert(j < i);
        for (unsigned k = i; k < m_m; k++) {
            m_H[k][j] -= u * m_H[k][i];
        }
    }

    void replace_column_j_by_j_minus_u_col_i_U(unsigned i, unsigned j, const mpq & u) {
        
        lp_assert(j < i);
        for (unsigned k = 0; k < m_n; k++) {
            m_U[k][j] -= u * m_U[k][i];
        }
        // Here we multiply from m_U from the right by the matrix ( 1,  0)
        //                                                        ( -u, 1).
        // To adjust the reverse we multiply it from the left by (1, 0)
        //                                                       (u, 1)

        for (unsigned k = 0; k < m_n; k++) {
            m_U_reverse[i][k] += u * m_U_reverse[j][k];
        }
       
        
    }

    void work_on_columns_less_than_i_in_the_triangle(unsigned i) {
        for (unsigned j = 0; j < i; j++) {
            mpq u = ceil(m_H[i][j]/m_H[i][i]);
            if (is_zero(u))
                return;
            replace_column_j_by_j_minus_u_col_i_H(i, j, u);
            replace_column_j_by_j_minus_u_col_i_U(i, j, u);

        }
    }
    
    void process_row(unsigned i) {
        lp_assert(is_correct());
        for (unsigned j = i + 1; j < m_n; j++) {
            process_row_column(i, j);
        }
        if (i >= m_n) {
            lp_assert(m_H == m_A_orig * m_U);
            return;
        }
        if (is_neg(m_H[i][i]))
            switch_sign_for_column(i);
        work_on_columns_less_than_i_in_the_triangle(i);
        std::cout << "H = " << std::endl;
        m_H.print(std::cout);
        auto product = m_A_orig * m_U;
        std::cout << "m_A_orig * m_U = \n"; product.print(std::cout);
        
        lp_assert(m_H == product);
        lp_assert(is_correct());
    }
    
    void calculate() {
        std::cout << "A orig\n";
        m_H.print(std::cout);
        std::cout << "working" << std::endl;
        for (unsigned i = 0; i < m_m; i++) {
            std::cout << "process_row " << i << std::endl;
            process_row(i);
        }
    }

    void prepare_U() {
        auto & v = m_U.m_data;
        v.resize(m_H.column_count());
        for (auto & row: v)
            row.resize(m_H.column_count());
        for (unsigned i = 0; i < m_U.column_count(); i++)
            m_U[i][i] = 1;

        m_U_reverse = m_U;
        
        lp_assert(m_H == m_A_orig * m_U);
    }

    bool row_is_correct_form(unsigned i) const {
        if (i >= m_n)
            return true;
        const mpq& hii = m_H[i][i];
        if (is_neg(hii))
            return false;
        for (unsigned j = 0; j < i; j++) {
            const mpq & hij = m_H[i][j];
            if (is_pos(hij))
                return false;
            if (- hij >= hii)
                return false;
        }
        
        return true;
    }
    
    bool is_correct_form() const {
        for (unsigned i = 0; i < m_m; i++)
            if (!row_is_correct_form(i))
                return false;
        return true;
    }


    bool is_unit_matrix(const M& u) const {
        unsigned m = u.row_count();
        unsigned n = u.column_count();
        if (m != n) return false;
        for (unsigned i = 0; i < m; i ++)
            for (unsigned j = 0; j < n; j++) {
                if (i == j) {
                    if (one_of_type<mpq>() != u[i][j])
                        return false;
                } else {
                    if (!is_zero(u[i][j]))
                        return false;
                }
            }
        return true;
    }
    
    bool is_correct() const {
        return m_H == m_A_orig * m_U && is_unit_matrix(m_U * m_U_reverse);
    }

    bool is_correct_final() const {
        return is_correct() && is_correct_form();
    }

    
public:
    hnf(M & A) : m_H(A),
                 m_A_orig(A),
                 m_buffer(std::max(A.row_count(), A.column_count())),
                 m_m(m_H.row_count()),
                 m_n(m_H.column_count()) {
        prepare_U();
        calculate();
        lp_assert(is_correct_final());
    }
};
}

/*
  Copyright (c) 2013 Microsoft Corporation. All rights reserved.
  Released under Apache 2.0 license as described in the file LICENSE.

  Author: Lev Nachmanson
*/
#include <utility>
#include <memory>
#include <string>
#include <vector>
#include "util/lp/lu.cpp"
template double lp::dot_product<double, double>(std::vector<double, std::allocator<double> > const&, std::vector<double, std::allocator<double> > const&, unsigned int);
template void lp::lu<double, double>::change_basis(unsigned int, unsigned int);
template lp::lu<double, double>::lu(lp::static_matrix<double, double> const&, std::vector<unsigned int, std::allocator<unsigned int> >&, std::vector<int, std::allocator<int> >&, lp::lp_settings&, std::vector<unsigned int, std::allocator<unsigned int> >&);
template void lp::lu<double, double>::push_matrix_to_tail(lp::tail_matrix<double, double>*);
template void lp::lu<double, double>::replace_column(unsigned int, double, lp::indexed_vector<double>&);
template void lp::lu<double, double>::restore_basis_change(unsigned int, unsigned int);
template void lp::lu<double, double>::solve_Bd(unsigned int, std::vector<double, std::allocator<double> >&, lp::indexed_vector<double>&);
template lp::lu<double, double>::~lu();
template void lp::lu<lp::mpq, lp::mpq>::change_basis(unsigned int, unsigned int);
template void lp::lu<lp::mpq, lp::mpq>::push_matrix_to_tail(lp::tail_matrix<lp::mpq, lp::mpq>*);
template void lp::lu<lp::mpq, lp::mpq>::restore_basis_change(unsigned int, unsigned int);
template void lp::lu<lp::mpq, lp::mpq>::solve_Bd(unsigned int, std::vector<lp::mpq, std::allocator<lp::mpq> >&, lp::indexed_vector<lp::mpq>&);
template lp::lu<lp::mpq, lp::mpq>::~lu();
template void lp::lu<lp::mpq, lp::numeric_pair<lp::mpq> >::change_basis(unsigned int, unsigned int);
template void lp::lu<lp::mpq, lp::numeric_pair<lp::mpq> >::push_matrix_to_tail(lp::tail_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >*);
template void lp::lu<lp::mpq, lp::numeric_pair<lp::mpq> >::restore_basis_change(unsigned int, unsigned int);
template void lp::lu<lp::mpq, lp::numeric_pair<lp::mpq> >::solve_Bd(unsigned int, std::vector<lp::mpq, std::allocator<lp::mpq> >&, lp::indexed_vector<lp::mpq>&);
template lp::lu<lp::mpq, lp::numeric_pair<lp::mpq> >::~lu();
template lp::mpq lp::dot_product<lp::mpq, lp::mpq>(std::vector<lp::mpq, std::allocator<lp::mpq> > const&, std::vector<lp::mpq, std::allocator<lp::mpq> > const&, unsigned int);
template void lp::init_factorization<double, double>(lp::lu<double, double>*&, lp::static_matrix<double, double>&, std::vector<unsigned int, std::allocator<unsigned int> >&, std::vector<int, std::allocator<int> >&, lp::lp_settings&, std::vector<unsigned int, std::allocator<unsigned int> >&);
template void lp::init_factorization<lp::mpq, lp::mpq>(lp::lu<lp::mpq, lp::mpq>*&, lp::static_matrix<lp::mpq, lp::mpq>&, std::vector<unsigned int, std::allocator<unsigned int> >&, std::vector<int, std::allocator<int> >&, lp::lp_settings&, std::vector<unsigned int, std::allocator<unsigned int> >&);
template void lp::init_factorization<lp::mpq, lp::numeric_pair<lp::mpq> >(lp::lu<lp::mpq, lp::numeric_pair<lp::mpq> >*&, lp::static_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >&, std::vector<unsigned int, std::allocator<unsigned int> >&, std::vector<int, std::allocator<int> >&, lp::lp_settings&, std::vector<unsigned int, std::allocator<unsigned int> >&);
#ifdef LEAN_DEBUG
template void lp::print_matrix<double, double>(lp::sparse_matrix<double, double>&, std::ostream & out);
template void lp::print_matrix<double, double>(lp::static_matrix<double, double>&, std::ostream & out);
template bool lp::lu<double, double>::is_correct();
template lp::dense_matrix<double, double> lp::get_B<double, double>(lp::lu<double, double>&);
#endif
template void lp::lu<lp::mpq, lp::mpq>::solve_yB(std::vector<lp::mpq, std::allocator<lp::mpq> >&);
template void lp::lu<double, double>::solve_yB(std::vector<double, std::allocator<double> >&);
template void lp::lu<lp::mpq, lp::mpq>::solve_By(std::vector<lp::mpq, std::allocator<lp::mpq> >&);
template void lp::lu<double, double>::solve_By(std::vector<double, std::allocator<double> >&);
template void lp::lu<lp::mpq, lp::mpq>::replace_column(unsigned int, lp::mpq, lp::indexed_vector<lp::mpq>&);
template void lp::lu<lp::mpq, lp::numeric_pair<lp::mpq> >::replace_column(unsigned int, lp::mpq, lp::indexed_vector<lp::mpq>&);
template void lp::lu<lp::mpq, lp::numeric_pair<lp::mpq> >::solve_yB(std::vector<lp::mpq, std::allocator<lp::mpq> >&);
template void lp::lu<lp::mpq, lp::numeric_pair<lp::mpq> >::solve_By(std::vector<lp::numeric_pair<lp::mpq>, std::allocator<lp::numeric_pair<lp::mpq> > >&);
template void lp::lu<double, double>::init_vector_w(unsigned int, lp::indexed_vector<double>&);
template lp::numeric_pair<lp::mpq> lp::dot_product<lp::mpq, lp::numeric_pair<lp::mpq> >(std::vector<lp::mpq, std::allocator<lp::mpq> > const&, std::vector<lp::numeric_pair<lp::mpq>, std::allocator<lp::numeric_pair<lp::mpq> > > const&, unsigned int);
template bool lp::lu<double, double>::pivot_the_row(int); // NOLINT

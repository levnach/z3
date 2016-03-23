/*
  Copyright (c) 2013 Microsoft Corporation. All rights reserved.
  Released under Apache 2.0 license as described in the file LICENSE.

  Author: Lev Nachmanson
*/
#include <memory>
#include <vector>
#include "util/lp/permutation_matrix.cpp"
#include "util/lp/numeric_pair.h"
template void lp::permutation_matrix<double, double>::apply_from_right(std::vector<double, std::allocator<double> >&);
template void lp::permutation_matrix<double, double>::init(unsigned int);
template bool lp::permutation_matrix<double, double>::is_identity() const;
template void lp::permutation_matrix<double, double>::multiply_by_permutation_from_left(lp::permutation_matrix<double, double>&);
template void lp::permutation_matrix<double, double>::multiply_by_permutation_reverse_from_left(lp::permutation_matrix<double, double>&);
template void lp::permutation_matrix<double, double>::multiply_by_reverse_from_right(lp::permutation_matrix<double, double>&);
template lp::permutation_matrix<double, double>::permutation_matrix(unsigned int, std::vector<unsigned int, std::allocator<unsigned int> > const&);
template void lp::permutation_matrix<double, double>::transpose_from_left(unsigned int, unsigned int);

template void lp::permutation_matrix<lp::mpq, lp::mpq>::apply_from_right(std::vector<lp::mpq, std::allocator<lp::mpq> >&);
template bool lp::permutation_matrix<lp::mpq, lp::mpq>::is_identity() const;
template void lp::permutation_matrix<lp::mpq, lp::mpq>::multiply_by_permutation_from_left(lp::permutation_matrix<lp::mpq, lp::mpq>&);
template void lp::permutation_matrix<lp::mpq, lp::mpq>::multiply_by_permutation_from_right(lp::permutation_matrix<lp::mpq, lp::mpq>&);
template void lp::permutation_matrix<lp::mpq, lp::mpq>::multiply_by_permutation_reverse_from_left(lp::permutation_matrix<lp::mpq, lp::mpq>&);
template void lp::permutation_matrix<lp::mpq, lp::mpq>::multiply_by_reverse_from_right(lp::permutation_matrix<lp::mpq, lp::mpq>&);
template lp::permutation_matrix<lp::mpq, lp::mpq>::permutation_matrix(unsigned int);
template void lp::permutation_matrix<lp::mpq, lp::mpq>::transpose_from_left(unsigned int, unsigned int);
template void lp::permutation_matrix<lp::mpq, lp::mpq>::transpose_from_right(unsigned int, unsigned int);
template void lp::permutation_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >::apply_from_right(std::vector<lp::mpq, std::allocator<lp::mpq> >&);
template bool lp::permutation_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >::is_identity() const;
template void lp::permutation_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >::multiply_by_permutation_from_left(lp::permutation_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >&);
template void lp::permutation_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >::multiply_by_permutation_from_right(lp::permutation_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >&);
template void lp::permutation_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >::multiply_by_permutation_reverse_from_left(lp::permutation_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >&);
template void lp::permutation_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >::multiply_by_reverse_from_right(lp::permutation_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >&);
template lp::permutation_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >::permutation_matrix(unsigned int);
template void lp::permutation_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >::transpose_from_left(unsigned int, unsigned int);
template void lp::permutation_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >::transpose_from_right(unsigned int, unsigned int);
template void lp::permutation_matrix<double, double>::apply_from_left_perm<double>(lp::indexed_vector<double>&, lp::lp_settings&);
template void lp::permutation_matrix<double, double>::apply_from_left_perm<double>(std::vector<double, std::allocator<double> >&);
template void lp::permutation_matrix<double, double>::apply_reverse_from_left<double>(lp::indexed_vector<double>&);
template void lp::permutation_matrix<double, double>::apply_reverse_from_left<double>(std::vector<double, std::allocator<double> >&);
template void lp::permutation_matrix<double, double>::apply_reverse_from_right<double>(std::vector<double, std::allocator<double> >&);
template void lp::permutation_matrix<double, double>::transpose_from_right(unsigned int, unsigned int);
template void lp::permutation_matrix<lp::mpq, lp::mpq>::apply_from_left_perm<lp::mpq>(lp::indexed_vector<lp::mpq>&, lp::lp_settings&);
template void lp::permutation_matrix<lp::mpq, lp::mpq>::apply_from_left_perm<lp::mpq>(std::vector<lp::mpq, std::allocator<lp::mpq> >&);
template void lp::permutation_matrix<lp::mpq, lp::mpq>::apply_reverse_from_left<lp::mpq>(lp::indexed_vector<lp::mpq>&);
template void lp::permutation_matrix<lp::mpq, lp::mpq>::apply_reverse_from_left<lp::mpq>(std::vector<lp::mpq, std::allocator<lp::mpq> >&);
template void lp::permutation_matrix<lp::mpq, lp::mpq>::apply_reverse_from_right<lp::mpq>(std::vector<lp::mpq, std::allocator<lp::mpq> >&);
template void lp::permutation_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >::apply_from_left_perm<lp::mpq>(lp::indexed_vector<lp::mpq>&, lp::lp_settings&);
template void lp::permutation_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >::apply_from_left_perm<lp::numeric_pair<lp::mpq> >(std::vector<lp::numeric_pair<lp::mpq>, std::allocator<lp::numeric_pair<lp::mpq> > >&);
template void lp::permutation_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >::apply_reverse_from_left<lp::mpq>(lp::indexed_vector<lp::mpq>&);
template void lp::permutation_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >::apply_reverse_from_left<lp::mpq>(std::vector<lp::mpq, std::allocator<lp::mpq> >&);
template void lp::permutation_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >::apply_reverse_from_left<lp::numeric_pair<lp::mpq> >(std::vector<lp::numeric_pair<lp::mpq>, std::allocator<lp::numeric_pair<lp::mpq> > >&);
template void lp::permutation_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >::apply_reverse_from_right<lp::mpq>(std::vector<lp::mpq, std::allocator<lp::mpq> >&);
template void lp::permutation_matrix<double, double>::multiply_by_permutation_from_right(lp::permutation_matrix<double, double>&);

#ifdef LEAN_DEBUG
template bool lp::permutation_generator<double, double>::move_next();
template lp::permutation_generator<double, double>::permutation_generator(unsigned int);
#endif
template lp::permutation_matrix<double, double>::permutation_matrix(unsigned int);

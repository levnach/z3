/*
  Copyright (c) 2013 Microsoft Corporation. All rights reserved.
  Released under Apache 2.0 license as described in the file LICENSE.

  Author: Lev Nachmanson
*/
#include <memory>
#include <vector>
#include "util/lp/numeric_pair.h"
#include "util/lp/eta_matrix.cpp"
#ifdef LEAN_DEBUG
template double lp::eta_matrix<double, double>::get_elem(unsigned int, unsigned int) const;
template lp::mpq lp::eta_matrix<lp::mpq, lp::mpq>::get_elem(unsigned int, unsigned int) const;
template lp::mpq lp::eta_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >::get_elem(unsigned int, unsigned int) const;
#endif
template void lp::eta_matrix<double, double>::apply_from_left(std::vector<double, std::allocator<double> >&, lp::lp_settings&);
template void lp::eta_matrix<double, double>::apply_from_right(std::vector<double, std::allocator<double> >&);
template void lp::eta_matrix<double, double>::conjugate_by_permutation(lp::permutation_matrix<double, double>&);
template void lp::eta_matrix<lp::mpq, lp::mpq>::apply_from_left(std::vector<lp::mpq, std::allocator<lp::mpq> >&, lp::lp_settings&);
template void lp::eta_matrix<lp::mpq, lp::mpq>::apply_from_right(std::vector<lp::mpq, std::allocator<lp::mpq> >&);
template void lp::eta_matrix<lp::mpq, lp::mpq>::conjugate_by_permutation(lp::permutation_matrix<lp::mpq, lp::mpq>&);
template void lp::eta_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >::apply_from_left(std::vector<lp::numeric_pair<lp::mpq>, std::allocator<lp::numeric_pair<lp::mpq> > >&, lp::lp_settings&);
template void lp::eta_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >::apply_from_right(std::vector<lp::mpq>&);
template void lp::eta_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >::conjugate_by_permutation(lp::permutation_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >&);
template void lp::eta_matrix<double, double>::apply_from_left_local<double>(lp::indexed_vector<double>&, lp::lp_settings&);
template void lp::eta_matrix<lp::mpq, lp::mpq>::apply_from_left_local<lp::mpq>(lp::indexed_vector<lp::mpq>&, lp::lp_settings&);
template void lp::eta_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >::apply_from_left_local<lp::mpq>(lp::indexed_vector<lp::mpq>&, lp::lp_settings&);

/*
  Copyright (c) 2013 Microsoft Corporation. All rights reserved.
  Released under Apache 2.0 license as described in the file LICENSE.

  Author: Lev Nachmanson
*/
#ifdef LEAN_DEBUG
#include <vector>
#include "util/lp/dense_matrix.cpp"
template lp::dense_matrix<double, double> lp::operator*<double, double>(lp::matrix<double, double>&, lp::matrix<double, double>&);
template void lp::dense_matrix<double, double>::apply_from_left(std::vector<double> &);
template lp::dense_matrix<double, double>::dense_matrix(lp::matrix<double, double> const&);
template lp::dense_matrix<double, double>::dense_matrix(unsigned int, unsigned int);
template lp::dense_matrix<double, double>& lp::dense_matrix<double, double>::operator=(lp::dense_matrix<double, double> const&);
#endif

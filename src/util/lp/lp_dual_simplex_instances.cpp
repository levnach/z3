/*
  Copyright (c) 2013 Microsoft Corporation. All rights reserved.
  Released under Apache 2.0 license as described in the file LICENSE.

  Author: Lev Nachmanson
*/
#include "util/lp/lp_dual_simplex.cpp"
template lp::mpq lp::lp_dual_simplex<lp::mpq, lp::mpq>::get_current_cost() const;
template void lp::lp_dual_simplex<lp::mpq, lp::mpq>::find_maximal_solution();
template double lp::lp_dual_simplex<double, double>::get_current_cost() const;
template void lp::lp_dual_simplex<double, double>::find_maximal_solution();

/*
  Copyright (c) 2013 Microsoft Corporation. All rights reserved.
  Released under Apache 2.0 license as described in the file LICENSE.

  Author: Lev Nachmanson
*/
#include <utility>
#include <memory>
#include <string>
#include <vector>
#include <functional>
#include "util/lp/lp_primal_simplex.cpp"
template bool lp::lp_primal_simplex<double, double>::bounds_hold(std::unordered_map<std::string, double, std::hash<std::string>, std::equal_to<std::string>, std::allocator<std::pair<std::string const, double> > > const&);
template bool lp::lp_primal_simplex<double, double>::row_constraints_hold(std::unordered_map<std::string, double, std::hash<std::string>, std::equal_to<std::string>, std::allocator<std::pair<std::string const, double> > > const&);
template double lp::lp_primal_simplex<double, double>::get_current_cost() const;
template double lp::lp_primal_simplex<double, double>::get_column_value(unsigned int) const;
template lp::lp_primal_simplex<double, double>::~lp_primal_simplex();
template lp::lp_primal_simplex<lp::mpq, lp::mpq>::~lp_primal_simplex();
template lp::mpq lp::lp_primal_simplex<lp::mpq, lp::mpq>::get_current_cost() const;
template lp::mpq lp::lp_primal_simplex<lp::mpq, lp::mpq>::get_column_value(unsigned int) const;
template void lp::lp_primal_simplex<double, double>::find_maximal_solution();
template void lp::lp_primal_simplex<lp::mpq, lp::mpq>::find_maximal_solution();

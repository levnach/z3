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
#include "util/lp/lp_dual_core_solver.cpp"
template lp::lp_dual_core_solver<lp::mpq, lp::mpq>::lp_dual_core_solver(lp::static_matrix<lp::mpq, lp::mpq>&, std::vector<bool, std::allocator<bool> >&, std::vector<lp::mpq, std::allocator<lp::mpq> >&, std::vector<lp::mpq, std::allocator<lp::mpq> >&, std::vector<unsigned int, std::allocator<unsigned int> >&, std::vector<lp::mpq, std::allocator<lp::mpq> >&, std::vector<lp::column_type, std::allocator<lp::column_type> >&, std::vector<lp::mpq, std::allocator<lp::mpq> >&, std::vector<lp::mpq, std::allocator<lp::mpq> >&, lp::lp_settings&, std::unordered_map<unsigned int, std::string, std::hash<unsigned int>, std::equal_to<unsigned int>, std::allocator<std::pair<unsigned int const, std::string> > > const&);
template void lp::lp_dual_core_solver<lp::mpq, lp::mpq>::start_with_initial_basis_and_make_it_dual_feasible();
template void lp::lp_dual_core_solver<lp::mpq, lp::mpq>::solve();
template lp::lp_dual_core_solver<double, double>::lp_dual_core_solver(lp::static_matrix<double, double>&, std::vector<bool, std::allocator<bool> >&, std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&, std::vector<unsigned int, std::allocator<unsigned int> >&, std::vector<double, std::allocator<double> >&, std::vector<lp::column_type, std::allocator<lp::column_type> >&, std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&, lp::lp_settings&, std::unordered_map<unsigned int, std::string, std::hash<unsigned int>, std::equal_to<unsigned int>, std::allocator<std::pair<unsigned int const, std::string> > > const&);
template void lp::lp_dual_core_solver<double, double>::start_with_initial_basis_and_make_it_dual_feasible();
template void lp::lp_dual_core_solver<double, double>::solve();
template void lp::lp_dual_core_solver<lp::mpq, lp::mpq>::restore_non_basis();
template void lp::lp_dual_core_solver<double, double>::restore_non_basis();
template void lp::lp_dual_core_solver<double, double>::revert_to_previous_basis();
template void lp::lp_dual_core_solver<lp::mpq, lp::mpq>::revert_to_previous_basis();

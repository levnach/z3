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
#include "util/lp/lar_core_solver.cpp"
template lp::lar_core_solver<lp::mpq, lp::numeric_pair<lp::mpq> >::lar_core_solver(std::vector<lp::numeric_pair<lp::mpq>, std::allocator<lp::numeric_pair<lp::mpq> > >&, std::vector<lp::column_type, std::allocator<lp::column_type> >&, std::vector<lp::numeric_pair<lp::mpq>, std::allocator<lp::numeric_pair<lp::mpq> > >&, std::vector<lp::numeric_pair<lp::mpq>, std::allocator<lp::numeric_pair<lp::mpq> > >&, std::vector<unsigned int, std::allocator<unsigned int> >&, lp::static_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >&, lp::lp_settings&, std::unordered_map<unsigned int, std::string, std::hash<unsigned int>, std::equal_to<unsigned int>, std::allocator<std::pair<unsigned int const, std::string> > >&, std::vector<lp::numeric_pair<lp::mpq>, std::allocator<lp::numeric_pair<lp::mpq> > >&, std::vector<lp::mpq, std::allocator<lp::mpq> >&);
template void lp::lar_core_solver<lp::mpq, lp::numeric_pair<lp::mpq> >::solve();
template void lp::lar_core_solver<lp::mpq, lp::numeric_pair<lp::mpq> >::prefix();
template void lp::lar_core_solver<lp::mpq, lp::numeric_pair<lp::mpq> >::print_column_info(unsigned int, std::ostream & out);
#ifdef LEAN_DEBUG
template lp::numeric_pair<lp::mpq> lp::lar_core_solver<lp::mpq, lp::numeric_pair<lp::mpq> >::get_deb_inf();
#endif
template bool lp::lar_core_solver<lp::mpq, lp::numeric_pair<lp::mpq> >::is_empty() const;

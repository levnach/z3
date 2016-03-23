/*
  Copyright (c) 2013 Microsoft Corporation. All rights reserved.
  Released under Apache 2.0 license as described in the file LICENSE.

  Author: Lev Nachmanson
*/
#include "util/lp/numeric_pair.h"
#include "util/lp/core_solver_pretty_printer.cpp"
template lp::core_solver_pretty_printer<double, double>::core_solver_pretty_printer(lp::lp_core_solver_base<double, double> &, std::ostream & out);
template void lp::core_solver_pretty_printer<double, double>::print();
template lp::core_solver_pretty_printer<double, double>::~core_solver_pretty_printer();
template lp::core_solver_pretty_printer<lp::mpq, lp::mpq>::core_solver_pretty_printer(lp::lp_core_solver_base<lp::mpq, lp::mpq> &, std::ostream & out);
template void lp::core_solver_pretty_printer<lp::mpq, lp::mpq>::print();
template lp::core_solver_pretty_printer<lp::mpq, lp::mpq>::~core_solver_pretty_printer();
template lp::core_solver_pretty_printer<lp::mpq, lp::numeric_pair<lp::mpq> >::core_solver_pretty_printer(lp::lp_core_solver_base<lp::mpq, lp::numeric_pair<lp::mpq> > &, std::ostream & out);
template lp::core_solver_pretty_printer<lp::mpq, lp::numeric_pair<lp::mpq> >::~core_solver_pretty_printer();
template void lp::core_solver_pretty_printer<lp::mpq, lp::numeric_pair<lp::mpq> >::print();

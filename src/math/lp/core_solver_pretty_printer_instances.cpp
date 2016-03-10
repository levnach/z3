/*
  Copyright (c) 2013 Microsoft Corporation. All rights reserved.
  Released under Apache 2.0 license as described in the file LICENSE.

  Author: Lev Nachmanson
*/
#include "numeric_pair.h"
#include "core_solver_pretty_printer.cpp"
template lean::core_solver_pretty_printer<double, double>::core_solver_pretty_printer(lean::lp_core_solver_base<double, double> &, std::ostream & out);
template void lean::core_solver_pretty_printer<double, double>::print();
template lean::core_solver_pretty_printer<double, double>::~core_solver_pretty_printer();
template lean::core_solver_pretty_printer<mpq, mpq>::core_solver_pretty_printer(lean::lp_core_solver_base<mpq, mpq> &, std::ostream & out);
template void lean::core_solver_pretty_printer<mpq, mpq>::print();
template lean::core_solver_pretty_printer<mpq, mpq>::~core_solver_pretty_printer();
template lean::core_solver_pretty_printer<mpq, lean::numeric_pair<mpq> >::core_solver_pretty_printer(lean::lp_core_solver_base<mpq, lean::numeric_pair<mpq> > &, std::ostream & out);
template lean::core_solver_pretty_printer<mpq, lean::numeric_pair<mpq> >::~core_solver_pretty_printer();
template void lean::core_solver_pretty_printer<mpq, lean::numeric_pair<mpq> >::print();

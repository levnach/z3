/*
  Copyright (c) 2017 Microsoft Corporation
  Author: Nikolaj Bjorner, Lev Nachmanson
*/
#include "util/lp/cut_solver_def.h"
namespace lp {
template <typename T>
T cut_solver<T>::m_local_zero = zero_of_type<T>();
template <> int cut_solver<int>::m_local_zero = 0;
template <> mpq cut_solver<mpq>::m_local_zero = zero_of_type<mpq>();
template lbool cut_solver<mpq>::check();
template void cut_solver<mpq>::print_state(std::ostream&) const;
template void cut_solver<mpq>::pop(unsigned int);
template int cut_solver<mpq>::propagate();
template void cut_solver<mpq>::push();
template cut_solver<int>::cut_solver(std::function<std::string (unsigned)> var_name_function,
                                     std::function<void (unsigned, std::ostream &)> print_constraint_function, lp_settings&);
template cut_solver<mpq>::cut_solver(std::function<std::string (unsigned)> var_name_function,
                                     std::function<void (unsigned, std::ostream &)> print_constraint_function, lp_settings&);
    

template bool cut_solver<mpq>::consistent(const ineq & i) const;
template bool cut_solver<int>::consistent(const ineq & i) const;
template unsigned cut_solver<mpq>::add_ineq(const std::vector<cut_solver<mpq>::monomial> & lhs,
                      const mpq& free_coeff,
                      vector<constraint_index> explanation);
template unsigned cut_solver<int>::add_ineq(const std::vector<cut_solver<int>::monomial> & lhs,
                      const int& free_coeff,
                      vector<constraint_index> explanation);
template void cut_solver<int>::print_literal_bound(std::ostream & o, const cut_solver<int>::literal & t) const;
template void cut_solver<mpq>::print_literal_bound(std::ostream & o, const cut_solver<mpq>::literal & t) const;

}


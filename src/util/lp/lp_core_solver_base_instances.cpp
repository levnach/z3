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
#include "util/lp/lp_core_solver_base.cpp"
template bool lp::lp_core_solver_base<double, double>::A_mult_x_is_off();
template bool lp::lp_core_solver_base<double, double>::basis_heading_is_correct();
template void lp::lp_core_solver_base<double, double>::calculate_pivot_row_of_B_1(unsigned int);
template void lp::lp_core_solver_base<double, double>::calculate_pivot_row_when_pivot_row_of_B1_is_ready();
template bool lp::lp_core_solver_base<double, double>::column_is_dual_feasible(unsigned int) const;
template void lp::lp_core_solver_base<double, double>::fill_reduced_costs_from_m_y_by_rows();
template bool lp::lp_core_solver_base<double, double>::find_x_by_solving();
template lp::non_basic_column_value_position lp::lp_core_solver_base<double, double>::get_non_basic_column_value_position(unsigned int);
template void lp::lp_core_solver_base<double, double>::init_reduced_costs_for_one_iteration();
template lp::lp_core_solver_base<double, double>::lp_core_solver_base(lp::static_matrix<double, double>&, std::vector<double, std::allocator<double> >&, std::vector<unsigned int, std::allocator<unsigned int> >&, std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&, lp::lp_settings&, std::unordered_map<unsigned int, std::string, std::hash<unsigned int>, std::equal_to<unsigned int>, std::allocator<std::pair<unsigned int const, std::string> > > const&, std::vector<lp::column_type, std::allocator<lp::column_type> >&, std::vector<double, std::allocator<double> >&, std::vector<double, std::allocator<double> >&);

template bool lp::lp_core_solver_base<double, double>::print_statistics_with_iterations_and_nonzeroes_and_cost_and_check_that_the_time_is_over(std::string, unsigned int);
template void lp::lp_core_solver_base<double, double>::restore_x(unsigned int, double const&);
template void lp::lp_core_solver_base<double, double>::set_non_basic_x_to_correct_bounds();
template void lp::lp_core_solver_base<double, double>::snap_xN_to_bounds_and_free_columns_to_zeroes();
template void lp::lp_core_solver_base<double, double>::solve_Ax_eq_b();
template void lp::lp_core_solver_base<double, double>::solve_Bd(unsigned int);
template void lp::lp_core_solver_base<double, double>::solve_yB(std::vector<double, std::allocator<double> >&);
template bool lp::lp_core_solver_base<double, double>::update_basis_and_x(int, int, double const&);
template void lp::lp_core_solver_base<double, double>::update_x(unsigned int, double);
template bool lp::lp_core_solver_base<lp::mpq, lp::mpq>::A_mult_x_is_off();
template bool lp::lp_core_solver_base<lp::mpq, lp::mpq>::basis_heading_is_correct();
template void lp::lp_core_solver_base<lp::mpq, lp::mpq>::calculate_pivot_row_of_B_1(unsigned int);
template void lp::lp_core_solver_base<lp::mpq, lp::mpq>::calculate_pivot_row_when_pivot_row_of_B1_is_ready();
template bool lp::lp_core_solver_base<lp::mpq, lp::mpq>::column_is_dual_feasible(unsigned int) const;
template void lp::lp_core_solver_base<lp::mpq, lp::mpq>::fill_reduced_costs_from_m_y_by_rows();
template bool lp::lp_core_solver_base<lp::mpq, lp::mpq>::find_x_by_solving();
template void lp::lp_core_solver_base<lp::mpq, lp::mpq>::init_reduced_costs_for_one_iteration();
template bool lp::lp_core_solver_base<lp::mpq, lp::mpq>::print_statistics_with_iterations_and_nonzeroes_and_cost_and_check_that_the_time_is_over(std::string, unsigned int);
template void lp::lp_core_solver_base<lp::mpq, lp::mpq>::restore_x(unsigned int, lp::mpq const&);
template void lp::lp_core_solver_base<lp::mpq, lp::mpq>::set_non_basic_x_to_correct_bounds();
template void lp::lp_core_solver_base<lp::mpq, lp::mpq>::solve_Ax_eq_b();
template void lp::lp_core_solver_base<lp::mpq, lp::mpq>::solve_Bd(unsigned int);
template void lp::lp_core_solver_base<lp::mpq, lp::mpq>::solve_yB(std::vector<lp::mpq, std::allocator<lp::mpq> >&);
template bool lp::lp_core_solver_base<lp::mpq, lp::mpq>::update_basis_and_x(int, int, lp::mpq const&);
template void lp::lp_core_solver_base<lp::mpq, lp::mpq>::update_x(unsigned int, lp::mpq);
template void lp::lp_core_solver_base<lp::mpq, lp::numeric_pair<lp::mpq> >::calculate_pivot_row_of_B_1(unsigned int);
template void lp::lp_core_solver_base<lp::mpq, lp::numeric_pair<lp::mpq> >::calculate_pivot_row_when_pivot_row_of_B1_is_ready();
template void lp::lp_core_solver_base<lp::mpq, lp::numeric_pair<lp::mpq> >::init();
template void lp::lp_core_solver_base<lp::mpq, lp::numeric_pair<lp::mpq> >::init_basis_heading();
template void lp::lp_core_solver_base<lp::mpq, lp::numeric_pair<lp::mpq> >::init_reduced_costs_for_one_iteration();
template lp::lp_core_solver_base<lp::mpq, lp::numeric_pair<lp::mpq> >::lp_core_solver_base(lp::static_matrix<lp::mpq, lp::numeric_pair<lp::mpq> >&, std::vector<lp::numeric_pair<lp::mpq>, std::allocator<lp::numeric_pair<lp::mpq> > >&, std::vector<unsigned int, std::allocator<unsigned int> >&, std::vector<lp::numeric_pair<lp::mpq>, std::allocator<lp::numeric_pair<lp::mpq> > >&, std::vector<lp::mpq, std::allocator<lp::mpq> >&, lp::lp_settings&, std::unordered_map<unsigned int, std::string, std::hash<unsigned int>, std::equal_to<unsigned int>, std::allocator<std::pair<unsigned int const, std::string> > > const&, std::vector<lp::column_type, std::allocator<lp::column_type> >&, std::vector<lp::numeric_pair<lp::mpq>, std::allocator<lp::numeric_pair<lp::mpq> > >&, std::vector<lp::numeric_pair<lp::mpq>, std::allocator<lp::numeric_pair<lp::mpq> > >&);
template bool lp::lp_core_solver_base<lp::mpq, lp::numeric_pair<lp::mpq> >::print_statistics_with_cost_and_check_that_the_time_is_over(unsigned int, lp::numeric_pair<lp::mpq>);
template void lp::lp_core_solver_base<lp::mpq, lp::numeric_pair<lp::mpq> >::snap_xN_to_bounds();
template void lp::lp_core_solver_base<lp::mpq, lp::numeric_pair<lp::mpq> >::solve_Bd(unsigned int);
template bool lp::lp_core_solver_base<lp::mpq, lp::numeric_pair<lp::mpq> >::update_basis_and_x(int, int, lp::numeric_pair<lp::mpq> const&);
template void lp::lp_core_solver_base<lp::mpq, lp::numeric_pair<lp::mpq> >::update_x(unsigned int, lp::numeric_pair<lp::mpq>);
template lp::lp_core_solver_base<lp::mpq, lp::mpq>::lp_core_solver_base(lp::static_matrix<lp::mpq, lp::mpq>&, std::vector<lp::mpq, std::allocator<lp::mpq> >&, std::vector<unsigned int, std::allocator<unsigned int> >&, std::vector<lp::mpq, std::allocator<lp::mpq> >&, std::vector<lp::mpq, std::allocator<lp::mpq> >&, lp::lp_settings&, std::unordered_map<unsigned int, std::string, std::hash<unsigned int>, std::equal_to<unsigned int>, std::allocator<std::pair<unsigned int const, std::string> > > const&, std::vector<lp::column_type, std::allocator<lp::column_type> >&, std::vector<lp::mpq, std::allocator<lp::mpq> >&, std::vector<lp::mpq, std::allocator<lp::mpq> >&);
template bool lp::lp_core_solver_base<lp::mpq, lp::numeric_pair<lp::mpq> >::print_statistics_with_iterations_and_check_that_the_time_is_over(unsigned int);
template std::string lp::lp_core_solver_base<double, double>::column_name(unsigned int) const;
template void lp::lp_core_solver_base<double, double>::pretty_print(std::ostream & out);
template void lp::lp_core_solver_base<double, double>::restore_state(double*, double*);
template void lp::lp_core_solver_base<double, double>::save_state(double*, double*);
template std::string lp::lp_core_solver_base<lp::mpq, lp::mpq>::column_name(unsigned int) const;
template void lp::lp_core_solver_base<lp::mpq, lp::mpq>::pretty_print(std::ostream & out);
template void lp::lp_core_solver_base<lp::mpq, lp::mpq>::restore_state(lp::mpq*, lp::mpq*);
template void lp::lp_core_solver_base<lp::mpq, lp::mpq>::save_state(lp::mpq*, lp::mpq*);
template std::string lp::lp_core_solver_base<lp::mpq, lp::numeric_pair<lp::mpq> >::column_name(unsigned int) const;
template void lp::lp_core_solver_base<lp::mpq, lp::numeric_pair<lp::mpq> >::pretty_print(std::ostream & out);
template void lp::lp_core_solver_base<lp::mpq, lp::numeric_pair<lp::mpq> >::restore_state(lp::mpq*, lp::mpq*);
template void lp::lp_core_solver_base<lp::mpq, lp::numeric_pair<lp::mpq> >::save_state(lp::mpq*, lp::mpq*);
template void lp::lp_core_solver_base<lp::mpq, lp::numeric_pair<lp::mpq> >::solve_yB(std::vector<lp::mpq, std::allocator<lp::mpq> >&);
template void lp::lp_core_solver_base<double, double>::init_lu();
template void lp::lp_core_solver_base<lp::mpq, lp::mpq>::init_lu();
template int lp::lp_core_solver_base<double, double>::pivots_in_column_and_row_are_different(int, int) const;
template int lp::lp_core_solver_base<lp::mpq, lp::mpq>::pivots_in_column_and_row_are_different(int, int) const;

/*
  Copyright (c) 2017 Microsoft Corporation
  Author: Nikolaj Bjorner, Lev Nachmanson
*/
#pragma once
#include "util/vector.h"
#include <utility>
#include "util/debug.h"
#include "util/buffer.h"
#include <unordered_map>
#include <unordered_set>
#include <string>
#include "util/lp/lar_constraints.h"
#include <functional>
#include "util/lp/lar_core_solver.h"
#include <algorithm>
#include "util/lp/numeric_pair.h"
#include "util/lp/scaler.h"
#include "util/lp/lp_primal_core_solver.h"
#include "util/lp/random_updater.h"
#include <stack>
#include "util/lp/stacked_value.h"
#include "util/lp/stacked_vector.h"
#include "util/lp/stacked_unordered_set.h"
#include "util/lp/iterator_on_pivot_row.h"
#include "util/lp/implied_bound.h"
#include "util/lp/bound_analyzer_on_row.h"
#include "util/lp/iterator_on_term_with_basis_var.h"
#include "util/lp/iterator_on_row.h"
#include "util/lp/quick_xplain.h"
#include "util/lp/conversion_helper.h"
#include "util/lp/int_solver.h"
#include "util/lp/nra_solver.h"

namespace lp {


class lar_solver : public column_namer {
    class ext_var_info {
        unsigned m_ext_j; // the external index
        bool m_is_integer;
    public:
        ext_var_info(unsigned j): ext_var_info(j, false) {}
        ext_var_info(unsigned j , bool is_int) : m_ext_j(j), m_is_integer(is_int) {}
        unsigned ext_j() const { return m_ext_j;}
        bool is_integer() const {return m_is_integer;}
    };
    //////////////////// fields //////////////////////////
    lp_settings m_settings;
    lp_status m_status;
    stacked_value<simplex_strategy_enum> m_simplex_strategy;
    std::unordered_map<unsigned, ext_var_info> m_ext_vars_to_columns;
    vector<unsigned> m_columns_to_ext_vars_or_term_indices;
    stacked_vector<ul_pair> m_columns_to_ul_pairs;
    vector<lar_base_constraint*> m_constraints;
public :
    const vector<lar_base_constraint*>& constraints() const {
        return m_constraints;
    }
private:
    stacked_value<unsigned> m_constraint_count;
    // the set of column indices j such that bounds have changed for j
    int_set m_columns_with_changed_bound;
    int_set m_rows_with_changed_bounds;
    int_set m_basic_columns_with_changed_cost;
    stacked_value<int> m_infeasible_column_index; // such can be found at the initialization step
    stacked_value<unsigned> m_term_count;
    vector<lar_term*> m_terms;
    const var_index m_terms_start_index;
    indexed_vector<mpq> m_column_buffer;
public:
    lar_core_solver m_mpq_lar_core_solver;
private:
    std::function<void (unsigned)> m_tracker_of_x_change;
    int_solver * m_int_solver;
public:
    void set_int_solver(int_solver * int_slv) {
        m_int_solver = int_slv;
    }
    int_solver * get_int_solver() {
        return m_int_solver;
    }
    unsigned constraint_count() const;
    const lar_base_constraint& get_constraint(unsigned ci) const;
    int_set m_inf_int_set; 
    ////////////////// methods ////////////////////////////////
    static_matrix<mpq, numeric_pair<mpq>> & A_r();
    static_matrix<mpq, numeric_pair<mpq>> const & A_r() const;
    static_matrix<double, double> & A_d();
    static_matrix<double, double > const & A_d() const;
    
    static bool valid_index(unsigned j){ return static_cast<int>(j) >= 0;}

    bool column_is_int(unsigned j) const;
    bool column_value_is_int(unsigned j) const {
        return m_mpq_lar_core_solver.m_r_x[j].is_int();
    }

    const impq& get_column_value(unsigned j) const {
        return m_mpq_lar_core_solver.m_r_x[j];
    }
    bool is_term(var_index j) const;
    bool column_is_fixed(unsigned j) const;
public:

    // init region
    bool strategy_is_undecided() const;

    var_index add_var(unsigned ext_j, bool is_integer);

    void register_new_ext_var_index(unsigned ext_v, bool is_int);

    bool term_is_int(const lar_term * t) const;

    bool var_is_int(var_index v) const;

    bool ext_var_is_int(var_index ext_var) const;
    
    void add_non_basic_var_to_core_fields(unsigned ext_j, bool is_int);

    void add_new_var_to_core_fields_for_doubles(bool register_in_basis);

    void add_new_var_to_core_fields_for_mpq(bool register_in_basis);


    var_index add_term_undecided(const vector<std::pair<mpq, var_index>> & coeffs,
                                 const mpq &m_v);

    // terms
    var_index add_term(const vector<std::pair<mpq, var_index>> & coeffs,
                       const mpq &m_v);

    void add_row_for_term(const lar_term * term, unsigned term_ext_index);

    void add_row_from_term_no_constraint(const lar_term * term, unsigned term_ext_index);

    void add_basic_var_to_core_fields();

    constraint_index add_var_bound(var_index j, lconstraint_kind kind, const mpq & right_side) ;

    void update_column_type_and_bound(var_index j, lconstraint_kind kind, const mpq & right_side, constraint_index constr_index);

    void add_var_bound_on_constraint_for_term(var_index j, lconstraint_kind kind, const mpq & right_side, constraint_index ci);


    void add_constraint_from_term_and_create_new_column_row(unsigned term_j, const lar_term* term,
                                                            lconstraint_kind kind, const mpq & right_side);

    void decide_on_strategy_and_adjust_initial_state();

    void adjust_initial_state();

    void adjust_initial_state_for_lu();

    void adjust_initial_state_for_tableau_rows();

    // this fills the last row of A_d and sets the basis column: -1 in the last column of the row
    void fill_last_row_of_A_d(static_matrix<double, double> & A, const lar_term* ls);

    void update_free_column_type_and_bound(var_index j, lconstraint_kind kind, const mpq & right_side, constraint_index constr_ind);

    void update_upper_bound_column_type_and_bound(var_index j, lconstraint_kind kind, const mpq & right_side, constraint_index ci);
    
    void update_boxed_column_type_and_bound(var_index j, lconstraint_kind kind, const mpq & right_side, constraint_index ci);
    void update_lower_bound_column_type_and_bound(var_index j, lconstraint_kind kind, const mpq & right_side, constraint_index ci);

    void update_fixed_column_type_and_bound(var_index j, lconstraint_kind kind, const mpq & right_side, constraint_index ci);
    //end of init region
    lp_settings & settings();

    lp_settings const & settings() const;

    void clear();
    lar_solver();
    void set_track_pivoted_rows(bool v);

    bool get_track_pivoted_rows() const;
    
    virtual ~lar_solver();

    unsigned adjust_term_index(unsigned j) const;


    bool use_lu() const;
    
    bool sizes_are_correct() const;
 
    void print_implied_bound(const implied_bound& be, std::ostream & out) const;
    
    bool implied_bound_is_correctly_explained(implied_bound const & be, const vector<std::pair<mpq, unsigned>> & explanation) const;
    
    void analyze_new_bounds_on_row(
        unsigned row_index,
        bound_propagator & bp);

    void analyze_new_bounds_on_row_tableau(
        unsigned row_index,
        bound_propagator & bp
                                           );

    
    void substitute_basis_var_in_terms_for_row(unsigned i);
    
    void calculate_implied_bounds_for_row(unsigned i, bound_propagator & bp);

  
    linear_combination_iterator<mpq> * create_new_iter_from_term(unsigned term_index) const;

    unsigned adjust_column_index_to_term_index(unsigned j) const;
    
    void propagate_bounds_on_a_term(const lar_term& t, bound_propagator & bp, unsigned term_offset);


    void explain_implied_bound(implied_bound & ib, bound_propagator & bp);


    bool term_is_used_as_row(unsigned term) const;
    
    void propagate_bounds_on_terms(bound_propagator & bp);


    // goes over touched rows and tries to induce bounds
    void propagate_bounds_for_touched_rows(bound_propagator & bp);

    lp_status get_status() const;

    void set_status(lp_status s);

    lp_status find_feasible_solution();   
    
    lp_status solve();

    void fill_explanation_from_infeasible_column(explanation_t & evidence) const;

    
    unsigned get_total_iterations() const;
    // see http://research.microsoft.com/projects/z3/smt07.pdf
    // This method searches for a feasible solution with as many different values of variables, reverenced in vars, as it can find
    // Attention, after a call to this method the non-basic variables don't necesserarly stick to their bounds anymore
    vector<unsigned> get_list_of_all_var_indices() const;
    void push();

    static void clean_popped_elements(unsigned n, int_set& set);

    static void shrink_inf_set_after_pop(unsigned n, int_set & set);

    
    void pop(unsigned k);
    
    vector<constraint_index> get_all_constraint_indices() const;

    bool maximize_term_on_tableau(const vector<std::pair<mpq, var_index>> & term,
                                  impq &term_max);

    bool costs_are_zeros_for_r_solver() const;
    bool reduced_costs_are_zeroes_for_r_solver() const;
    
    void set_costs_to_zero(const vector<std::pair<mpq, var_index>> & term);

    void prepare_costs_for_r_solver(const vector<std::pair<mpq, var_index>> & term);
    
    bool maximize_term_on_corrected_r_solver(const vector<std::pair<mpq, var_index>> & term,
                                             impq &term_max);    
    // starting from a given feasible state look for the maximum of the term
    // return true if found and false if unbounded
    bool maximize_term(const vector<std::pair<mpq, var_index>> & term,
                       impq &term_max);
    

    
    const lar_term &  get_term(unsigned j) const;

    void pop_core_solver_params();

    void pop_core_solver_params(unsigned k);


    void set_upper_bound_witness(var_index j, constraint_index ci);

    void set_lower_bound_witness(var_index j, constraint_index ci);


    void substitute_terms_in_linear_expression( const vector<std::pair<mpq, var_index>>& left_side_with_terms,
                                                vector<std::pair<mpq, var_index>> &left_side, mpq & free_coeff) const;


    void detect_rows_of_bound_change_column_for_nbasic_column(unsigned j);


    
    void detect_rows_of_bound_change_column_for_nbasic_column_tableau(unsigned j);

    bool use_tableau() const;

    bool use_tableau_costs() const;
    
    void detect_rows_of_column_with_bound_change(unsigned j);

    void adjust_x_of_column(unsigned j);

    bool row_is_correct(unsigned i) const;
    
    bool ax_is_correct() const;

    bool tableau_with_costs() const;

    bool costs_are_used() const;
    
    void change_basic_columns_dependend_on_a_given_nb_column(unsigned j, const numeric_pair<mpq> & delta);

    void update_x_and_inf_costs_for_column_with_changed_bounds(unsigned j);

    
    void detect_rows_with_changed_bounds_for_column(unsigned j);
    
    void detect_rows_with_changed_bounds();

    void update_x_and_inf_costs_for_columns_with_changed_bounds();

    void update_x_and_inf_costs_for_columns_with_changed_bounds_tableau();

    
    void solve_with_core_solver();

    
    numeric_pair<mpq> get_basic_var_value_from_row_directly(unsigned i);
    
    
    
    numeric_pair<mpq> get_basic_var_value_from_row(unsigned i);

    template <typename K, typename L>
    void add_last_rows_to_lu(lp_primal_core_solver<K,L> & s);
    
    bool x_is_correct() const;

    bool var_is_registered(var_index vj) const;

    unsigned constraint_stack_size() const;

    void fill_last_row_of_A_r(static_matrix<mpq, numeric_pair<mpq>> & A, const lar_term * ls);

    template <typename U, typename V>
    void create_matrix_A(static_matrix<U, V> & matr);

    template <typename U, typename V>
    void copy_from_mpq_matrix(static_matrix<U, V> & matr);


    bool try_to_set_fixed(column_info<mpq> & ci);

    column_type get_column_type(const column_info<mpq> & ci);

    std::string get_column_name(unsigned j) const;

    bool all_constrained_variables_are_registered(const vector<std::pair<mpq, var_index>>& left_side);

    constraint_index add_constraint(const vector<std::pair<mpq, var_index>>& left_side_with_terms, lconstraint_kind kind_par, const mpq& right_side_parm);
    bool all_constraints_hold() const;
    bool constraint_holds(const lar_base_constraint & constr, std::unordered_map<var_index, mpq> & var_map) const;
    bool the_relations_are_of_same_type(const vector<std::pair<mpq, unsigned>> & evidence, lconstraint_kind & the_kind_of_sum) const;

    static void register_in_map(std::unordered_map<var_index, mpq> & coeffs, const lar_base_constraint & cn, const mpq & a);
    static void register_monoid_in_map(std::unordered_map<var_index, mpq> & coeffs, const mpq & a, unsigned j);


    bool the_left_sides_sum_to_zero(const vector<std::pair<mpq, unsigned>> & evidence) const;

    bool the_right_sides_do_not_sum_to_zero(const vector<std::pair<mpq, unsigned>> & evidence);

    bool explanation_is_correct(const vector<std::pair<mpq, unsigned>>& explanation) const;

    bool inf_explanation_is_correct() const;

    mpq sum_of_right_sides_of_explanation(const vector<std::pair<mpq, unsigned>> & explanation) const;

    bool has_lower_bound(var_index var, constraint_index& ci, mpq& value, bool& is_strict);
    
    bool has_upper_bound(var_index var, constraint_index& ci, mpq& value, bool& is_strict);


    void get_infeasibility_explanation(vector<std::pair<mpq, constraint_index>> & explanation) const;

    void get_infeasibility_explanation_for_inf_sign(
        vector<std::pair<mpq, constraint_index>> & explanation,
        const vector<std::pair<mpq, unsigned>> & inf_row,
        int inf_sign) const;



    void get_model(std::unordered_map<var_index, mpq> & variable_values) const;

    void get_model_do_not_care_about_diff_vars(std::unordered_map<var_index, mpq> & variable_values) const;

    std::string get_variable_name(var_index vi) const;

    // ********** print region start
    void print_constraint(constraint_index ci, std::ostream & out) const;

    void print_constraints(std::ostream& out) const ;

    void print_terms(std::ostream& out) const;

    void print_left_side_of_constraint(const lar_base_constraint * c, std::ostream & out) const;

    void print_term(lar_term const& term, std::ostream & out) const;

    void print_term_as_indices(lar_term const& term, std::ostream & out) const;

    mpq get_left_side_val(const lar_base_constraint &  cns, const std::unordered_map<var_index, mpq> & var_map) const;

    void print_constraint(const lar_base_constraint * c, std::ostream & out) const;

    void fill_var_set_for_random_update(unsigned sz, var_index const * vars, vector<unsigned>& column_list);

    void random_update(unsigned sz, var_index const * vars);
    void pivot_fixed_vars_from_basis();
    void pop();
    bool column_represents_row_in_tableau(unsigned j);
    void make_sure_that_the_bottom_right_elem_not_zero_in_tableau(unsigned i, unsigned j);
    void remove_last_row_and_column_from_tableau(unsigned j);
    void remove_last_column_from_A();

    void remove_last_column_from_basis_tableau(unsigned j);
    void remove_last_column_from_tableau();
    void pop_tableau();
    void clean_inf_set_of_r_solver_after_pop();
    void shrink_explanation_to_minimum(vector<std::pair<mpq, constraint_index>> & explanation) const;

    
    
    bool column_value_is_integer(unsigned j) const {
        return get_column_value(j).is_int();
    }

    bool column_is_real(unsigned j) const {
        return !column_is_int(j);
    }	
	
    bool model_is_int_feasible() const;

    const impq & column_lower_bound(unsigned j) const {
        return m_mpq_lar_core_solver.lower_bound(j);
    }

    const impq & column_upper_bound(unsigned j) const {
        return m_mpq_lar_core_solver.upper_bound(j);
    }

    bool column_is_bounded(unsigned j) const {
        return m_mpq_lar_core_solver.column_is_bounded(j);
    }

    void get_bound_constraint_witnesses_for_column(unsigned j, constraint_index & lc, constraint_index & uc) const {
        const ul_pair & ul = m_columns_to_ul_pairs[j];
        lc = ul.lower_bound_witness();
        uc = ul.upper_bound_witness();
    }
    indexed_vector<mpq> & get_column_in_lu_mode(unsigned j) {
        m_column_buffer.clear();
        m_column_buffer.resize(A_r().row_count());
        m_mpq_lar_core_solver.m_r_solver.solve_Bd(j, m_column_buffer);
        return m_column_buffer;
    }
    
    bool bound_is_integer_for_integer_column(unsigned j, const mpq & right_side) const;
    linear_combination_iterator<mpq> * get_iterator_on_row(unsigned i) {
        return m_mpq_lar_core_solver.m_r_solver.get_iterator_on_row(i);
    }

    unsigned get_base_column_in_row(unsigned row_index) const {
        return m_mpq_lar_core_solver.m_r_solver.get_base_column_in_row(row_index);
    }

    constraint_index get_column_upper_bound_witness(unsigned j) const {
        return m_columns_to_ul_pairs()[j].upper_bound_witness();
    }

    constraint_index get_column_lower_bound_witness(unsigned j) const {
        return m_columns_to_ul_pairs()[j].lower_bound_witness();
    }

    void subs_term_columns(lar_term& t) {
        vector<std::pair<mpq, unsigned> > pol;
        for (const auto & m : t.m_coeffs) {
            pol.push_back(std::make_pair(m.second, adjust_column_index_to_term_index(m.first)));
        }
        mpq v = t.m_v;
        vector<std::pair<mpq, unsigned>> pol_after_subs;
        // todo : remove the call to substitute_terms_in_linear_expression, when theory_lra handles the terms indices
        substitute_terms_in_linear_expression(pol, pol_after_subs, v);
        t.clear();
        t = lar_term(pol_after_subs, v);
    }

    bool inf_int_set_is_correct_for_column(unsigned j) const {
        if (m_inf_int_set.contains(j) != (column_is_int(j) && (!column_value_is_integer(j)))) {
            TRACE("arith_int",
                  tout << "j= " << j <<
                  " inf_int_set().contains(j) = " << m_inf_int_set.contains(j) <<
                  ", column_is_int(j) = "   << column_is_int(j) <<
                  "\n column_value_is_integer(j) = " << column_value_is_integer(j) <<
                  ", val = " << get_column_value(j) << std::endl;); 
            return false;
        }
        return true;
    }
    
    bool inf_int_set_is_correct() const {
        if (!has_int_var())
            return true;
        for (unsigned j = 0; j < A_r().column_count(); j++) {
            if (inf_int_set_is_correct_for_column(j) == false)
                return false;
        }
        return true;
    }
    bool has_int_var() const;
    void call_assignment_tracker(unsigned j) {
        if (!var_is_int(j)) {
            lp_assert(m_inf_int_set.contains(j) == false);
            return;
        }
        if (m_mpq_lar_core_solver.m_r_x[j].is_int())
            m_inf_int_set.erase(j);
        else
            m_inf_int_set.insert(j);
    }

    lar_core_solver & get_core_solver() { return m_mpq_lar_core_solver; }
    bool column_corresponds_to_term(unsigned) const;
    void catch_up_in_updating_int_solver();
};
}

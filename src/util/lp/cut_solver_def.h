/*
  Copyright (c) 2017 Microsoft Corporation
  Author: Nikolaj Bjorner, Lev Nachmanson
*/
#include "util/lp/cut_solver.h"
namespace lp {

template <typename T>
lbool cut_solver<T>::check() {
    init_search();
    while (true) {
        TRACE("cs_ch", tout << "inside loop\n";);
        lbool r = bounded_search();
        if (r != lbool::l_undef)
            return r;
        restart();
        simplify_problem();
        if (check_inconsistent()) return lbool::l_false;
        gc();
    }
}

template <typename T>
void cut_solver<T>::init_search() {
    lp_assert(m_explanation.size() == 0);
    lp_assert(m_changed_vars.size() == 0);
}


// returns true if there is no conflict and false otherwise
template <typename T>
bool cut_solver<T>::propagate_inequality(unsigned i) {
    TRACE("ba_int", trace_print_ineq(tout, i););
    const ineq & in = m_ineqs[i];
    if (in.is_simple())
        return true;
    // consider a special case for inequalities with just two variables
    unsigned the_only_unlim;
    int r = lower_analize(in, the_only_unlim);
    if (r == 0) {
        T b;
        lower(in.m_poly, b);
        if (is_pos(b)) {
            TRACE("cs_inconsistent", trace_print_ineq(tout, i););
            return false;
        } else {
            propagate_inequality_on_lower(i, b);
        }
    } else if (r == 1) {
        propagate_inequality_only_one_unlim(i, the_only_unlim);
    }
    return true;
}
template <typename T>
void cut_solver<T>::print_trail(std::ostream & out) const {
    for (const auto & l : m_trail) {
        print_literal(out, l);
        out << "\n";
    }
}
template <typename T>
void cut_solver<T>::print_state(std::ostream & out) const {
    out << "ineqs:\n";
    for (const auto & i: m_ineqs) {
        print_ineq(out, i);
    }
    out << "end of ineqs\n";
    out << "trail\n";
    print_trail(out);
    out << "end of trail\n";
    out << "var_infos\n";
    for (const auto & v: m_var_infos) {
        print_var_info(out, v);
    }
    out << "end of var_infos\n";
    out << "end of state dump" << std::endl;
}
template <typename T>
lbool cut_solver<T>::bounded_search() {
    lbool is_sat = propagate_and_backjump_step();
    if (is_sat != lbool::l_undef)
        return is_sat; //start here

    gc();

    if (!decide()) {
        lbool is_sat = final_check();
        if (is_sat != lbool::l_undef) {
            return is_sat;
        }
    }
    return lbool::l_undef;
}

template <typename T>
void cut_solver<T>::print_literal_bound(std::ostream & o, const literal & t) const {
    o << "BOUND: ";
    o << get_column_name(t.m_var_index) << " ";
    if (t.m_is_lower)
        o << ">= ";
    else
        o << "<= ";
    o << t.m_bound;
    if (t.m_tight_explanation_ineq_index >= 0) {
        o << "  tight ineq ";
        print_ineq(o, t.m_tight_explanation_ineq_index);
    }
    if (t.m_ineq_index >= 0) {
        if (t.m_ineq_index == t.m_tight_explanation_ineq_index) {
            o << " ineq is the same as the tight one\n";
        } else {
            o << "  ineq ";
            print_ineq(o, t.m_ineq_index);
        }
    }
}


template <typename T>
bool cut_solver<T>::decide() {
    int j = find_non_fixed_var();
    if (j < 0)
        return false;
    TRACE("cs_dec", tout << "decided " << var_name(j) << " var index = " << j << "\n";);
    decide_var_on_bound(j, flip_coin());
    return true;
}

template <typename T>
lbool cut_solver<T>::propagate_and_backjump_step() {
    int incostistent_ineq = propagate();
    TRACE("cs_dec", tout << "trail = \n"; print_trail(tout); tout << "end of trail\n";);

    if (incostistent_ineq >= 0) {
        if (at_base_lvl()) {
            fill_conflict_explanation(incostistent_ineq, m_trail.size());
            return lbool::l_false;
        }
        resolve_conflict(incostistent_ineq);
    }

    if(!all_vars_are_fixed())
        return lbool::l_undef;

    return lbool::l_true;
}

// returns -1 if consistent, or the index of an incostistent inequality
template <typename T>
int cut_solver<T>::inconsistent() const {
    for (int i=0;i<m_ineqs.size();i++)
        if (lower(i) > 0)
            return i;

    return -1;
}
template <typename T>
bool cut_solver<T>::decision_is_redundant_for_ineq(const ineq& i, const literal & l) const {
    const T & coeff = i.coeff(l.m_var_index);
    if (is_zero(coeff))
        return true;
    if (is_pos(coeff)) {
        return l.m_is_lower;
    }
    return !l.m_is_lower;
}
// i is an incostistent inequality
template <typename T>
bool cut_solver<T>::resolve_conflict_on_trail_bound(const ineq & i, const literal & l) {
    lp_assert(lower_val(i) > 0);
    lp_assert(l.m_tag == literal_type::BOUND);
    if (l.is_decided()) {
        if (decision_is_redundant_for_ineq(i, l))
            pop(); // skip decision
        else
            lp_assert(false); // not implemented
        
    } else {
        lp_assert(false);// not implemented
    }
        
    return true;
}

template <typename T>
bool cut_solver<T>::resolve_conflict(int inconstistent_ineq) {
    lp_assert(!at_base_lvl());
    ineq & i = m_ineqs[inconstistent_ineq];
    const auto & l = m_trail.back();
    TRACE("cut_solver_state_inconsistent", print_literal(tout, l););
    if (l.m_tag == literal_type::BOUND)
        return resolve_conflict_on_trail_bound(i, l);
    else 
        lp_assert(false);
    
    return true; 
    /*
      while (true) {
      bool r = resolve_conflict_core();
      // after pop, clauses are reinitialized, 
      // this may trigger another conflict.
      if (!r)
      return false;
      if (!inconsistent())
      return true;
      }*/
}

template <typename T>
void cut_solver<T>::print_scope(std::ostream& out) const {
    out << "trail_size = " << m_scope().m_trail_size << ", ineqs_size = " << m_scope().m_ineqs_size <<
        ", div_constraint size = " << m_scope().m_div_constraints_size << "\n";
}

template <typename T>
void cut_solver<T>::pop(unsigned k) {
    TRACE("trace_push_pop_in_cut_solver", tout << "before pop\n";print_state(tout););
    m_scope.pop(k);
    TRACE("trace_push_pop_in_cut_solver", tout << "scope = ";print_scope(tout); tout << "\n";);
    m_trail.resize(m_scope().m_trail_size);
    pop_ineqs();
    pop_div_constraints();
    pop_var_domains(k);
    m_decision_has_been_made.pop(k);
    TRACE("trace_push_pop_in_cut_solver", tout << "after pop\n";print_state(tout););
}

template <typename T>
void cut_solver<T>::push() {
    TRACE("trace_push_pop_in_cut_solver", print_state(tout););
    m_scope = scope(m_trail.size(), m_ineqs.size(), m_div_constraints.size());
    m_scope.push();
    push_var_domains();
    m_decision_has_been_made.push();
}

template <typename T>
cut_solver<T>::cut_solver(std::function<std::string (unsigned)> var_name_function,
                          std::function<void (unsigned, std::ostream &)> print_constraint_function,
                          lp_settings & settings
                          ) : m_var_name_function(var_name_function),
                              m_print_constraint_function(print_constraint_function),
                              m_decision_has_been_made(false),
                              m_settings(settings)
{}

// returns -1 if there is no conflict and the index of the conflict inequality otherwise
template <typename T>
int cut_solver<T>::propagate() {
    propagate_simple_ineqs();
    return propagate_ineqs_for_changed_vars();
}

template <typename T>
void cut_solver<T>::decide_var_on_bound(unsigned j, bool decide_on_lower) {
    push();    
    vector<constraint_index> explanation; // empty in this case
    T b;
    std::vector<monomial> lhs;
    if (decide_on_lower) {
        m_var_infos[j].m_domain.get_lower_bound(b);
        lhs.push_back(monomial(one_of_type<T>(), j));
        b = -b;
    }
    else {
        m_var_infos[j].m_domain.get_upper_bound(b);
        lhs.push_back(monomial(-one_of_type<T>(), j));
    }
    m_decision_has_been_made = true;
    add_ineq_in_local_vars(lhs, b, explanation);
}


template <typename T>
void cut_solver<T>::decide_var() {
    int j = find_non_fixed_var();
    lp_assert(j >= 0);
    decide_var_on_bound(j, flip_coin());
}

template <typename T>
bool cut_solver<T>::propagate_simple_ineq(unsigned ineq_index) {
    const ineq & t = m_ineqs[ineq_index];
    TRACE("cut_solver_state_simpe_ineq",   print_ineq(tout, t); tout << std::endl;);
    var_index j = t.m_poly.m_coeffs[0].var();
        
    bound_result br = bound(t, j);
    TRACE("cut_solver_state_simpe_ineq", tout << "bound result = {"; br.print(tout); tout << "}\n";
          tout << "domain of " << get_column_name(j) << " = "; m_var_infos[j].m_domain.print(tout);
          tout << "\n";
          );
        
    if (improves(j, br)) {
        literal l(j, br.m_type == bound_type::LOWER, br.m_bound, ineq_index);
        l.m_tight_explanation_ineq_index = ineq_index;
        push_literal_to_trail(l);
        restrict_var_domain_with_bound_result(j, br);
        TRACE("cut_solver_state_simpe_ineq", tout <<"improved domain = ";
              m_var_infos[j].m_domain.print(tout);
              tout<<"\n";
              tout << "literal = "; print_literal(tout, l);
              tout <<"\n";
              );
        if (j >= m_changed_vars.data_size())
            m_changed_vars.resize(j + 1);
        m_changed_vars.insert(j);
        return true;
    }
        
    TRACE("cut_solver_state", tout <<"no improvement\n";);
    return false;
}
    

template <typename T>
bool cut_solver<T>::propagate_simple_ineqs() {
    bool ret = false;
    for (unsigned i = 0; i < m_ineqs.size(); i++) {
        if (m_ineqs[i].is_simple() && propagate_simple_ineq(i)){
            ret = true;
        }
    }
    return ret;
}

template <typename T>
bool cut_solver<T>::consistent(const ineq & i) const {
    // an option could be to check that upper(i.m_poly) <= 0
    bool ret = value(i.m_poly) <= zero_of_type<T>();
    if (!ret) {
        TRACE("cut_solver_state_inconsistent", 
              tout << "inconsistent ineq "; print_ineq(tout,i); tout <<"\n";
              tout << "value = " << value(i.m_poly) << '\n';
              );
        }
    return ret;
}

template <typename T>
int cut_solver<T>::find_non_fixed_var() const {
    // it is a very non efficient implementation for now.
    // the current limitation is that we only deal with bounded vars.
    // the search should be randomized.
    for (unsigned j = 0; j < m_var_infos.size(); j++) {
        const auto & d = m_var_infos[j].m_domain;
        lp_assert(d.lower_bound_exists() && d.upper_bound_exists());
        if (!d.is_fixed())
            return j;
    }
    return -1;
}

template <typename T>
bool cut_solver<T>::all_vars_are_fixed() const {
    return find_non_fixed_var() == -1;
}

template <typename T>
bool cut_solver<T>::consistent() const {
    if (find_non_fixed_var() != -1) {
        // this check could be removed if we use upper bound to check if an inequality holds
        return false; // ignore the variables values and only return true if every variable is fixed
    }
    
    for (const auto & i : m_ineqs) {
        if (!consistent(i))
            return false;
    }
    return true;
}

template <typename T>
unsigned cut_solver<T>::add_ineq_in_local_vars(const std::vector<monomial> & lhs,
                                 const T& free_coeff,
                                 vector<constraint_index> explanation) {
    lp_assert(lhs_is_int(lhs));
    lp_assert(is_int(free_coeff));
    unsigned ineq_index = m_ineqs.size();
        
    m_ineqs.push_back(ineq(lhs, free_coeff, explanation));
    TRACE("ba_int",
          tout << "explanation :";
          for (auto i: explanation) {
              m_print_constraint_function(i, tout);
              tout << "\n";
          });

    for (auto & p : lhs)
        m_var_infos[p.var()].add_dependent_ineq(ineq_index);
        
    return ineq_index;
}

template <typename T>
unsigned cut_solver<T>::add_ineq(const std::vector<monomial> & lhs,
                                 const T& free_coeff,
                                 vector<constraint_index> explanation) {
    lp_assert(lhs_is_int(lhs));
    lp_assert(is_int(free_coeff));
    std::vector<monomial> local_lhs;
    unsigned ineq_index = m_ineqs.size();
    for (auto & p : lhs)
        local_lhs.push_back(monomial(p.coeff(), add_var(p.var())));
        
    m_ineqs.push_back(ineq(local_lhs, free_coeff, explanation));
    TRACE("ba_int",
          tout << "explanation :";
          for (auto i: explanation) {
              m_print_constraint_function(i, tout);
              tout << "\n";
          });

    for (auto & p : local_lhs)
        m_var_infos[p.var()].add_dependent_ineq(ineq_index);
        
    return ineq_index;
}

}

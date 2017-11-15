/*
  Copyright (c) 2017 Microsoft Corporation
  Author: Nikolaj Bjorner, Lev Nachmanson
*/
#pragma once
#include "util/vector.h"
#include "util/trace.h"
#include "util/lp/lp_settings.h"
#include "util/lp/column_namer.h"
#include "util/lp/integer_domain.h"
#include "util/lp/lp_utils.h"
#include <functional>
#include "util/lp/int_set.h"
#include "util/lp/linear_combination_iterator_on_std_vector.h"
#include "util/lp/stacked_vector.h"
#include "util/lp/linear_combination_iterator_on_std_vector.h"
namespace lp {
enum
class lbool { l_false, l_true, l_undef };
template <typename T>
class cut_solver : public column_namer {
public: // for debugging
    class monomial {
        T m_coeff; // the coefficient of the monomial
        var_index m_var; // the variable index
    public:
        monomial(const T& coeff, var_index var) : m_coeff(coeff), m_var(var) {}
        // copy constructor
        monomial(const monomial& m) : monomial(m.coeff(), m.var()) {}
        const T & coeff() const { return m_coeff; }
        T & coeff() { return m_coeff; }
        var_index var() const { return m_var; }
        std::pair<T, var_index> to_pair() const { return std::make_pair(coeff(), var());}
    };

    std::vector<std::pair<T, var_index>> to_pairs(const std::vector<monomial>& ms) const {
        std::vector<std::pair<T, var_index>> ret;
        for (const auto p : ms)
            ret.push_back(p.to_pair());
        return ret;
    }
    
    struct polynomial {
        // the polynomial evaluates to m_coeffs + m_a
        std::vector<monomial> m_coeffs;
        T m_a; // the free coefficient
        polynomial(const std::vector<monomial>& p, const T & a) : m_coeffs(p), m_a(a) {}
        polynomial(const std::vector<monomial>& p) : polynomial(p, 0) {}
        polynomial(): m_a(zero_of_type<T>()) {}
        polynomial(const polynomial & p) : m_coeffs(p.m_coeffs), m_a(p.m_a) {} 
            
        const T & coeff(var_index j) const {
            for (const auto & t : m_coeffs) {
                if (j == t.var()) {
                    return t.coeff();
                }
            }
            return cut_solver::m_local_zero;
        }

        std::vector<monomial> copy_coeff_but_one(var_index j) const {
            std::vector<monomial> ret;
            for (const auto & t : m_coeffs)
                if (t.var() != j)
                    ret.push_back(std::make_pair(t.coeff(), t.var()));

            return ret;
        }

        polynomial &  operator+=(const polynomial & p) {
            m_a += p.m_a;
            for (const auto & t: p.m_coeffs)
                *this += monomial(t.coeff(), t.var());
	    return *this;
	}

        void clear() {
            m_coeffs.clear();
            m_a = zero_of_type<T>();
        }
        
        bool is_zero() const { return m_coeffs.size() == 0 && numeric_traits<T>::is_zero(m_a); }

        unsigned number_of_monomials() const { return m_coeffs.size();}
        
        polynomial & operator+=(const monomial &m ){
            for (unsigned k = 0; k < m_coeffs.size(); k++) {
                auto & l = m_coeffs[k];
                if (m.var() == l.var()) {
                    l.coeff() += m.coeff();
                    if (l.coeff() == 0)
                        m_coeffs.erase(m_coeffs.begin() + k);
                    return *this;
                }
            }
            m_coeffs.push_back(m);
            lp_assert(is_correct());
            return *this;
        }

        bool is_correct() const {
            std::unordered_set<var_index> s;
            for (auto & l : m_coeffs) {
                if (l.coeff() == 0)
                    return false;
                s.insert(l.var());
            }
            return m_coeffs.size() == s.size();
        }
        bool is_tight(unsigned j) const {
            const T & a = coeff(j);
            return a == 1 || a == -1;
        }
    };


        
    
    struct constraint { // we only have less or equal, which is enough for integral variables
        bool m_is_lemma;
        bool m_is_ineq;
        polynomial m_poly;
        T m_d; // the divider for the case of a divisibility constraint
        svector<constraint_index> m_origins; // these indices come from the client
        constraint(const std::vector<monomial>& term, const T& a, bool is_lemma,
                   const svector<constraint_index> & origin):
            m_is_lemma(is_lemma),
            m_is_ineq(true),
            m_poly(term, a),
            m_origins(origin)
        {}

        constraint(const std::vector<monomial>& term, const T& a, const T & divider,
                   bool is_lemma,
                   const svector<constraint_index> & origin):
            m_is_lemma(is_lemma),
            m_is_ineq(false), // it is a division constraint
            m_poly(term, a),
            m_d(divider),
            m_origins(origin)
        {}


        constraint() {}

        bool contains(var_index j) const {
            return m_poly.contains(j);
        }

        const T & coeff(var_index j) const {
            return m_poly.coeff(j);
        }

        const std::vector<monomial>& coeffs() const { return m_poly.m_coeffs;}
        
        void clear() { m_poly.clear(); }

        bool is_simple() const {
            return m_poly.m_coeffs.size() == 1 &&
                (m_poly.m_coeffs[0].coeff() == one_of_type<T>()
                 || m_poly.m_coeffs[0].coeff() == -one_of_type<T>());
        }

        bool is_tight(unsigned j) const {
            const T & a = m_poly.coeff(j);
            return a == 1 || a == -1;
        }
        void add_monomial(const T & t, var_index j) {
            m_poly += monomial(t, j);
        }

        bool is_ineq() const { return m_is_ineq;}
    };

    struct literal {
        bool m_is_decided;
        unsigned m_var;        
        bool m_is_lower;
        T m_bound;
        int m_expl_index; // if m_is_decided, then m_expl_index points to trail, otherwise it points to m_constraints
        polynomial m_tight_ineq;
        literal(bool is_decided, unsigned var_index, bool is_lower, const T & bound, int expl_index):
            m_is_decided(is_decided),
            m_var(var_index),
            m_is_lower(is_lower),
            m_bound(bound),
            m_expl_index(expl_index)
        {
            m_tight_ineq.m_a = zero_of_type<T>();
        }
        literal() {}

        bool is_decided() const { return m_is_decided; }

        bool is_implied() const { return !m_is_decided;}

        static literal make_implied_literal(unsigned var_index, bool is_lower, const T & bound, int expl_index) {
            return literal(false, var_index, is_lower, bound, expl_index);
        }

        static literal make_decided_literal(unsigned var_index, bool is_lower, const T & bound, int expl_index) {
            return literal(true, var_index, is_lower, bound, expl_index);
        }
        unsigned var() const { return m_var; }
    };    

    enum class bound_type {
        LOWER, UPPER, UNDEF
            };
    struct bound_result {
        T m_bound;
        bound_type m_type;
        
        bound_result(const T & b, bound_type bt): m_bound(b), m_type(bt) {}
        bound_result() : m_type(bound_type::UNDEF) {
        }
        void print( std::ostream & out) const {
            if (m_type == bound_type::LOWER) {
                out << "lower bound = ";
            }
            else if (m_type == bound_type::UPPER) {
                out << "upper bound = ";
            }
            else {
                out << "undef";
                return;
            }
            out << m_bound;
        }
    };
    
    struct var_info {
        unsigned m_user_var_index;
        var_info(unsigned user_var_index) : m_user_var_index(user_var_index) {}
        std::vector<unsigned> m_literals; // point to m_trail
        integer_domain<T> m_domain;
        bool is_fixed() const { return m_domain.is_fixed();}
        std::unordered_set<int> m_dependent_constraints; // the set of constraintualities involving the var
        void add_dependent_constraint(unsigned i) {
            m_dependent_constraints.insert(i);
        }
        void remove_depended_constraint(unsigned i) {
            m_dependent_constraints.erase(i);
        }
    };

    std::vector<var_info> m_var_infos;
    
    bool lhs_is_int(const std::vector<monomial> & lhs) const {
        for (auto & p : lhs) {
            if (numeric_traits<T>::is_int(p.coeff()) == false) return false;
        }
        return true;
    }
    
public:
    std::string get_column_name(unsigned j) const {
        return m_var_name_function(m_var_infos[j].m_user_var_index);
    }

    constraint & get_constraint(unsigned i) {
        return m_constraints[i];
    }
    
    const constraint & get_constraint(unsigned i) const {
        return m_constraints[i];
    }

    std::vector<constraint> m_constraints;
    std::vector<T> m_v; // the values of the variables
    std::function<std::string (unsigned)> m_var_name_function;
    std::function<void (unsigned, std::ostream &)> m_print_constraint_function;
    stacked_value<bool> m_decision_has_been_made;  // tracks the number of case splits
    int_set m_changed_vars;
    std::vector<literal>          m_trail;
    lp_settings & m_settings;
    struct scope {
        unsigned m_trail_size;
        unsigned m_constraints_size;
        scope() {}
        scope(const scope& s) : m_trail_size(s.m_trail_size),
                                m_constraints_size(s.m_constraints_size) {}
        scope(unsigned trail_size,
              unsigned constraints_size) : m_trail_size(trail_size),
                                               m_constraints_size(constraints_size) {}

    };

    stacked_value<scope>          m_scope;
    std::unordered_map<unsigned, unsigned> m_user_vars_to_cut_solver_vars;
    static T m_local_zero;
    std::unordered_set<constraint_index> m_explanation; // if this collection is not empty we have a conflict 

    unsigned add_var(unsigned user_var_index) {
        unsigned ret;
        if (try_get_value(m_user_vars_to_cut_solver_vars, user_var_index, ret))
            return ret;
        unsigned j = m_var_infos.size();
        m_var_infos.push_back(var_info(user_var_index));
        return m_user_vars_to_cut_solver_vars[user_var_index] = j;      
    }

    bool is_lower_bound(literal & l) const {
        return l.m_is_lower;
    }
    
    bool lower_for_var(unsigned j, T & lower) const {
        bool ret = false;
        for (unsigned i : m_var_infos[j].m_literals)
            if (is_lower_bound(m_trail[i])) {
                if (ret == false) {
                    ret = true;
                    lower = get_bound(m_trail[i]);
                } else {
                    lower = std::max(lower, get_bound(m_trail[i]));
                }
            }
        return ret;
    }

    bool is_upper_bound(literal & l) const {
        return !l.m_is_upper;
    }
    
    bool upper_for_var(unsigned j, T & upper) const {
        bool ret = false;
        for (unsigned i : m_var_infos[j].m_literals)
            if (is_upper_bound(m_trail[i])) {
                if (ret == false) {
                    ret = true;
                    upper = get_bound(m_trail[i]);
                } else {
                    upper = std::min(upper, get_bound(m_trail[i]));
                }
            }
        return ret;
    }

    // used for testing only
    void add_lower_bound_for_user_var(unsigned user_var_index, const T& bound) {
        unsigned j = m_user_vars_to_cut_solver_vars[user_var_index];
        auto & vi = m_var_infos[j];
        vi.m_domain.intersect_with_lower_bound(bound);
    }

    // used for testing only
    void add_upper_bound_for_user_var(unsigned user_var_index, const T& bound) {
        unsigned j = m_user_vars_to_cut_solver_vars[user_var_index];
        auto & vi = m_var_infos[j];
        vi.m_domain.intersect_with_upper_bound(bound);
    }

    
    bool  at_base_lvl() const { return !m_decision_has_been_made; }

    void simplify_problem() {
        // no-op
    }

    void gc() {
        // no-op
    }

    void restart() {
        // no-op for now
    }

    bool check_inconsistent() {
        // TBD
        return false;
    }

    void cleanup() {  }

    lbool final_check() {
        // there are no more case splits, and all clauses are satisfied.
        // prepare the model for external consumption.
        return lbool::l_true;
    }
    
    bool resolve_conflict_core() {
        // this is where the main action is.
        return true;
    }

    void restrict_var_domain_with_bound_result(var_index j, const bound_result & br) {
        auto & d = m_var_infos[j].m_domain;
        if (br.m_type == bound_type::UPPER) {
            d.intersect_with_upper_bound(br.m_bound);
        } else {
            d.intersect_with_lower_bound(br.m_bound);
        }
        if (d.is_fixed()) {
            if (j >= m_v.size())
                m_v.resize(j + 1);
            d.get_upper_bound(m_v[j]);
        }
    }

    void push_literal_to_trail(const literal & l) {
        unsigned literal_index = m_trail.size();
        m_trail.push_back(l);
        TRACE("ba_int", print_literal(tout, l););
        m_var_infos[l.m_var].m_literals.push_back(literal_index);
    }

    unsigned find_large_enough_j(unsigned i) {
        unsigned r = 0;
        for (const auto & p : m_constraints[i].m_poly.m_coeffs) {
            r = std::max(r, p.var() + 1);
        }
        return r;
    }

    std::string var_name(unsigned j) const {
        return get_column_name(j);
    }

    
    void trace_print_domain_change(std::ostream& out,unsigned j, const T& v, const monomial & p, unsigned constraint_index) const {
        out << "change in domain of " << var_name(j) << ", v = " << v << ", domain becomes ";
        print_var_domain(out, j);
        T lb;
        bool r = lower(constraint_index, lb);
        if (r)
            out << "lower_of_constraint = " << lb << "\n";
        else
            out << "no lower bound for constraint\n";
    }
    
    void propagate_monomial_on_lower(const monomial & p, const T& lower_val, unsigned constraint_index) {
        unsigned j = p.var();
        if (is_pos(p.coeff())) {
            T m;
            get_var_lower_bound(p.var(), m);
            T v = floor(- lower_val / p.coeff()) + m;
            bool change = m_var_infos[j].m_domain.intersect_with_upper_bound(v);
            if (change) {
                TRACE("ba_int", trace_print_domain_change(tout, j, v, p, constraint_index););
                add_bound(v, j, false, constraint_index);
            }
        } else {
            T m;
            get_var_upper_bound(p.var(), m);
            T v = ceil( - lower_val / p.coeff()) + m;
            bool change = m_var_infos[j].m_domain.intersect_with_lower_bound(v);
            if (change) {
                TRACE("ba_int", trace_print_domain_change(tout, j, v, p, constraint_index););
                add_bound(v, j, true , constraint_index);
            }
        }
    }

    void propagate_monomial_on_right_side(const monomial & p, const T& rs, unsigned constraint_index) {
        unsigned j = p.var();
        if (is_pos(p.coeff())) {
            T m;
            T v = floor(rs / p.coeff());
            bool change = m_var_infos[j].m_domain.intersect_with_upper_bound(v);
            if (change) {
                TRACE("ba_int", trace_print_domain_change(tout, j, v, p, constraint_index););
                add_bound(v, j, false, constraint_index);
            }
        } else {
            T v = ceil(rs / p.coeff());
            bool change = m_var_infos[j].m_domain.intersect_with_lower_bound(v);
            if (change) {
                TRACE("ba_int", trace_print_domain_change(tout, j, v, p, constraint_index););
                add_bound(v, j, true , constraint_index);
            }
        }
    }

    void print_var_info(std::ostream & out, const var_info & vi) const {
        out << m_var_name_function(vi.m_user_var_index) << " ";
        print_var_domain(out, vi);
    }

    
    void print_var_domain(std::ostream & out, unsigned j) const {
        m_var_infos[j].m_domain.print(out);
    }

    void print_var_domain(std::ostream & out, const var_info & vi) const {
        vi.m_domain.print(out);
    }

    // b is the value of lower
    void propagate_constraint_on_lower(unsigned i, const T & b) {
        for (const auto & p: m_constraints[i].coeffs()) {
            propagate_monomial_on_lower(p, b, i);
        }
    }

    unsigned find_lower_bound_literal(bool is_lower, unsigned j, unsigned & upper_end_of_trail) const {
        TRACE("ba_int", tout << get_column_name(j) << "\n"; tout << "literal's size = " << m_var_infos[j].m_literals.size() << "\n";);
        for (unsigned k = m_var_infos[j].m_literals.size(); k--;) {
            unsigned literal_index = m_var_infos[j].m_literals[k];
            if (literal_index >= upper_end_of_trail)
                continue;
            const literal& l = m_trail[literal_index];
            if (l.m_var == j && l.m_is_lower == is_lower) {
                TRACE("ba_int",
                      tout << "found lower bound expl\n";
                      print_literal(tout, l); tout << "\n";);
                return literal_index;
            }
        }
        lp_assert(false); // unreachable
        return 0;// to avoid the warning
    }

    void add_constraint_origins_to_explanation(unsigned constraint_index) {
        for (auto j : m_constraints[constraint_index].m_origins)
            m_explanation.insert(j);
    }

    void fill_conflict_explanation(unsigned constraint_index, unsigned upper_end_of_trail) {
        // it is a depth search in the DAG of constraint: the chidlren of a constraint are those constraints that provide its lower bound
        add_constraint_origins_to_explanation(constraint_index);
        const constraint& in = m_constraints[constraint_index];
        TRACE("ba_int", print_constraint(tout, in););
        for (const auto & p: in.coeffs()){
            unsigned literal_index = find_lower_bound_literal(is_pos(p.coeff()), p.var(), upper_end_of_trail);
            unsigned l_constraint_index = m_trail[literal_index].m_expl_index;
            if (!m_constraints[l_constraint_index].is_simple()) 
                fill_conflict_explanation(l_constraint_index, literal_index);
            else
                add_constraint_origins_to_explanation(l_constraint_index);
        }
    }

    void trace_print_constraint(std::ostream& out, unsigned i) {
        print_constraint(out, i); out << "\n";
        unsigned j;
        auto pairs = to_pairs(m_constraints[i].m_poly.m_coeffs);
        auto it = linear_combination_iterator_on_std_vector<T>(pairs);
        while (it.next(j)) {
            out << "domain of " << var_name(j) << " = ";
            print_var_domain(out, j);
        }
    }

   
    void propagate_constraint_only_one_unlim(unsigned constraint_index, unsigned the_only_unlim) {
        const constraint& i = m_constraints[constraint_index];
        T rs = - i.m_poly.m_a;
        for (unsigned j = 0; j < i.m_poly.m_coeffs.size(); j++) {
            if (j == the_only_unlim) continue;
            T m;
            lower_monomial(i.m_poly.m_coeffs[j], m);
            rs -= m;
        }

        // we cannot get a conflict here because the monomial i.m_poly.m_coeffs[the_only_unlim]
        // is unlimited from below and we are adding an upper bound for it
        propagate_monomial_on_right_side(i.m_poly.m_coeffs[the_only_unlim], rs, constraint_index);
    }

    bool conflict() const { return m_explanation.size() > 0; }

    // returns -1 if there is no conflict and the index of the conflict constraint otherwise
    int  propagate_on_constraints_of_var(var_index j) {
        for (unsigned i : m_var_infos[j].m_dependent_constraints) {
            if (!propagate_constraint(i))
                return i;
        }
        return -1;
    }
    
    // returns -1 if there is no conflict and the index of the conflict constraint otherwise
    int propagate_constraints_for_changed_vars() {
        TRACE("cut_solver_state", tout << "changed vars size = " << m_changed_vars.size() << "\n";);
        while (!m_changed_vars.is_empty()) {
            unsigned j = m_changed_vars.m_index.back();
            int i = propagate_on_constraints_of_var(j);
            if (i >= 0) {
                m_changed_vars.clear();
                return i;
            }
            m_changed_vars.erase(j);
        }
        return -1;
    }

    bool var_lower_bound_exists(var_info i) const {
        const var_info & v = m_var_infos[i];
        return v.m_domain.lower_bound_exists();
    }
    
    bool get_var_lower_bound(var_index i, T & bound) const {
        const var_info & v = m_var_infos[i];
        return v.m_domain.get_lower_bound(bound);
    }

    bool get_var_upper_bound(var_index i, T & bound) const {
        const var_info & v = m_var_infos[i];
        return v.m_domain.get_upper_bound(bound);
    }

    bool lower_monomial_exists(const monomial & p) const {
        lp_assert(p.coeff() != 0);

        if (p.coeff() > 0) {
            if (!m_var_infos[p.var()].m_domain.lower_bound_exists())
                return false;
        }
        else {
            if (!m_var_infos[p.var()].m_domain.upper_bound_exists())
                return false;
        }
        return true;
    }

    bool upper_monomial_exists(const monomial & p) const {
        lp_assert(p.coeff() != 0);
        if (p.coeff() > 0) {
            if (!m_var_infos[p.var()].m_domain.upper_bound_exists())
                return false;
        }
        else {
            if (!m_var_infos[p.var()].m_domain.lower_bound_exists())
                return false;
        }
        return true;
    }

    
    // finds the lower bound of the monomial,
    // otherwise returns false
    bool lower_monomial(const monomial & p, T & lb) const {
        lp_assert(p.coeff() != 0);
        T var_bound;
        if (p.coeff() > 0) {
            if (!get_var_lower_bound(p.var(), var_bound))
                return false;
            lb = p.coeff() * var_bound;
        }
        else {
            if (!get_var_upper_bound(p.var(), var_bound))
                return false;
            lb = p.coeff() * var_bound;
        }
        return true;
    }

    bool upper_monomial(const monomial & p, T & lb) const {
        lp_assert(p.coeff() != 0);
        T var_bound;
        if (p.coeff() > 0) {
            if (!get_var_upper_bound(p.var(), var_bound))
                return false;
        }
        else {
            if (!get_var_lower_bound(p.var(), var_bound))
                return false;
        }
        lb = p.coeff() * var_bound;
        return true;
    }

    
    
    // returns false if not limited from below
    // otherwise the answer is put into lb
    bool lower(unsigned constraint_index, T & lb) const {
        return lower(m_constraints[constraint_index].m_poly, lb);
    }

    // returns false if not limited from below
    // otherwise the answer is put into lb
    T lower_val(unsigned constraint_index) const {
        return lower_val(m_constraints[constraint_index]);
    }
    
    T lower_val(const constraint & i) const {
        T lb;
#if Z3DEBUG
        bool r =
#endif
            lower(i.m_poly, lb);
        lp_assert(r);
        return lb;
    }

    void pop() { pop(1); }
        
    // returns false if not limited from below
    // otherwise the answer is put into lb
    bool lower(const polynomial & f, T & lb) const {
        lb = f.m_a;
        T lm;
        for (const auto & p : f.m_coeffs) {
            if (lower_monomial(p, lm)) {
                lb += lm;
            } else {
                return false;
            }
        }
        return true;
    }


    // returns the number of lower unlimited monomials - 0, 1, >=2
    // if there is only one lower unlimited then the index is put into the_only_unlimited
    int lower_analize(const constraint & f, unsigned & the_only_unlimited) const {
        int ret = 0;
        for (unsigned j = 0; j < f.m_poly.m_coeffs.size(); j++) {
            if (!lower_monomial_exists(f.m_poly.m_coeffs[j])) {
                if (ret == 1)
                    return 2;
                ret ++;
                the_only_unlimited = j;
            }
        }
        return ret;
    }


 
    bound_result lower_without(const polynomial & p, var_index j) const {
        for (const auto & t:  p.m_coeffs) {
            if (t.var() == j)
                continue;
            if (!lower_monomial_exists(t)) {
                return bound_result();
            }
        }
        // if we are here then there is a lower bound for p
        T bound = p.m_a;
        for (const auto & t:  p.m_coeffs) {
            if (t.var() == j)
                continue;

            T l;
            lower_monomial(t, l);
            bound += l;
        }
        return bound_result(bound,bound_type::LOWER);
    }

    bound_result upper_without(const polynomial & p, var_index j) const {
        for (const auto & t:  p.m_coeffs) {
            if (t.var() == j)
                continue;
            if (!upper_monomial_exists(t)) {
                return bound_result();
            }
        }
        // if we are here then there is an upper bound for p
        T bound = p.m_a;
        for (const auto & t:  p.m_coeffs) {
            if (t.var() == j)
                continue;
            T b;
            upper_monomial(t, b);
            bound += b;
        }
        return bound_result(bound, bound_type::UPPER);
    }

    bool upper(const polynomial & p, T b) const {
        for (const auto & t:  p.m_coeffs) {
            if (!upper_monomial_exists(t)) {
                return false;
            }
        }
        // if we are here then there is an upper bound for p
        b = p.m_a;
        T bb;
        for (const auto & t:  p.m_coeffs) {
            upper_monomial(t, bb);
            b += bb;
        }
        return true;
    }

    
    
    
    // a is the coefficient before j
    bound_result bound_on_polynomial(const polynomial & p, const T& a, var_index j) const {
        lp_assert(!is_zero(a));
        if (numeric_traits<T>::is_pos(a)) {
            bound_result r = upper_without(p, j);
            if (r.m_type == bound_type::UNDEF)
                return r;
            lp_assert(is_int(r.m_bound));
            r.m_bound = - ceil_ratio(r.m_bound , a);
            r.m_type = bound_type::UPPER;
            return r;
        }
        else {
            bound_result r = lower_without(p, j);
            if (r.m_type == bound_type::UNDEF)
                return r;
            r.m_bound = -floor_ratio(r.m_bound, a);
            r.m_type = bound_type::LOWER;
            return r;
        }
    }
    

    
    bound_result bound(const constraint & q, var_index j) const {
        const T& a = q.m_poly.coeff(j);
        return bound_on_polynomial(q.m_poly, a, j);
    }

    bound_result bound(unsigned constraint_index, var_index j) const {
        return bound(m_constraints[constraint_index], j);
    }

    
    void print_constraint(std::ostream & out, unsigned i)  const {
        print_constraint(out, m_constraints[i]);
    }

    void print_constraint(std::ostream & out, const constraint & i ) const {
        out << (i.m_is_lemma? "lemma ": "assert ");
        if (!i.m_is_ineq) {
            out << i.m_d << " | ";
        }
        print_polynomial(out, i.m_poly) ;
        if (i.m_is_ineq) {
            out << " <= 0";
        }
        out << "\n";
    }

    void print_literal(std::ostream & out, const literal & t) const {
        out << (t.m_is_decided? "decided ": "implied ");
        out << get_column_name(t.m_var) << " ";
        if (t.m_is_lower)
            out << ">= ";
        else
            out << "<= ";
        out << t.m_bound;

        if (t.m_is_decided == false) {
            out << " by constraint " << t.m_expl_index << " "; print_constraint(out, m_constraints[t.m_expl_index]);
        }
        
     }
    

    
    void print_polynomial(std::ostream & out, const polynomial & p) const {
        std::vector<std::pair<T, unsigned>> pairs = to_pairs(p.m_coeffs);
        this->print_linear_combination_of_column_indices_std(pairs, out);
        if (!is_zero(p.m_a)) {
            if (p.m_a < 0) {
                out << " - " << -p.m_a;
            } else {
                out << " + " << p.m_a;
            }
        }
    }

    // Trying to improve constraint "ie" by eliminating var j by using a  tight inequality 
    // for j. The left side of the inequality is passed as a parameter.
    bool resolve(const polynomial & ie, unsigned j, bool sign_j_in_ti_is_pos, const polynomial & ti, polynomial & result) const {
        lp_assert(ti.is_tight(j));
        lp_assert(result.is_zero());
        result.clear();
        const auto &coeffs = ie.m_coeffs;
        // todo: implement a more efficient version
        bool found = false;
        T a;
        for (const auto & c : coeffs) {
            if (c.var() == j) {
                a = c.coeff();
                found = true;
            }
            else {
                result.m_coeffs.push_back(c);
            }
        }

        if (!found) {
            result.m_a = ie.m_a;
            return false;
        }
        
        if (is_neg(a)) {
            a = -a;
            if (!sign_j_in_ti_is_pos) return false;
        } else {
            if (sign_j_in_ti_is_pos) return false;
        }

        for (auto & c : ti.m_coeffs) {
            if (c.var() != j)
                result += monomial(a * c.coeff(), c.var());
            else {
                lp_assert(c.coeff() == (sign_j_in_ti_is_pos? one_of_type<T>():-one_of_type<T>()));
            }
        }
        result.m_a = ie.m_a + a * ti.m_a;
        return true;
    }

    // returns true iff p imposes a better bound on j
    bool improves(var_index j, const constraint & p) const {
        auto a = p.coeff(j);
        if (is_zero(a))
            return false;
        const auto& dom = m_var_infos[j].m_domain;
        if (dom.is_empty())
            return false;
        if (is_pos(a)) {
            bound_result new_upper = bound(p, j);
            if (new_upper.m_type == bound_type::UNDEF)
                return false;
            T b;
            bool upper_bound_exists = get_var_upper_bound(j, b);
            return (!upper_bound_exists || new_upper.m_bound < b) &&
                dom.intersection_with_upper_bound_is_empty(new_upper.m_bound);
        }

        lp_assert(is_neg(a));
        bound_result new_lower = bound(p, j);
        if (new_lower.m_type == bound_type::UNDEF)
            return false;
        T b;
        bool lower_bound_exists = get_var_lower_bound(j, b);
        return (!lower_bound_exists || new_lower.m_bound > b) &&
            dom.intersection_with_lower_bound_is_empty(new_lower.m_bound);
    }


    // returns true iff br imposes a better bound on j
    bool improves(var_index j, const bound_result & br) const {
        if (br.m_type == bound_type::UNDEF)
            return false;
        const auto& dom = m_var_infos[j].m_domain;
        if (dom.is_empty())
            return false;
        if (br.m_type == bound_type::UPPER) {
            T b;
            bool j_has_upper_bound = get_var_upper_bound(j, b);
            return (!j_has_upper_bound || br.m_bound < b) &&
                !dom.intersection_with_upper_bound_is_empty(br.m_bound);
        }

        if (br.m_type == bound_type::UNDEF)
            return false;
        T b;
        bool lower_bound_exists = get_var_lower_bound(j, b);
        return (!lower_bound_exists || br.m_bound > b) &&
            !dom.intersection_with_lower_bound_is_empty(br.m_bound);
    }


    void add_bound(T v, unsigned j, bool is_lower, unsigned constraint_index) {
        push_literal_to_trail(literal::make_implied_literal(j, is_lower, v, constraint_index));
    }
    
    bool literal_is_correct(const literal &t ) const {
        if (t.is_decided())
            return true;
        auto & i = m_constraints[t.m_constraint_index];
        int sign_should_be = t.m_is_lower? -1: 1;
        const T &a = i.coeff(t.m_var);
        int sign = a > 0? 1: -1;
        return sign == sign_should_be;
    }

    T value(const polynomial & p) const {
        T ret= p.m_a;
        for (const auto & t:p.m_coeffs)
            ret += t.coeff() * m_v[t.var()];
        return ret;
    }

    void pop_constraints() {
        for (unsigned j = m_constraints.size(); j-- > m_scope().m_constraints_size;) {
            const constraint & i = m_constraints[j];
            for (const auto & p: i.m_poly.m_coeffs) {
                m_var_infos[p.var()].remove_depended_constraint(j);
            }
            m_constraints.pop_back();
        }
        lp_assert(m_constraints.size() == m_scope().m_constraints_size);
    }

    void pop_var_domains(unsigned k) {
        for (auto & v : m_var_infos) {
            v.m_domain.pop(k);
        }
    }

    void push_var_domains() {
        for (auto & v : m_var_infos) {
            v.m_domain.push();
        }
    }

    bool flip_coin() {
        return m_settings.random_next() % 2 == 0;
    }
    

    lbool check() {
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

    void init_search() {
        lp_assert(m_explanation.size() == 0);
        lp_assert(m_changed_vars.size() == 0);
    }


    // returns true if there is no conflict and false otherwise
    bool propagate_constraint(unsigned i) {
        TRACE("ba_int", trace_print_constraint(tout, i););
        const constraint & in = m_constraints[i];
        if (in.is_simple())
            return true;
        // consider a special case for constraint with just two variables
        unsigned the_only_unlim;
        int r = lower_analize(in, the_only_unlim);
        if (r == 0) {
            T b;
            lower(in.m_poly, b);
            if (is_pos(b)) {
                TRACE("cs_inconsistent", trace_print_constraint(tout, i););
                return false;
            } else {
                propagate_constraint_on_lower(i, b);
            }
        } else if (r == 1) {
            propagate_constraint_only_one_unlim(i, the_only_unlim);
        }
        return true;
    }
    void print_trail(std::ostream & out) const {
        for (const auto & l : m_trail) {
            print_literal(out, l);
            out << "\n";
        }
    }
    void print_state(std::ostream & out) const {
        out << "constraints:\n";
        for (const auto & i: m_constraints) {
            print_constraint(out, i);
        }
        out << "end of constraints\n";
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
    lbool bounded_search() {
        lbool is_sat = propagate_and_backjump_step();
        if (is_sat != lbool::l_undef)
            return is_sat;

        gc();

        if (!decide()) {
            lbool is_sat = final_check();
            if (is_sat != lbool::l_undef) {
                return is_sat;
            }
        }
        return lbool::l_undef;
    }

    bool decide() {
        int j = find_non_fixed_var();
        if (j < 0)
            return false;
        TRACE("decide", tout << "going to decide " << var_name(j) << " var index = " << j << "\n";
              tout << "domain = "; print_var_domain(tout, j); tout << "\n";);
        m_changed_vars.insert(j);
        decide_var_on_bound(j, flip_coin());
        return true;
    }

    lbool propagate_and_backjump_step() {
        int incostistent_constraint = propagate();
        TRACE("cs_dec", tout << "trail = \n"; print_trail(tout); tout << "end of trail\n";);

        if (incostistent_constraint >= 0) {
            if (at_base_lvl()) {
                fill_conflict_explanation(incostistent_constraint, m_trail.size());
                return lbool::l_false;
            }
            resolve_conflict(incostistent_constraint);
        }

        if(!all_vars_are_fixed())
            return lbool::l_undef;

        return lbool::l_true;
    }

    // returns -1 if consistent, or the index of an incostistent constraint
    int inconsistent() const {
        for (int i=0;i<m_constraints.size();i++)
            if (lower(i) > 0)
                return i;

        return -1;
    }
    
    bool decision_is_redundant_for_constraint(const polynomial& i, const literal & l) const {
        const T & coeff = i.coeff(l.m_var);
        if (is_zero(coeff))
            return true;
        if (is_pos(coeff)) {
            return l.m_is_lower;
        }
        return !l.m_is_lower;
    }

    bool is_divizible(const T & a, const T & b) const {
        lp_assert(!is_zero(b));
        return is_zero(a % b);
    }
    
    void create_div_ndiv_parts_for_tightening(const polynomial & p, const T & coeff, polynomial & div_part, polynomial & ndiv_part) {
        for (const auto &m : p.m_coeffs) {
            if (is_divizible(m.coeff(), coeff)){
                div_part.m_coeffs.push_back(m);
            } else {
                ndiv_part.m_coeffs.push_back(m);
            }
        }

        TRACE("tight",
              tout << "div_part = ";
              print_polynomial(tout, div_part);
              tout << "\nndiv_part = ";
              print_polynomial(tout, ndiv_part););
    }

    void decide_lower_neg(polynomial &ndiv_part, T & c, literal &l) {
        ndiv_part += -c * m_trail[l.m_expl_index].m_tight_ineq;
        lp_assert(is_zero(ndiv_part.coeff(l.var())));
        TRACE("tight", tout << "Decided-Lower-Neg, "; tout << "ndiv_part = ";
              print_polynomial(tout, ndiv_part); tout << "\n";);
    }

    void decided_lower(const T & j_coeff, const T & c, polynomial &div_part, polynomial &ndiv_part, const literal &l) {
        T k = is_pos(j_coeff)?ceil( c / j_coeff): floor(c / j_coeff);
        ndiv_part += monomial(-c, l.var()); // it will be moved to div_part
        T j_coeff_by_k = j_coeff * k;
        T m = j_coeff_by_k - c;
        TRACE("tight", tout << "c = " << c << ", j_coeff = " << j_coeff <<
              ", c / j_coeff = " << c/j_coeff << ", k = " <<
              k << ", j_coeff * k = " << j_coeff * k << ", m = " << m << "\n"; );  
        lp_assert(is_pos(m));
        for (const monomial & n : l.m_tight_ineq.m_coeffs) {
            if (n.var() == l.var()) {
                lp_assert(n.coeff() == one_of_type<T>());
                div_part += monomial(j_coeff_by_k, l.var());
            } else {
                ndiv_part += monomial(m*n.coeff(), n.var());
            }
        }
        TRACE("tight", tout << "Decided-Lower ";
              tout << "div_part = ";
              print_polynomial(tout, div_part);
              tout << "\n";
              tout << "ndiv_part = ";
              print_polynomial(tout, ndiv_part);
              tout << "\n";);
                    
        lp_assert(false); // just to stop here
    }
    
    void tighten_on_prev_literal(polynomial & p, const T & j_coeff, polynomial & div_part, polynomial &ndiv_part, int trail_index) {
        TRACE("tight",
              tout << "trail_index = " << trail_index << ", ";
              print_literal(tout, m_trail[trail_index]););
        literal & l = m_trail[trail_index];
        if (l.m_tight_ineq.number_of_monomials() == 0) {
            create_tight_ineq_under_literal(trail_index);
        }
        if (l.is_implied()) { // Resolve-Implied
            polynomial result;
            resolve(ndiv_part, l.m_var, !l.m_is_lower, l.m_tight_ineq, result);
            ndiv_part = result;
            TRACE("tight", tout << "after resolve ndiv_part = "; print_polynomial(tout, ndiv_part);
                  tout << "\n";);
        } else { 
            lp_assert(l.is_decided());
            create_tight_ineq_under_literal(l.m_expl_index);
            T c = ndiv_part.coeff(l.var());
            if (l.m_is_lower) {
                if (is_neg(c)) {
                    decide_lower_neg(ndiv_part, c, l);
                } else { 
                    decided_lower(j_coeff, c, div_part, ndiv_part, l);
                }
            } else {
                lp_assert(!l.m_is_lower);
                if (is_pos(c)) { // Decided-Upper-Pos
                    lp_assert(false);
                } else { // Decided-Upper
                    lp_assert(false);
                }
            }
        }
        
    }
    
    // see page 88
    void tighten(polynomial & p, unsigned j_of_var, const T& j_coeff, unsigned trail_index) {
        polynomial div_part, ndiv_part;
        ndiv_part.m_a = p.m_a;
        TRACE("tight",
              tout << "trail_index = " << trail_index;
              tout << ", var name = " << var_name(j_of_var) << ";";
              tout << "p = ";
              print_polynomial(tout, p););
        create_div_ndiv_parts_for_tightening(p, j_coeff, div_part, ndiv_part);
        int k = trail_index - 1;
        lp_assert(k >= 0);
        while (ndiv_part.number_of_monomials() > 0) {
            tighten_on_prev_literal(p, j_coeff, div_part, ndiv_part, k--);
        }
        T abs_j_coeff = abs(j_coeff);
        p.clear();
        for (const auto & m : div_part.m_coeffs) {
            p.m_coeffs.push_back(monomial(m.coeff() / abs_j_coeff, m.var()));
        }
        p.m_a = ceil(ndiv_part.m_a / abs_j_coeff);
        TRACE("tight", tout << "trail_index = " << trail_index << ", got tight p = "; print_polynomial(tout, p); tout << "\n";);
    }
    
    void create_tight_ineq_under_literal(unsigned trail_index) {
        TRACE("tight", tout << "trail_index = " << trail_index << "\n";);
        literal & l = m_trail[trail_index];
        if (l.m_tight_ineq.number_of_monomials() > 0)
            return;
        const constraint & c = m_constraints[l.m_expl_index];
        lp_assert(c.m_is_ineq);
        polynomial &p = l.m_tight_ineq = c.m_poly;
        unsigned j = l.m_var;
        const T& a = p.coeff(j);
        lp_assert(!is_zero(a));
        if (a == one_of_type<T>() || a == - one_of_type<T>()) {
            TRACE("tight", tout << "created = "; print_polynomial(tout, l.m_tight_ineq););
            return;
        }
        tighten(p, j, a, trail_index);
    }

    bool resolve_conflict_for_inequality_on_trail(polynomial & p, unsigned trail_index) {
        const literal & l = m_trail[trail_index];
        if (l.is_decided()) {
            if (decision_is_redundant_for_constraint(p, l))
                pop(); // skip decision
            else
                lp_assert(false); // not implemented
            return false;
        } else { // the literal is implied
            create_tight_ineq_under_literal(trail_index);
        }
        return false;
    }
    
    bool resolve_conflict_for_inequality(polynomial & p) {
#if Z3DEBUG
        T b;
        lp_assert(lower(p, b) && is_pos(b));
#endif
        bool done = false;
        unsigned j = m_trail.size() - 1;
        while (!done) {
            done = resolve_conflict_for_inequality_on_trail(p, j--);
            if (j >= m_trail.size()) {
                lp_assert(m_trail.size());
                j = m_trail.size() - 1;
            }
        }
        return false;
    }

    bool resolve_conflict(int inconstistent_constraint) {
        lp_assert(!at_base_lvl());
        constraint & i = m_constraints[inconstistent_constraint];
        if (i.is_ineq()) {
            return resolve_conflict_for_inequality(i.m_poly);
        } else {
            lp_assert(false); // not implemented
            return false;
        }
    
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

    void print_scope(std::ostream& out) const {
        out << "trail_size = " << m_scope().m_trail_size << ", constraints_size = " << m_scope().m_constraints_size << "\n";
    }

    void pop(unsigned k) {
        TRACE("trace_push_pop_in_cut_solver", tout << "before pop\n";print_state(tout););
        m_scope.pop(k);
        TRACE("trace_push_pop_in_cut_solver", tout << "scope = ";print_scope(tout); tout << "\n";);
        m_trail.resize(m_scope().m_trail_size);
        pop_constraints();
        pop_var_domains(k);
        m_decision_has_been_made.pop(k);
        TRACE("trace_push_pop_in_cut_solver", tout << "after pop\n";print_state(tout););
    }

    void push() {
        TRACE("trace_push_pop_in_cut_solver", print_state(tout););
        m_scope = scope(m_trail.size(), m_constraints.size());
        m_scope.push();
        push_var_domains();
        m_decision_has_been_made.push();
    }

    cut_solver(std::function<std::string (unsigned)> var_name_function,
               std::function<void (unsigned, std::ostream &)> print_constraint_function,
               lp_settings & settings
               ) : m_var_name_function(var_name_function),
                   m_print_constraint_function(print_constraint_function),
                   m_decision_has_been_made(false),
                   m_settings(settings)
    {}

    // returns -1 if there is no conflict and the index of the conflict constraint otherwise
    int propagate() {
        propagate_simple_constraints();
        return propagate_constraints_for_changed_vars();
    }

    // walk the trail backward and find the last implied bound on j (of the right kind)
    unsigned find_literal_index(unsigned j, bool is_lower) const {
        for (unsigned k = m_trail.size(); k-- > 0;){
            const auto & l = m_trail[k];
            if (!l.m_is_decided && l.m_var == j && l.m_is_lower == is_lower)
                return k;
        }
        TRACE("find_literal", tout << "cannot find deciding literal for " << var_name(j)<<  " j = " << j << " is_lower = " << is_lower << std::endl;);
        lp_assert(false); // unreacheable
        return 0;
    }
    
    void decide_var_on_bound(unsigned j, bool decide_on_lower) {
        push();    
        T b;
        std::vector<monomial> lhs;
        if (decide_on_lower) {
            m_var_infos[j].m_domain.get_lower_bound(b);
            m_var_infos[j].m_domain.intersect_with_upper_bound(b);
        }
        else {
            m_var_infos[j].m_domain.get_upper_bound(b);
            m_var_infos[j].m_domain.intersect_with_lower_bound(b);
        }
        m_decision_has_been_made = true;
        m_trail.push_back(literal::make_decided_literal(j, !decide_on_lower, b, find_literal_index(j, decide_on_lower)));
    }

    bool propagate_simple_constraint(unsigned constraint_index) {
        const constraint & t = m_constraints[constraint_index];
        TRACE("cut_solver_state_simpe_constraint",   print_constraint(tout, t); tout << std::endl;);
        var_index j = t.m_poly.m_coeffs[0].var();
        
        bound_result br = bound(t, j);
        TRACE("cut_solver_state_simpe_constraint", tout << "bound result = {"; br.print(tout); tout << "}\n";
              tout << "domain of " << get_column_name(j) << " = "; m_var_infos[j].m_domain.print(tout);
              tout << "\n";
              );
        
        if (improves(j, br)) {
            literal l = literal::make_implied_literal(j, br.m_type == bound_type::LOWER, br.m_bound, constraint_index);
            push_literal_to_trail(l);
            restrict_var_domain_with_bound_result(j, br);
            TRACE("cut_solver_state_simpe_constraint", tout <<"improved domain = ";
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
    

    bool propagate_simple_constraints() {
        bool ret = false;
        for (unsigned i = 0; i < m_constraints.size(); i++) {
            if (m_constraints[i].is_simple() && propagate_simple_constraint(i)){
                ret = true;
            }
        }
        return ret;
    }

    bool consistent(const constraint & i) const {
        // an option could be to check that upper(i.m_poly) <= 0
        bool ret = value(i.m_poly) <= zero_of_type<T>();
        if (!ret) {
            TRACE("cut_solver_state_inconsistent", 
                  tout << "inconsistent constraint "; print_constraint(tout,i); tout <<"\n";
                  tout << "value = " << value(i.m_poly) << '\n';
                  );
        }
        return ret;
    }

    int find_non_fixed_var() const {
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

    bool all_vars_are_fixed() const {
        return find_non_fixed_var() == -1;
    }

    bool consistent() const {
        if (find_non_fixed_var() != -1) {
            // this check could be removed if we use upper bound to check if an constraint holds
            return false; // ignore the variables values and only return true if every variable is fixed
        }
    
        for (const auto & i : m_constraints) {
            if (!consistent(i))
                return false;
        }
        return true;
    }

    unsigned add_ineq(const std::vector<monomial> & lhs,
                      const T& free_coeff,
                      svector<constraint_index> origins) {
        lp_assert(lhs_is_int(lhs));
        lp_assert(is_int(free_coeff));
        std::vector<monomial> local_lhs;
        unsigned constraint_index = m_constraints.size();
        for (auto & p : lhs)
            local_lhs.push_back(monomial(p.coeff(), add_var(p.var())));
        
        m_constraints.push_back(constraint(local_lhs, free_coeff, false, origins)); // false means it is not a lemma
        TRACE("ba_int",
              tout << "explanation :";
              for (auto i: origins) {
                  m_print_constraint_function(i, tout);
                  tout << "\n";
              });

        for (auto & p : local_lhs)
            m_var_infos[p.var()].add_dependent_constraint(constraint_index);
        
        return constraint_index;
    }
};
template <typename T>
inline typename cut_solver<T>::polynomial operator*(const T & a, const typename cut_solver<T>::polynomial & p) {
    typename cut_solver<T>::polynomial ret;
    ret.m_a = p.m_a * a;
    
    for (const auto & t: p.m_coeffs)
        ret.m_coeffs.emplace_back(a * t.coeff(), t.var());
    
    return ret;
}

}

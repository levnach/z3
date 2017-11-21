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

class cut_solver : public column_namer {
public: // for debugging
    struct monomial {
        mpq           m_coeff; // the coefficient of the monomial
        var_index   m_var; // the variable index
    public:
        monomial(const mpq& coeff, var_index var) : m_coeff(coeff), m_var(var) {}
        // copy constructor
        monomial(const monomial& m) : monomial(m.coeff(), m.var()) {}
        monomial(var_index var) : monomial(one_of_type<mpq>(), var) {}
        const mpq & coeff() const { return m_coeff; }
        mpq & coeff() { return m_coeff; }
        var_index var() const { return m_var; }
        std::pair<mpq, var_index> to_pair() const { return std::make_pair(coeff(), var());}
    };

    std::vector<std::pair<mpq, var_index>> to_pairs(const std::vector<monomial>& ms) const {
        std::vector<std::pair<mpq, var_index>> ret;
        for (const auto p : ms)
            ret.push_back(p.to_pair());
        return ret;
    }
    
    struct polynomial {
        // the polynomial evaluates to m_coeffs + m_a
        std::vector<monomial> m_coeffs;
        mpq                     m_a; // the free coefficient
        polynomial(const std::vector<monomial>& p, const mpq & a) : m_coeffs(p), m_a(a) {}
        polynomial(const std::vector<monomial>& p) : polynomial(p, 0) {}
        polynomial(): m_a(zero_of_type<mpq>()) {}
        polynomial(const polynomial & p) : m_coeffs(p.m_coeffs), m_a(p.m_a) {} 
            
        const mpq & coeff(var_index j) const {
            for (const auto & t : m_coeffs) {
                if (j == t.var()) {
                    return t.coeff();
                }
            }
            return m_local_zero;
        }

        polynomial &  operator+=(const polynomial & p) {
            m_a += p.m_a;
            for (const auto & t: p.m_coeffs)
                *this += monomial(t.coeff(), t.var());
            return *this;
        }

        void add(const mpq & c, const polynomial &p) {
            m_a += p.m_a * c;
            
            for (const auto & t: p.m_coeffs)
                *this += monomial(c * t.coeff(), t.var());
        }
        
        void clear() {
            m_coeffs.clear();
            m_a = zero_of_type<mpq>();
        }
        
        bool is_empty() const { return m_coeffs.size() == 0 && numeric_traits<mpq>::is_zero(m_a); }

        unsigned number_of_monomials() const { return m_coeffs.size();}

        void add(const monomial &m ){
            if (is_zero(m.coeff())) return;
            for (unsigned k = 0; k < m_coeffs.size(); k++) {
                auto & l = m_coeffs[k];
                if (m.var() == l.var()) {
                    l.coeff() += m.coeff();
                    if (l.coeff() == 0)
                        m_coeffs.erase(m_coeffs.begin() + k);
                    return;
                }
            }
            m_coeffs.push_back(m);
            lp_assert(is_correct());
        }
        
        polynomial & operator+=(const monomial &m ){
            add(m);
            return *this;
        }

        polynomial & operator+=(const mpq &c ){
            m_a += c;
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
            const mpq & a = coeff(j);
            return a == 1 || a == -1;
        }
    };
    
    class constraint { // we only have less or equal for the inequality sign, which is enough for integral variables
        bool                        m_is_lemma;
        bool                        m_is_ineq;
        const polynomial            m_poly;
        mpq                         m_d; // the divider for the case of a divisibility constraint
        const svector<constraint_index>   m_origins; // these indices come from the client
    public :
        const polynomial & poly() const { return m_poly; }
        const svector<constraint_index> & origins() const { return m_origins;}
        bool is_lemma() const { return m_is_lemma; }
        bool is_ineq() const { return m_is_ineq; }
        const mpq & divider() const { return m_d; }

    public :
        static constraint make_ineq_assert(
            const std::vector<monomial>& term,
            const mpq& a,
            const svector<constraint_index> & origin) {
            return constraint(origin, polynomial(term, a), false, true);
        }
    
    public :
        static constraint make_ineq_lemma(const polynomial &p) {
            return constraint(svector<constraint_index>(), p, true, true);
        }

        static constraint make_div_lemma(const polynomial &p, const mpq & div) {
            auto c =  constraint(svector<constraint_index>(), p, true, true);
            c.m_d = div;
            return c;
        }

        
        constraint(const svector<constraint_index> & origin,
                   const polynomial & p,
                   bool is_lemma,
                   bool is_ineq):
            m_is_lemma(is_lemma),
            m_is_ineq(is_ineq),
            m_poly(p),
            m_origins(origin)
        {}
    public:
        constraint() {}

        const mpq & coeff(var_index j) const {
            return m_poly.coeff(j);
        }
        const std::vector<monomial>& coeffs() const { return m_poly.m_coeffs;}
        bool is_simple() const {
            return m_poly.m_coeffs.size() == 1 &&
                (m_poly.m_coeffs[0].coeff() == one_of_type<mpq>()
                 || m_poly.m_coeffs[0].coeff() == -one_of_type<mpq>());
        }

        bool is_tight(unsigned j) const {
            const mpq & a = m_poly.coeff(j);
            return a == 1 || a == -1;
        }
        
    };

    class literal {
        bool          m_is_decided;
        unsigned      m_var;        
        bool          m_is_lower;
        mpq           m_bound;
        int           m_expl_index; // if m_is_decided, then m_expl_index points to trail, otherwise it points to m_constraints
        polynomial    m_tight_ineq;
        literal(bool is_decided, unsigned var_index, bool is_lower, const mpq & bound, int expl_index):
            m_is_decided(is_decided),
            m_var(var_index),
            m_is_lower(is_lower),
            m_bound(bound),
            m_expl_index(expl_index)
        {
            m_tight_ineq.m_a = zero_of_type<mpq>();
        }
    public:
        const polynomial & tight_ineq () const { return m_tight_ineq; }
        polynomial& tight_ineq () { return m_tight_ineq; }
        const mpq & bound() const { return m_bound; }
        bool is_lower() const { return m_is_lower; }
        int expl_index() const { return m_expl_index; }
        literal() {}

        bool is_decided() const { return m_is_decided; }

        bool is_implied() const { return !m_is_decided;}

        unsigned var() const { return m_var; }
        static literal make_implied_literal(unsigned var_index, bool is_lower, const mpq & bound, int expl_index) {
            return literal(false, var_index, is_lower, bound, expl_index);
        }
        static literal make_decided_literal(unsigned var_index, bool is_lower, const mpq & bound, int expl_index) {
            return literal(true, var_index, is_lower, bound, expl_index);
        }
    };    



    
    enum class bound_type {
        LOWER, UPPER, UNDEF
    };
    struct bound_result {
        mpq m_bound;
        bound_type m_type;
        
        bound_result(const mpq & b, bound_type bt): m_bound(b), m_type(bt) {}
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
        const mpq & bound() const { return m_bound; }
    };
    
    class var_info {
        unsigned m_user_var_index;
        svector<unsigned> m_literals; // point to m_trail
        integer_domain<mpq> m_domain;
        std::unordered_set<int> m_dependent_constraints; // the set of constraintualities involving the var
    public:
        var_info(unsigned user_var_index) : m_user_var_index(user_var_index) {}
        const integer_domain<mpq> & domain() const { return m_domain; }
        unsigned user_var() const {
            return m_user_var_index;
        }
        void add_dependent_constraint(unsigned i) {
            m_dependent_constraints.insert(i);
        }
        void remove_depended_constraint(unsigned i) {
            m_dependent_constraints.erase(i);
        }
        bool intersect_with_lower_bound(const mpq & b) {

            return m_domain.intersect_with_lower_bound(b);
        }
        bool intersect_with_upper_bound(const mpq & b) {
            return m_domain.intersect_with_upper_bound(b);
        }
        bool is_fixed() const { return m_domain.is_fixed();}
        bool get_upper_bound(mpq & b) const { return m_domain.get_upper_bound(b); }
        bool get_lower_bound(mpq & b) const { return m_domain.get_lower_bound(b); }
        void print_var_domain(std::ostream &out) const { m_domain.print(out); }
        const svector<unsigned> & literals() const { return m_literals; }
        const std::unordered_set<int> & dependent_constraints() const { return m_dependent_constraints; }
        bool lower_bound_exists() const { return m_domain.lower_bound_exists();}
        bool upper_bound_exists() const { return m_domain.upper_bound_exists();}
        void push() {
            m_domain.push();
        }
        void pop() {
            m_domain.pop();
        }
        void add_literal(unsigned j) { m_literals.push_back(j); }
        void remove_literal(unsigned literal_index) {
            auto it = std::find (m_literals.begin(), m_literals.end(), literal_index);
            lp_assert(it != m_literals.end());
            m_literals.erase(it);
        }
        
    }; // end of var_info

    std::vector<var_info> m_var_infos;
    
    bool lhs_is_int(const std::vector<monomial> & lhs) const {
        for (auto & p : lhs) {
            if (numeric_traits<mpq>::is_int(p.coeff()) == false) return false;
        }
        return true;
    }

public:
    std::string get_column_name(unsigned j) const {
        return m_var_name_function(m_var_infos[j].user_var());
    }

    constraint & get_constraint(unsigned i) {
        return m_constraints[i];
    }
    
    const constraint & get_constraint(unsigned i) const {
        return m_constraints[i];
    }

    vector<constraint> m_constraints;
    std::vector<mpq> m_v; // the values of the variables
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
    static mpq m_local_zero;
    std::unordered_set<constraint_index> m_explanation; // if this collection is not empty we have a conflict 

    unsigned add_var(unsigned user_var_index) {
        unsigned ret;
        if (try_get_value(m_user_vars_to_cut_solver_vars, user_var_index, ret))
            return ret;
        unsigned j = m_var_infos.size();
        m_var_infos.push_back(var_info(user_var_index));
        return m_user_vars_to_cut_solver_vars[user_var_index] = j;      
    }

    bool is_lower_bound(const literal & l) const {
        return l.is_lower();
    }
    
    // bool lower_for_var(unsigned j, mpq & lower) const {
    //     bool ret = false;
    //     for (unsigned i : m_var_infos[j].m_literals)
    //         if (is_lower_bound(m_trail[i])) {
    //             if (ret == false) {
    //                 ret = true;
    //                 lower = get_bound(m_trail[i]);
    //             } else {
    //                 lower = std::max(lower, get_bound(m_trail[i]));
    //             }
    //         }
    //     return ret;
    // }

    bool is_upper_bound(const literal & l) const {
        return !l.is_lower();
    }
    
    // used for testing only
    void add_lower_bound_for_user_var(unsigned user_var_index, const mpq& bound) {
        unsigned j = m_user_vars_to_cut_solver_vars[user_var_index];
        m_var_infos[j].intersect_with_lower_bound(bound);
    }

    // used for testing only
    void add_upper_bound_for_user_var(unsigned user_var_index, const mpq& bound) {
        unsigned j = m_user_vars_to_cut_solver_vars[user_var_index];
        m_var_infos[j].intersect_with_upper_bound(bound);
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
        // mpqBD
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
        auto & d = m_var_infos[j];
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
        m_var_infos[l.var()].add_literal(literal_index);
    }

    unsigned find_large_enough_j(unsigned i) {
        unsigned r = 0;
        for (const auto & p : m_constraints[i].poly().m_coeffs) {
            r = std::max(r, p.var() + 1);
        }
        return r;
    }

    std::string var_name(unsigned j) const {
        return get_column_name(j);
    }

    
    void trace_print_domain_change(std::ostream& out,unsigned j, const mpq& v, const monomial & p, unsigned constraint_index) const {
        out << "trail.size() = " << m_trail.size() << "\n";
        out << "change in domain of " << var_name(j) << ", v = " << v << ", domain becomes ";
        print_var_domain(out, j);
        mpq lb;
        bool r = lower(constraint_index, lb);
        if (r)
            out << "lower_of_constraint = " << lb << "\n";
        else
            out << "no lower bound for constraint\n";
    }
    
    void propagate_monomial_on_lower(const monomial & p, const mpq& lower_val, unsigned constraint_index) {
        unsigned j = p.var();
        if (is_pos(p.coeff())) {
            mpq m;
            get_var_lower_bound(p.var(), m);
            mpq v = floor(- lower_val / p.coeff()) + m;
            bool change = m_var_infos[j].intersect_with_upper_bound(v);
            if (change) {
                TRACE("ba_int", trace_print_domain_change(tout, j, v, p, constraint_index););
                add_bound(v, j, false, constraint_index);
            }
        } else {
            mpq m;
            get_var_upper_bound(p.var(), m);
            mpq v = ceil( - lower_val / p.coeff()) + m;
            bool change = m_var_infos[j].intersect_with_lower_bound(v);
            if (change) {
                TRACE("ba_int", trace_print_domain_change(tout, j, v, p, constraint_index););
                add_bound(v, j, true , constraint_index);
            }
        }
    }

    void propagate_monomial_on_right_side(const monomial & p, const mpq& rs, unsigned constraint_index) {
        unsigned j = p.var();
        if (is_pos(p.coeff())) {
            mpq m;
            mpq v = floor(rs / p.coeff());
            bool change = m_var_infos[j].intersect_with_upper_bound(v);
            if (change) {
                TRACE("ba_int", trace_print_domain_change(tout, j, v, p, constraint_index););
                add_bound(v, j, false, constraint_index);
            }
        } else {
            mpq v = ceil(rs / p.coeff());
            bool change = m_var_infos[j].intersect_with_lower_bound(v);
            if (change) {
                TRACE("ba_int", trace_print_domain_change(tout, j, v, p, constraint_index););
                add_bound(v, j, true , constraint_index);
            }
        }
    }

    void print_var_info(std::ostream & out, const var_info & vi) const {
        out << m_var_name_function(vi.user_var()) << " ";
        print_var_domain(out, vi);
    }


    void print_var_info(std::ostream & out, unsigned j) const {
        print_var_info(out, m_var_infos[j]);
    }

    
    void print_var_domain(std::ostream & out, unsigned j) const {
        m_var_infos[j].print_var_domain(out);
    }

    void print_var_domain(std::ostream & out, const var_info & vi) const {
        vi.print_var_domain(out);
    }

    // b is the value of lower
    void propagate_constraint_on_lower(unsigned i, const mpq & b) {
        for (const auto & p: m_constraints[i].coeffs()) {
            propagate_monomial_on_lower(p, b, i);
        }
    }

    unsigned find_lower_bound_literal(bool is_lower, unsigned j, unsigned & upper_end_of_trail) const {
        TRACE("ba_int", tout << get_column_name(j) << "\n"; tout << "literal's size = " << m_var_infos[j].literals().size() << "\n";);
        for (unsigned k = m_var_infos[j].literals().size(); k--;) {
            unsigned literal_index = m_var_infos[j].literals()[k];
            if (literal_index >= upper_end_of_trail)
                continue;
            const literal& l = m_trail[literal_index];
            if (l.var() == j && l.is_lower() == is_lower) {
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
        for (auto j : m_constraints[constraint_index].origins())
            m_explanation.insert(j);
    }

    void fill_conflict_explanation(unsigned constraint_index, unsigned upper_end_of_trail) {
        // it is a depth search in the DAG of constraint: the chidlren of a constraint are those constraints that provide its lower bound
        add_constraint_origins_to_explanation(constraint_index);
        const constraint& in = m_constraints[constraint_index];
        TRACE("ba_int", print_constraint(tout, in););
        for (const auto & p: in.coeffs()){
            unsigned literal_index = find_lower_bound_literal(is_pos(p.coeff()), p.var(), upper_end_of_trail);
            unsigned l_constraint_index = m_trail[literal_index].expl_index();
            if (!m_constraints[l_constraint_index].is_simple()) 
                fill_conflict_explanation(l_constraint_index, literal_index);
            else
                add_constraint_origins_to_explanation(l_constraint_index);
        }
    }

    void trace_print_constraint(std::ostream& out, unsigned i) {
        print_constraint(out, i); out << "\n";
        unsigned j;
        auto pairs = to_pairs(m_constraints[i].poly().m_coeffs);
        auto it = linear_combination_iterator_on_std_vector<mpq>(pairs);
        while (it.next(j)) {
            out << "domain of " << var_name(j) << " = ";
            print_var_domain(out, j);
        }
    }

   
    void propagate_constraint_only_one_unlim(unsigned constraint_index, unsigned the_only_unlim) {
        const constraint& i = m_constraints[constraint_index];
        mpq rs = - i.poly().m_a;
        for (unsigned j = 0; j < i.poly().m_coeffs.size(); j++) {
            if (j == the_only_unlim) continue;
            mpq m;
            lower_monomial(i.poly().m_coeffs[j], m);
            rs -= m;
        }

        // we cannot get a conflict here because the monomial i.poly().m_coeffs[the_only_unlim]
        // is unlimited from below and we are adding an upper bound for it
        propagate_monomial_on_right_side(i.poly().m_coeffs[the_only_unlim], rs, constraint_index);
    }

    bool conflict() const { return m_explanation.size() > 0; }

    // returns -1 if there is no conflict and the index of the conflict constraint otherwise
    int  propagate_on_constraints_of_var(var_index j) {
        for (unsigned i : m_var_infos[j].dependent_constraints()) {
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

    bool get_var_lower_bound(var_index i, mpq & bound) const {
        const var_info & v = m_var_infos[i];
        return v.get_lower_bound(bound);
    }

    bool get_var_upper_bound(var_index i, mpq & bound) const {
        const var_info & v = m_var_infos[i];
        return v.get_upper_bound(bound);
    }

    bool lower_monomial_exists(const monomial & p) const {
        lp_assert(p.coeff() != 0);

        if (p.coeff() > 0) {
            if (!m_var_infos[p.var()].lower_bound_exists())
                return false;
        }
        else {
            if (!m_var_infos[p.var()].upper_bound_exists())
                return false;
        }
        return true;
    }

    bool upper_monomial_exists(const monomial & p) const {
        lp_assert(p.coeff() != 0);
        if (p.coeff() > 0) {
            if (!m_var_infos[p.var()].upper_bound_exists())
                return false;
        }
        else {
            if (!m_var_infos[p.var()].lower_bound_exists())
                return false;
        }
        return true;
    }

    
    // finds the lower bound of the monomial,
    // otherwise returns false
    bool lower_monomial(const monomial & p, mpq & lb) const {
        lp_assert(p.coeff() != 0);
        mpq var_bound;
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

    bool upper_monomial(const monomial & p, mpq & lb) const {
        lp_assert(p.coeff() != 0);
        mpq var_bound;
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
    bool lower(unsigned constraint_index, mpq & lb) const {
        return lower(m_constraints[constraint_index].poly(), lb);
    }

    // returns false if not limited from below
    // otherwise the answer is put into lb
    mpq lower_val(unsigned constraint_index) const {
        return lower_val(m_constraints[constraint_index]);
    }
    
    mpq lower_val(const constraint & i) const {
        mpq lb;
#if Z3DEBUG
        bool r =
#endif
            lower(i.poly(), lb);
        lp_assert(r);
        return lb;
    }

    void pop() { pop(1); }
        
    // returns false if not limited from below
    // otherwise the answer is put into lb
    bool lower(const polynomial & f, mpq & lb) const {
        lb = f.m_a;
        mpq lm;
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
        for (unsigned j = 0; j < f.poly().m_coeffs.size(); j++) {
            if (!lower_monomial_exists(f.poly().m_coeffs[j])) {
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
        mpq bound = p.m_a;
        for (const auto & t:  p.m_coeffs) {
            if (t.var() == j)
                continue;

            mpq l;
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
        mpq bound = p.m_a;
        for (const auto & t:  p.m_coeffs) {
            if (t.var() == j)
                continue;
            mpq b;
            upper_monomial(t, b);
            bound += b;
        }
        return bound_result(bound, bound_type::UPPER);
    }

    bool upper(const polynomial & p, mpq b) const {
        for (const auto & t:  p.m_coeffs) {
            if (!upper_monomial_exists(t)) {
                return false;
            }
        }
        // if we are here then there is an upper bound for p
        b = p.m_a;
        mpq bb;
        for (const auto & t:  p.m_coeffs) {
            upper_monomial(t, bb);
            b += bb;
        }
        return true;
    }

    
    
    
    // a is the coefficient before j
    bound_result bound_on_polynomial(const polynomial & p, const mpq& a, var_index j) const {
        lp_assert(!is_zero(a));
        if (numeric_traits<mpq>::is_pos(a)) {
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
        const mpq& a = q.poly().coeff(j);
        return bound_on_polynomial(q.poly(), a, j);
    }

    bound_result bound(unsigned constraint_index, var_index j) const {
        return bound(m_constraints[constraint_index], j);
    }

    
    void print_constraint(std::ostream & out, unsigned i)  const {
        print_constraint(out, m_constraints[i]);
    }

    void print_constraint(std::ostream & out, const constraint & i ) const {
        out << (i.is_lemma()? "lemma ": "assert ");
        if (!i.is_ineq()) {
            out << i.divider() << " | ";
        }
        print_polynomial(out, i.poly()) ;
        if (i.is_ineq()) {
            out << " <= 0";
        }
        out << "\n";
    }

    void print_literal(std::ostream & out, const literal & t) const {
        out << (t.is_decided()? "decided ": "implied ");
        out << get_column_name(t.var()) << " ";
        if (t.is_lower())
            out << ">= ";
        else
            out << "<= ";
        out << t.bound();

        if (t.is_decided() == false) {
            out << " by constraint " << t.expl_index() << " "; print_constraint(out, m_constraints[t.expl_index()]);
        }
        if (!t.tight_ineq().is_empty()) {
            out << " tight_ineq() = ";
            print_polynomial(out, t.tight_ineq());
            out << "\n";
        }
    }
    

    
    void print_polynomial(std::ostream & out, const polynomial & p) const {
        std::vector<std::pair<mpq, unsigned>> pairs = to_pairs(p.m_coeffs);
        this->print_linear_combination_of_column_indices_std(pairs, out);
        if (!is_zero(p.m_a)) {
            if (p.m_a < 0) {
                out << " - " << -p.m_a;
            } else {
                out << " + " << p.m_a;
            }
        }
    }

    // mpqrying to improve constraint "ie" by eliminating var j by using a  tight inequality 
    // for j. mpqhe left side of the inequality is passed as a parameter.
    void resolve(polynomial & ie, unsigned j, bool sign_j_in_ti_is_pos, const polynomial & ti) const {
        lp_assert(ti.is_tight(j));
        lp_assert(is_pos(ti.coeff(j)) == sign_j_in_ti_is_pos);
        auto &coeffs = ie.m_coeffs;
        // todo: implement a more efficient version
        bool found = false;
        mpq a;
        for (const auto & c : coeffs) {
            if (c.var() == j) {
                a = c.coeff();
                found = true;
                break;
            }
        }

        if (!found) {
            return;
        }

        if (is_neg(a)) {
            if (!sign_j_in_ti_is_pos)
                return;
            a = -a;
        } else {
            if (sign_j_in_ti_is_pos)
                return;
        }

        for (auto & c : ti.m_coeffs) {
            ie += monomial(a * c.coeff(), c.var());
        }

        ie.m_a += a * ti.m_a;
    }

    // returns true iff p imposes a better bound on j
    bool improves(var_index j, const constraint & p) const {
        auto a = p.coeff(j);
        if (is_zero(a))
            return false;
        const auto& dom = m_var_infos[j].domain();
        if (dom.is_empty())
            return false;
        if (is_pos(a)) {
            bound_result new_upper = bound(p, j);
            if (new_upper.m_type == bound_type::UNDEF)
                return false;
            mpq b;
            bool upper_bound_exists = get_var_upper_bound(j, b);
            return (!upper_bound_exists || new_upper.bound() < b) &&
                dom.intersection_with_upper_bound_is_empty(new_upper.bound());
        }

        lp_assert(is_neg(a));
        bound_result new_lower = bound(p, j);
        if (new_lower.m_type == bound_type::UNDEF)
            return false;
        mpq b;
        bool lower_bound_exists = get_var_lower_bound(j, b);
        return (!lower_bound_exists || new_lower.bound() > b) &&
            dom.intersection_with_lower_bound_is_empty(new_lower.bound());
    }


    // returns true iff br imposes a better bound on j
    bool improves(var_index j, const bound_result & br) const {
        if (br.m_type == bound_type::UNDEF)
            return false;
        const auto& dom = m_var_infos[j].domain();
        if (dom.is_empty())
            return false;
        if (br.m_type == bound_type::UPPER) {
            mpq b;
            bool j_has_upper_bound = get_var_upper_bound(j, b);
            return (!j_has_upper_bound || br.bound() < b) &&
                !dom.intersection_with_upper_bound_is_empty(br.bound());
        }

        if (br.m_type == bound_type::UNDEF)
            return false;
        mpq b;
        bool lower_bound_exists = get_var_lower_bound(j, b);
        return (!lower_bound_exists || br.bound() > b) &&
            !dom.intersection_with_lower_bound_is_empty(br.bound());
    }


    void add_bound(mpq v, unsigned j, bool is_lower, unsigned constraint_index) {
        push_literal_to_trail(literal::make_implied_literal(j, is_lower, v, constraint_index));
    }
    
    mpq value(const polynomial & p) const {
        mpq ret= p.m_a;
        for (const auto & t:p.m_coeffs)
            ret += t.coeff() * m_v[t.var()];
        return ret;
    }

    void pop_constraints() {
        for (unsigned j = m_constraints.size(); j-- > m_scope().m_constraints_size;) {
            const constraint & i = m_constraints[j];
            for (const auto & p: i.poly().m_coeffs) {
                m_var_infos[p.var()].remove_depended_constraint(j);
            }
            m_constraints.pop_back();
        }
        lp_assert(m_constraints.size() == m_scope().m_constraints_size);
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
            mpq b;
            lower(in.poly(), b);
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
        TRACE("decide_int", tout << "going to decide " << var_name(j) << " var index = " << j << "\n";
              tout << "domain = "; print_var_domain(tout, j); tout << "\n";);
        m_changed_vars.insert(j);
        decide_var_on_bound(j, flip_coin());
        return true;
    }

    lbool propagate_and_backjump_step() {
        do {
            int incostistent_constraint = propagate();
            TRACE("cs_dec", tout << "trail = \n"; print_trail(tout); tout << "end of trail\n";);

            if (incostistent_constraint >= 0) {
                if (at_base_lvl()) {
                    fill_conflict_explanation(incostistent_constraint, m_trail.size());
                    return lbool::l_false;
                }
                resolve_conflict(incostistent_constraint);
            }
        }
        while (m_changed_vars.size());
        
        if(!all_vars_are_fixed())
            return lbool::l_undef;

        return lbool::l_true;
    }

    bool decision_is_redundant_for_constraint(const polynomial& i, const literal & l) const {
        const mpq & coeff = i.coeff(l.var());
        if (is_zero(coeff))
            return true;
        return is_pos(coeff)? ! l.is_lower(): l.is_lower();
    }

    bool is_divizible(const mpq & a, const mpq & b) const {
        lp_assert(!is_zero(b));
        return is_zero(a % b);
    }
    
    void create_div_ndiv_parts_for_tightening(const polynomial & p, const mpq & coeff, polynomial & div_part, polynomial & ndiv_part) {
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

    void decide_lower_neg(polynomial &ndiv_part, const mpq & c, literal &l) {
        ndiv_part.add( -c, m_trail[l.expl_index()].tight_ineq());
        lp_assert(is_zero(ndiv_part.coeff(l.var())));
        TRACE("tight", tout << "Decided-Lower-Neg, "; tout << "ndiv_part = ";
              print_polynomial(tout, ndiv_part); tout << "\n";);
    }

    void decided_lower(const mpq & a, const mpq & c, polynomial &div_part, polynomial &ndiv_part, const literal &l) {
        mpq k = is_pos(a)?ceil( c / a): floor(c / a);
        ndiv_part += monomial(-c, l.var()); // it will be moved to div_part
        mpq a_k = a * k;
        mpq m = a_k - c;
        TRACE("tight", tout << "c = " << c << ", a = " << a <<
              ", c / a = " << c/a << ", k = " <<
              k << ", a * k = " << a * k << ", m = " << m << "\n"; );  
        lp_assert(!is_neg(m));
        create_tight_ineq_under_literal(l.expl_index());
        const literal & lex = m_trail[l.expl_index()];
        lp_assert(lex.var() == l.var());
        for (const monomial & n : lex.tight_ineq().m_coeffs) {
            if (n.var() == l.var()) {
                lp_assert(n.coeff() == one_of_type<mpq>());
                div_part += monomial(a_k, l.var());
            } else {
                ndiv_part += monomial(m * n.coeff(), n.var());
            }
        }
        ndiv_part += m * lex.tight_ineq().m_a;
        TRACE("tight", tout << "Decided-Lower ";
              tout << "div_part = ";
              print_polynomial(tout, div_part);
              tout << "\n";
              tout << "ndiv_part = ";
              print_polynomial(tout, ndiv_part);
              tout << "\n";);
    }

    void decided_upper(const mpq & a, const mpq & c, polynomial &div_part, polynomial &r, const literal &l) {
        // we would like to have c - ak > 0, or ak < c
        mpq k = is_pos(a)? floor( c / a): ceil(c / a);
        r += monomial(-c, l.var()); // it will be moved to div_part
        
        mpq a_k = a * k;
        mpq m = c - a_k;
        TRACE("tight", tout << "c = " << c << ", a = " << a <<
              ", c / a = " << c/a << ", k = " <<
              k << ", a * k = " << a * k << ", m = " << m << "\n"; );  
        lp_assert(!is_neg(m));
        create_tight_ineq_under_literal(l.expl_index());
        const literal & lex = m_trail[l.expl_index()];
        lp_assert(lex.var() == l.var());
        for (const monomial & n : lex.tight_ineq().m_coeffs) {
            if (n.var() == l.var()) {
                lp_assert(n.coeff() == -one_of_type<mpq>());
                div_part += monomial(a_k, l.var());
            } else {
                r += monomial(m * n.coeff(), n.var());
            }
        }
        r += m * lex.tight_ineq().m_a;
        TRACE("tight", tout << "Decided-Lower ";
              tout << "div_part = ";
              print_polynomial(tout, div_part);
              tout << "\n";
              tout << "r = ";
              print_polynomial(tout, r);
              tout << "\n";);
    }

    void tighten_on_literal(polynomial & p, const mpq & a, polynomial & div_part, polynomial &ndiv_part, int trail_index) {
        TRACE("tight",
              tout << "trail_index = " << trail_index << ", ";
              print_literal(tout, m_trail[trail_index]););
        literal & l = m_trail[trail_index];
        if (l.tight_ineq().number_of_monomials() == 0) {
            create_tight_ineq_under_literal(trail_index);
        }
        if (l.is_implied()) { // Resolve-Implied
            resolve(ndiv_part, l.var(), !l.is_lower(), l.tight_ineq());
            TRACE("tight", tout << "after resolve ndiv_part = "; print_polynomial(tout, ndiv_part);
                  tout << "\n";);
        } else { 
            lp_assert(l.is_decided());
            create_tight_ineq_under_literal(l.expl_index());
            mpq c = ndiv_part.coeff(l.var());
            if (l.is_lower()) {
                if (is_neg(c)) {
                    decide_lower_neg(ndiv_part, c, l);
                } else { 
                    decided_lower(a, c, div_part, ndiv_part, l);
                }
            } else {
                lp_assert(!l.is_lower());
                if (is_pos(c)) { // Decided-Upper-Pos
                    lp_assert(false);
                } else { 
                    decided_upper(a, c, div_part, ndiv_part, l);
                }
            }
        }
        
    }
    
    // see page 88 of Cutting to Chase
    void tighten(polynomial & p, unsigned j_of_var, const mpq& a, unsigned trail_index) {
        polynomial div_part, ndiv_part;
        ndiv_part.m_a = p.m_a;
        TRACE("tight",
              tout << "trail_index = " << trail_index;
              tout << ", var name = " << var_name(j_of_var) << ";";
              tout << "p = ";
              print_polynomial(tout, p););
        create_div_ndiv_parts_for_tightening(p, a, div_part, ndiv_part);
        int k = trail_index - 1;
        lp_assert(k >= 0);
        while (ndiv_part.number_of_monomials() > 0) {
            tighten_on_literal(p, a, div_part, ndiv_part, k--);
        }
        mpq abs_a = abs(a);
        p.clear();
        for (const auto & m : div_part.m_coeffs) {
            p.m_coeffs.push_back(monomial(m.coeff() / abs_a, m.var()));
        }
        p.m_a = ceil(ndiv_part.m_a / abs_a);
        lp_assert(p.m_a >= ndiv_part.m_a / abs_a);
        TRACE("tight", tout << "trail_index = " << trail_index << ", got tight p = "; print_polynomial(tout, p); tout << "\n";);
    }

    void create_tight_ineq_for_decided_literal(literal &l ) {
        if (!l.is_lower()) {
            l.tight_ineq().m_coeffs.push_back(monomial(l.var()));
            l.tight_ineq().m_a = -l.bound();
        } else {
            l.tight_ineq().m_coeffs.push_back(monomial(-one_of_type<mpq>(), l.var()));
            l.tight_ineq().m_a = l.bound();
        }
    }
    
    void create_tight_ineq_under_literal(unsigned trail_index) {
        literal & l = m_trail[trail_index];
        if (l.tight_ineq().number_of_monomials() > 0) {
            return;
        }
        if (l.is_decided()) {
            create_tight_ineq_for_decided_literal(l);
            return;
        }
        const constraint & c = m_constraints[l.expl_index()];
        lp_assert(c.is_ineq());
        polynomial &p = l.tight_ineq() = c.poly();
        unsigned j = l.var();
        const mpq& a = p.coeff(j);
        lp_assert(!is_zero(a));
        if (a == one_of_type<mpq>() || a == - one_of_type<mpq>()) {
            return;
        }
        tighten(p, j, a, trail_index);
    }

    bool backjump(polynomial &p,unsigned trail_index) {
        const literal &l = m_trail[trail_index];
        bound_result br = bound_on_polynomial(p,
                                              p.coeff(l.var()),
                                              l.var());
        TRACE("int_backjump", br.print(tout);
              tout << "; var info of " << get_column_name(l.var()) << ": ";
              print_var_info(tout, l.var());
              );
        unsigned var = l.var(); // l is goin away in the pop!
        while (m_trail.size() > trail_index)
            pop();
        lp_assert(m_trail.size() == trail_index);
        TRACE("int_backjump", tout << "var info after pop = ";  print_var_info(tout, var););
        add_lemma(p);
        add_bound(br.bound(), var, br.m_type == bound_type::LOWER, m_constraints.size() - 1);
        m_changed_vars.insert(var);
        restrict_var_domain_with_bound_result(var, br);
        TRACE("int_backjump", tout << "var info after restricton = ";
              print_var_info(tout, var);
              tout << "new literal = "; print_literal(tout, m_trail.back()););
        return true;  // we are done resolving
    }
    
    bool resolve_conflict_for_inequality_on_trail_element(polynomial & p, unsigned trail_index) {
        const literal & l = m_trail[trail_index];
        TRACE("int_resolve_confl", tout << "p = ";
              print_polynomial(tout, p);
              tout << "\nl = ";  print_literal(tout, l););
        if (l.is_decided()) {
            if (decision_is_redundant_for_constraint(p, l))
                pop(); // skip decision
            else {
                return backjump(p, trail_index);
            }
        } else { // the literal is implied
            create_tight_ineq_under_literal(trail_index);
            // applying Resolve rool
            resolve(p, l.var(), !l.is_lower(), l.tight_ineq());
            TRACE("int_resolve_confl",
                  tout << "p = ";
                  print_polynomial(tout, p););
        }
        return false;
    }
    
    bool resolve_conflict_for_inequality(const polynomial & confl_ineq) {
        polynomial p = confl_ineq;
#if Z3DEBUG
        mpq b;
        lp_assert(lower(p, b) && is_pos(b));
#endif
        bool done = false;
        unsigned j = m_trail.size() - 1;
        while (!done) {
            done = resolve_conflict_for_inequality_on_trail_element(p, j--);
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
        TRACE("int_resolve_confl", tout << "inconstistent_constraint = ";
              print_constraint(tout, i););
        if (i.is_ineq()) {
            return resolve_conflict_for_inequality(i.poly());
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

    void push_var_infos() {
        for (var_info & vi : m_var_infos)
            vi.push();
    }
    void pop_var_infos() {
        for (var_info & vi : m_var_infos)
            vi.pop();
    }
    
    void print_scope(std::ostream& out) const {
        out << "trail_size = " << m_scope().m_trail_size << ", constraints_size = " << m_scope().m_constraints_size << "\n";
    }

    void pop(unsigned k) {
        TRACE("trace_push_pop_in_cut_solver", tout << "before pop\n";print_state(tout););
        m_scope.pop(k);        
        TRACE("trace_push_pop_in_cut_solver", tout << "scope = ";print_scope(tout); tout << "\n";);
        pop_var_infos();

        for (unsigned j = m_trail.size(); j-- > m_scope().m_trail_size; ) {
            m_var_infos[m_trail[j].var()].remove_literal(j);
        }
            
        m_trail.resize(m_scope().m_trail_size);
        pop_constraints();
        m_decision_has_been_made.pop(k);
        TRACE("trace_push_pop_in_cut_solver", tout << "after pop\n";print_state(tout););
    }

    void push() {
        TRACE("trace_push_pop_in_cut_solver", print_state(tout););
        m_scope = scope(m_trail.size(), m_constraints.size());
        m_scope.push();
        m_decision_has_been_made.push();
        push_var_infos();
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
            if (!l.is_decided() && l.var() == j && l.is_lower() == is_lower)
                return k;
        }
        TRACE("find_literal", tout << "cannot find deciding literal for " << var_name(j)<<  " j = " << j << " is_lower = " << is_lower << std::endl;);
        lp_assert(false); // unreacheable
        return 0;
    }
    
    void decide_var_on_bound(unsigned j, bool decide_on_lower) {
        push();    
        mpq b;
        std::vector<monomial> lhs;
        if (decide_on_lower) {
            m_var_infos[j].domain().get_lower_bound(b);
            m_var_infos[j].intersect_with_upper_bound(b);
        }
        else {
            m_var_infos[j].domain().get_upper_bound(b);
            m_var_infos[j].intersect_with_lower_bound(b);
        }
        m_decision_has_been_made = true;
        push_literal_to_trail(literal::make_decided_literal(j, !decide_on_lower, b, find_literal_index(j, decide_on_lower)));
    }

    bool propagate_simple_constraint(unsigned constraint_index) {
        const constraint & t = m_constraints[constraint_index];
        TRACE("cut_solver_state_simpe_constraint",   print_constraint(tout, t); tout << std::endl;);
        var_index j = t.poly().m_coeffs[0].var();
        
        bound_result br = bound(t, j);
        TRACE("cut_solver_state_simpe_constraint", tout << "bound result = {"; br.print(tout); tout << "}\n";
              tout << "domain of " << get_column_name(j) << " = "; m_var_infos[j].domain().print(tout);
              tout << "\n";
              );
        
        if (improves(j, br)) {
            literal l = literal::make_implied_literal(j, br.m_type == bound_type::LOWER, br.bound(), constraint_index);
            l.tight_ineq() = t.poly();
            push_literal_to_trail(l);
            restrict_var_domain_with_bound_result(j, br);
            TRACE("cut_solver_state_simpe_constraint", tout <<"improved domain = ";
                  m_var_infos[j].domain().print(tout);
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
        // an option could be to check that upper(i.poly()) <= 0
        bool ret = value(i.poly()) <= zero_of_type<mpq>();
        if (!ret) {
            TRACE("cut_solver_state_inconsistent", 
                  tout << "inconsistent constraint "; print_constraint(tout,i); tout <<"\n";
                  tout << "value = " << value(i.poly()) << '\n';
                  );
        }
        return ret;
    }

    int find_non_fixed_var() const {
        // it is a very non efficient implementation for now.
        // the current limitation is that we only deal with bounded vars.
        // the search should be randomized.
        for (unsigned j = 0; j < m_var_infos.size(); j++) {
            const auto & d = m_var_infos[j].domain();
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


    void simplify_lemma(polynomial & p) {
        TRACE("simplify_lemma_int", tout << "p = "; print_polynomial(tout, p););
        auto & ms = p.m_coeffs;
        lp_assert(ms.size());
        mpq g;
        if (ms.size() == 1) {
            g = abs(ms[0].coeff());
        } else {
            g = gcd(ms[0].coeff(), ms[1].coeff());
            for (unsigned j = 2; j < ms.size(); j++) {
                g = gcd(g, ms[j].coeff());
            }
            lp_assert(is_pos(g));
        }
        if (g != one_of_type<mpq>()) {
            for (auto & m : ms)
                m.coeff() /= g;
            p.m_a = ceil(p.m_a /g);
        }
        TRACE("simplify_lemma_int", tout << "p = "; print_polynomial(tout, p););
    }
    
    void add_lemma(polynomial& p) {
        simplify_lemma(p);
        m_constraints.push_back(constraint::make_ineq_lemma(p));
    }
    
    unsigned add_ineq(const std::vector<monomial> & lhs,
                      const mpq& free_coeff,
                      svector<constraint_index> origins) {
        lp_assert(lhs_is_int(lhs));
        lp_assert(is_int(free_coeff));
        std::vector<monomial> local_lhs;
        unsigned constraint_index = m_constraints.size();
        for (auto & p : lhs)
            local_lhs.push_back(monomial(p.coeff(), add_var(p.var())));
        m_constraints.push_back(constraint::make_ineq_assert(local_lhs, free_coeff,origins)); 
        TRACE("ba_int",
              tout << "explanation :";
              for (auto i: origins) {
                  m_print_constraint_function(i, tout);
                  tout << "\n";
              });

        TRACE("add_ineq_int", tout << "m_constraints[" << m_constraints.size() - 1 << "] =  ";
              print_constraint(tout, m_constraints.back()); tout << "\n";);
        
        for (auto & p : local_lhs)
            m_var_infos[p.var()].add_dependent_constraint(constraint_index);
        
        return constraint_index;
    }
};

inline cut_solver::polynomial operator*(const mpq & a, cut_solver::polynomial & p) {
    cut_solver::polynomial ret;
    ret.m_a = p.m_a * a;
    
    for (const auto & t: p.m_coeffs)
        ret.m_coeffs.emplace_back(a * t.coeff(), t.var());
    
    return ret;
}

}

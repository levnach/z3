/*++
Copyright (c) 2016 Microsoft Corporation

Module Name:

    theory_lra.cpp

Abstract:

    <abstract>

Author:

    Lev Nachmanson (levnach) 2016-25-3
    Nikolaj Bjorner (nbjorner)

Revision History:


--*/
#include "util/lp/lp_solver.h"
#include "util/lp/lp_primal_simplex.h"
#include "util/lp/lp_dual_simplex.h"
#include "util/lp/indexed_value.h"
#include "util/lp/lar_solver.h"
#include "util/nat_set.h"
#include "util/optional.h"
#include "smt/smt_theory.h"
#include "smt/smt_context.h"
#include "smt/theory_lra.h"
#include "smt/proto_model/numeral_factory.h"
#include "smt/smt_model_generator.h"

namespace lp {
    enum bound_kind { lower_t, upper_t };

    //typedef std::pair<rational, rational> inf_numeral;
    typedef rational inf_numeral;

    std::ostream& operator<<(std::ostream& out, bound_kind const& k) {
        switch (k) {
        case lower_t: return out << "<=";
        case upper_t: return out << ">=";
        }
        return out;
    }

    class bound { 
    protected:
        smt::theory_var  m_var;
        inf_numeral      m_value; // tbd: inf_numeral
        bound_kind       m_bound_kind;
        bool             m_true;

    public:
        bound(smt::theory_var v, inf_numeral const & val, bound_kind k):
            m_var(v),
            m_value(val),
            m_bound_kind(k), 
            m_true(false) {
        }
        virtual ~bound() {}
        smt::theory_var get_var() const { return m_var; }
        bound_kind get_bound_kind() const { return m_true?m_bound_kind:static_cast<bound_kind>(1-m_bound_kind); }
        bool get_phase() const { return m_true; }
        void set_phase(bool t) { m_true = t; }
        inf_numeral const & get_value() const { return m_value; } // tbd: adjust with infinitesimals if negated.
        virtual std::ostream& display(std::ostream& out) const {
            return out << "v" << get_var() << "  " << get_bound_kind() << " " << get_value();
        }
    };

    std::ostream& operator<<(std::ostream& out, bound const& b) {
        return b.display(out);
    }

    struct stats {
        unsigned m_assert_lower;
        unsigned m_assert_upper;
        stats() { reset(); }
        void reset() {
            memset(this, 0, sizeof(*this));
        }
    };

    typedef optional<inf_numeral> opt_inf_numeral;

    struct replay_bound {
        smt::theory_var   m_v;
        opt_inf_numeral   m_bound;
        bound_kind        m_bound_kind;
        
        replay_bound(smt::theory_var v, opt_inf_numeral const& n, bound_kind k):
            m_v(v), m_bound(n), m_bound_kind(k)
        {}
    };


}

namespace smt {


    class theory_lra::imp {        

        struct scope {
            unsigned m_bounds_lim;
            unsigned m_asserted_qhead;            
            unsigned m_replay_lim;
            unsigned m_delayed_atom_lim;
        };

        struct delayed_atom {
            unsigned m_bv;
            bool     m_is_true;
            delayed_atom(unsigned b, bool t): m_bv(b), m_is_true(t) {}
        };


        theory_lra&         th;
        ast_manager&        m;
        arith_util          a;


        expr_ref_vector     m_terms;        // Internalization
        vector<rational>    m_coeffs;
        svector<theory_var> m_vars;
        unsigned            m_index;
        rational            m_coeff;
        vector<rational>    m_columns;

        svector<lean::var_index> m_theory_var2var_index;          // translate from theory variables to lar vars
        buffer<std::pair<rational, lean::var_index>>  m_left_side; // constraint left side
        std::unordered_map<lean::var_index, rational> m_variable_values; // current model
        svector<literal>                              m_inequalities;

        u_map<lp::bound*>      m_bool_var2bound;
        ptr_vector<lp::bound>  m_bounds;

        ptr_vector<lp::bound>  m_asserted_bounds;
        unsigned               m_asserted_qhead;

        vector<lp::replay_bound>   m_replay_bounds;

        vector<lp::opt_inf_numeral> m_lower, m_upper;

        bool                   m_delay_atoms;    // configuration
        svector<delayed_atom>  m_delayed_atoms;

        svector<scope>         m_scopes;
        lp::stats              m_stats;
        

        arith_factory*         m_factory;
        

        scoped_ptr<lean::lar_solver> m_solver;


        void found_non_arith(expr* n) {
            
        }

        bool decompose() {
            rational r;
            expr* n1, *n2;
            while (m_index < m_terms.size()) {
                expr* n = m_terms[m_index].get();
                if (a.is_add(n)) {
                    unsigned sz = to_app(n)->get_num_args();
                    for (unsigned i = 0; i < sz; ++i) {
                        m_terms.push_back(to_app(n)->get_arg(i));
                        m_coeffs.push_back(m_coeffs[m_index]);
                    }
                    m_terms[m_index] = m_terms.back();
                    m_coeffs[m_index] = m_coeffs.back();
                    m_terms.pop_back();
                    m_coeffs.pop_back();
                }
                else if (a.is_sub(n)) {
                    unsigned sz = to_app(n)->get_num_args();
                    m_terms[m_index] = to_app(n)->get_arg(0);                    
                    for (unsigned i = 1; i < sz; ++i) {
                        m_terms.push_back(to_app(n)->get_arg(i));
                        m_coeffs.push_back(-m_coeffs[m_index]);
                    }
                }
                else if (a.is_mul(n, n1, n2) && a.is_numeral(n1, r)) {
                    m_coeffs[m_index] *= r;
                    m_terms[m_index] = n2;
                }
                else if (a.is_mul(n, n1, n2) && a.is_numeral(n2, r)) {
                    m_coeffs[m_index] *= r;
                    m_terms[m_index] = n1;
                }
                else if (a.is_numeral(n, r)) {
                    m_coeff += r;
                    ++m_index;
                }
                else if (a.is_uminus(n, n1)) {
                    m_coeffs[m_index].neg();
                    m_terms[m_index] = n1;
                }
                else {
                    theory_var v = mk_var(n);
                    SASSERT(m_vars.size() == m_index);
                    m_vars.push_back(v);
                    ++m_index;
                }
            }
            return true;
        }

        theory_var mk_var(expr* n) {
            context& ctx = th.get_context();
            if (!ctx.e_internalized(n)) {
                ctx.internalize(n, false);
            }
            enode* e = ctx.get_enode(n);
            theory_var v;
            if (!th.is_attached_to_var(e)) {
                v = th.mk_var(e);                        
                ctx.attach_th_var(e, &th, v);
            }
            else {
                v = e->get_th_var(th.get_id());                
            }
            return v;
        }

        lean::var_index get_var_index(theory_var v) {
            lean::var_index result = UINT_MAX;
            if (m_theory_var2var_index.size() > static_cast<unsigned>(v)) {
                result = m_theory_var2var_index[v];
            }
            if (result == UINT_MAX) {
                std::ostringstream s;
                s << v;
                result = m_solver->add_var(s.str());
                m_theory_var2var_index.setx(v, result, UINT_MAX);
            }
            return result;
        }

        void init_left_side() {
            SASSERT(all_zeros(m_columns));
            for (unsigned i = 0; i < m_vars.size(); ++i) {
                theory_var column = m_vars[i];
                rational const& coeff = m_coeffs[i];
                if (m_columns.size() <= static_cast<unsigned>(column)) {
                    m_columns.setx(column, coeff, rational::zero());
                }
                else {
                    m_columns[column] += coeff;
                }                
            }
            m_left_side.reset();
            // reset the coefficients after they have been used.
            for (unsigned i = 0; i < m_vars.size(); ++i) {
                theory_var var = m_vars[i];
                rational const& r = m_columns[var];
                if (!r.is_zero()) {
                    m_left_side.push_back(std::make_pair(r, get_var_index(var)));
                    m_columns[var].reset();                    
                }
            }
            SASSERT(all_zeros(m_columns));
        }

        bool all_zeros(vector<rational> const& v) const {
            for (unsigned i = 0; i < v.size(); ++i) {
                if (!v[i].is_zero()) {
                    return false;
                }
            }
            return true;
        }

        void internalize_ineq(expr* atom, bool_var bv, bool is_true) {
            TRACE("arith", tout << mk_pp(atom, m) << " " << is_true << "\n";);
            expr* n1, *n2;
            rational right_side;
            lean::lconstraint_kind k = lean::EQ;
            reset_term();
            m_coeffs.push_back(rational(1));
            if (a.is_le(atom, n1, n2) && a.is_numeral(n2, right_side)) {
                m_terms.push_back(n1);
                k = is_true ? lean::LE : lean::GT;
            }
            else if (a.is_ge(atom, n1, n2) && a.is_numeral(n2, right_side)) {
                k = is_true ? lean::GE : lean::LT;
            }
            else {
                UNREACHABLE();
            }
            VERIFY(decompose());

            init_left_side();
            
            lean::constraint_index index = m_solver->add_constraint(m_left_side, k, right_side);
            
            m_inequalities.setx(index, literal(bv, !is_true), null_literal);
        }

        void internalize() {
            unsigned row = 0; // TBD m_solver->row_count();

            init_left_side();

            
            // TBD: m_solver->add_constraint(lean::Equal, rational::zero(), row);

        }

        void reset_term() {
            m_terms.reset();
            m_coeff.reset();
            m_coeffs.reset();
            m_vars.reset();
            m_index = 0;
        }

        void del_bounds(unsigned old_size) {
            for (unsigned i = m_bounds.size(); i > old_size; ) {
                --i;
                dealloc(m_bounds[i]);
            }
            m_bounds.shrink(old_size);
        }

        typedef lean::lp_dual_simplex<lp::inf_numeral, rational> dual_simplex;

    public:
        imp(theory_lra& th, ast_manager& m): 
            th(th), m(m), a(m), m_terms(m), 
            m_asserted_qhead(0), m_delay_atoms(true) {
            m_solver = alloc(lean::lar_solver); 
        }

        ~imp() {
            del_bounds(0);
        }

        bool internalize_atom(app * atom, bool gate_ctx) {
            context& ctx = th.get_context();
            bool_var bv = ctx.mk_bool_var(atom);
            ctx.set_var_theory(bv, th.get_id());
            if (m_delay_atoms) {
                return true;
            }
            expr* n1, *n2;
            rational r;
            lp::bound_kind k;
            theory_var v = null_theory_var;
            if (a.is_le(atom, n1, n2) && a.is_numeral(n2, r)) {
                v = internalize_term_core(n1);
                k = lp::upper_t;
            }
            else if (a.is_ge(atom, n1, n2) && a.is_numeral(n2, r)) {
                v = internalize_term_core(n1);
                k = lp::lower_t;
            }    
            else {
                TRACE("arith", tout << "Could not internalize " << mk_pp(atom, m) << "\n";);
                return false;
            }
            lp::inf_numeral _r(r);
            lp::bound* b = alloc(lp::bound, v, _r, k);
            m_bounds.push_back(b);
            m_bool_var2bound.insert(bv, b);
            TRACE("arith", tout << "Internalized " << mk_pp(atom, m) << " " << b << "\n";);
            return true;
        }

        theory_var internalize_term_core(expr* term) {
            reset_term();
            m_terms.push_back(term);
            m_coeffs.push_back(rational(1));
            if (!decompose()) {
                return null_theory_var;
            }
            switch (m_terms.size()) {
            case 0:
                break;
            case 1:
                // TBD optimize for single variable case
                if (m_coeff.is_zero() && m_coeffs[0].is_one()) {
                    return m_vars[0];
                }
                break;
            default:
                break;
            }
            theory_var v = mk_var(term);
            if (v == null_theory_var) {
                return v;
            }
            m_vars.push_back(v);
            m_coeffs.push_back(rational::minus_one());
            internalize();
            return v;
        }

        bool internalize_term(app * term) {
            return null_theory_var != internalize_term_core(term);
        }

        void internalize_eq_eh(app * atom, bool_var v) {

            
        }
        void assign_eh(bool_var v, bool is_true) {
            if (m_delay_atoms) {
                m_delayed_atoms.push_back(delayed_atom(v, is_true));
            }
            else {
                lp::bound* b = m_bool_var2bound[v];
                SASSERT(b);
                b->set_phase(is_true);
                m_asserted_bounds.push_back(b);
            }
        }

        void new_eq_eh(theory_var v1, theory_var v2) {

        }
        bool use_diseqs() const {
            return false;
        }
        void new_diseq_eh(theory_var v1, theory_var v2) {
            UNREACHABLE();
        }

        void push_scope_eh() {
            m_scopes.push_back(scope());
            scope& s = m_scopes.back();
            s.m_bounds_lim = m_bounds.size();
            s.m_asserted_qhead = m_asserted_qhead;
            s.m_replay_lim = m_replay_bounds.size();
            s.m_delayed_atom_lim = m_delayed_atoms.size();
        }

        void pop_scope_eh(unsigned num_scopes) {
            if (num_scopes > 0) {
                unsigned old_size = m_scopes.size() - num_scopes;
                del_bounds(m_scopes[old_size].m_bounds_lim);
                m_delayed_atoms.shrink(m_scopes[old_size].m_delayed_atom_lim);
                m_asserted_qhead = m_scopes[old_size].m_asserted_qhead;
                unsigned sz = m_scopes[old_size].m_replay_lim;
                for (unsigned i = m_replay_bounds.size(); i > sz; ) {
                    --i;
                    lp::replay_bound const& r = m_replay_bounds[i];
                    switch (r.m_bound_kind) {
                    case lp::lower_t:
                        m_lower[r.m_v] = r.m_bound;
                        // TBD: m_solver->unset_low_bound(r.m_v);
                        if (r.m_bound) {
                            // TBD: m_solver->set_low_bound(r.m_v, *r.m_bound);
                        }
                        break;
                    case lp::upper_t:
                        m_upper[r.m_v] = r.m_bound;
                        // TBD: m_solver->unset_upper_bound(r.m_v);
                        if (r.m_bound) {
                            // TBD: m_solver->set_upper_bound(r.m_v, *r.m_bound);
                        }
                        break;
                    }
                }
                m_replay_bounds.shrink(sz);
                m_scopes.resize(old_size);
            }
        }

        void restart_eh() {

        }
        void relevant_eh(app* e) {

        }
        void init_search_eh() {

        }
        final_check_status final_check_eh() {
            if (m_delayed_atoms.empty()) {
                return FC_DONE;
            }
            m_solver = alloc(lean::lar_solver); 
            context& ctx = th.get_context();
            for (unsigned i = 0; i < m_delayed_atoms.size(); ++i) {
                bool_var bv = m_delayed_atoms[i].m_bv;
                expr* atom = ctx.bool_var2expr(bv);
                internalize_ineq(atom, bv, m_delayed_atoms[i].m_is_true);
            }
            lbool is_sat = make_feasible();
            switch (is_sat) {
            case l_true:
                return FC_DONE;
            case l_false:
                failed();
                return FC_CONTINUE;
            case l_undef:
                return FC_GIVEUP;
            default:
                UNREACHABLE();
                break;
            }
            return FC_GIVEUP;

        }

        bool is_shared(theory_var v) const {
            return false;
        }
        bool can_propagate() {
            return m_asserted_bounds.size() > m_asserted_qhead;
        }
        void propagate() {
            while (m_asserted_qhead < m_asserted_bounds.size()) {
                lp::bound* b = m_asserted_bounds[m_asserted_qhead];
                if (!assert_bound(*b)) {
                    failed();
                    return;
                }
                ++m_asserted_qhead;
            }
            switch (make_feasible()) {
            case l_false:
                failed();
                return;
            case l_true:
                break;
            case l_undef:
                // canceled
                break;
            }
        }

        bool assert_bound(lp::bound& b) {
            lp::bound_kind k = b.get_bound_kind();
            theory_var v = b.get_var();
            lp::inf_numeral const& val = b.get_value();
            switch (k) {
            case lp::lower_t:
                m_replay_bounds.push_back(lp::replay_bound(v, m_lower[v], k));
                // TBD: m_solver->set_low_bound(v, val);
                m_lower[v] = val;
                ++m_stats.m_assert_lower;
                break;
            case lp::upper_t:
                m_replay_bounds.push_back(lp::replay_bound(v, m_upper[v], k));
                // TBD: m_solver->set_upper_bound(v, val);
                m_upper[v] = val;
                ++m_stats.m_assert_upper;
                break;
            }
            return true;
        }

        lbool make_feasible() {
            lean::lp_status status = m_solver->check();
            switch (status) {
            case lean::lp_status::INFEASIBLE:
                return l_false;
            case lean::lp_status::FEASIBLE:
            case lean::lp_status::OPTIMAL:
                return l_true;
            default:
                TRACE("arith", tout << "status treated as inconclusive: " << status << "\n";);
                // TENTATIVE_UNBOUNDED, UNBOUNDED, TENTATIVE_DUAL_UNBOUNDED, DUAL_UNBOUNDED, 
                // FLOATING_POINT_ERROR, TIME_EXAUSTED, ITERATIONS_EXHAUSTED, EMPTY, UNSTABLE
                return l_undef;
            }
        }

        void failed() {
            buffer<std::pair<rational, lean::constraint_index>> evidence;
            literal_vector core;
            m_solver->get_infeasibility_evidence(evidence);
            for (unsigned i = 0; i < evidence.size(); ++i) {
                if (!evidence[i].first.is_zero()) {
                    core.push_back(m_inequalities[evidence[i].second]);
                }
            }
            context& ctx = th.get_context();
            ctx.set_conflict(
                ctx.mk_justification(
                    theory_conflict_justification(th.get_id(), ctx.get_region(), core.size(), core.c_ptr(), 0, 0)));
        }

        justification * why_is_diseq(theory_var v1, theory_var v2) {
            return 0;
        }
        void reset_eh() {

        }
        void init_model(model_generator & mg) {
            m_variable_values.clear();
            m_solver->get_model(m_variable_values);
            m_factory = alloc(arith_factory, m);
            mg.register_factory(m_factory);
        }

        model_value_proc * mk_value(enode * n, model_generator & mg) {
            bool is_int = false;
            theory_var v = n->get_th_var(th.get_id());
            SASSERT(v != null_theory_var);
            SASSERT(v < static_cast<theory_var>(m_theory_var2var_index.size()));
            SASSERT(UINT_MAX != m_theory_var2var_index[v]);
            rational num = m_variable_values[m_theory_var2var_index[v]];
            return alloc(expr_wrapper_proc, m_factory->mk_value(num, is_int));
        }

        bool validate_eq_in_model(theory_var v1, theory_var v2, bool is_true) const {
            return false;
        }
        void display(std::ostream & out) const {
            if (m_solver) {
                m_solver->print_constraints(out);
            }
        }
        void collect_statistics(::statistics & st) const {

        }        
    };
    
    theory_lra::theory_lra(ast_manager& m):
        theory(m.get_family_id("arith")) {
        m_imp = alloc(imp, *this, m);
    }    
    theory_lra::~theory_lra() {
        dealloc(m_imp);
    }   
    theory* theory_lra::mk_fresh(context* new_ctx) {
        return alloc(theory_lra, new_ctx->get_manager());
    }
    void theory_lra::init(context * ctx) {
        theory::init(ctx);
    }
    bool theory_lra::internalize_atom(app * atom, bool gate_ctx) {
        return m_imp->internalize_atom(atom, gate_ctx);
    }
    bool theory_lra::internalize_term(app * term) {
        return m_imp->internalize_term(term);
    }
    void theory_lra::internalize_eq_eh(app * atom, bool_var v) {
        m_imp->internalize_eq_eh(atom, v);
    }
    void theory_lra::assign_eh(bool_var v, bool is_true) {
        m_imp->assign_eh(v, is_true);
    }
    void theory_lra::new_eq_eh(theory_var v1, theory_var v2) {
        m_imp->new_eq_eh(v1, v2);
    }
    bool theory_lra::use_diseqs() const {
        return m_imp->use_diseqs();
    }
    void theory_lra::new_diseq_eh(theory_var v1, theory_var v2) {
        m_imp->new_diseq_eh(v1, v2);
    }
    void theory_lra::push_scope_eh() {
        theory::push_scope_eh();
        m_imp->push_scope_eh();
    }
    void theory_lra::pop_scope_eh(unsigned num_scopes) {
        m_imp->pop_scope_eh(num_scopes);
        theory::pop_scope_eh(num_scopes);
    }
    void theory_lra::restart_eh() {
        m_imp->restart_eh();
    }
    void theory_lra::relevant_eh(app* e) {
        m_imp->relevant_eh(e);
    }
    void theory_lra::init_search_eh() {
        m_imp->init_search_eh();
    }
    final_check_status theory_lra::final_check_eh() {
        return m_imp->final_check_eh();
    }
    bool theory_lra::is_shared(theory_var v) const {
        return m_imp->is_shared(v);
    }
    bool theory_lra::can_propagate() {
        return m_imp->can_propagate();
    }
    void theory_lra::propagate() {
        m_imp->propagate();
    }
    justification * theory_lra::why_is_diseq(theory_var v1, theory_var v2) {
        return m_imp->why_is_diseq(v1, v2);
    }
    void theory_lra::reset_eh() {
        m_imp->reset_eh();
    }
    void theory_lra::init_model(model_generator & m) {
        m_imp->init_model(m);
    }
    model_value_proc * theory_lra::mk_value(enode * n, model_generator & mg) {
        return m_imp->mk_value(n, mg);
    }
    bool theory_lra::validate_eq_in_model(theory_var v1, theory_var v2, bool is_true) const {
        return m_imp->validate_eq_in_model(v1, v2, is_true);
    }
    void theory_lra::display(std::ostream & out) const {
        m_imp->display(out);
    }
    void theory_lra::collect_statistics(::statistics & st) const {
        m_imp->collect_statistics(st);
    }

}

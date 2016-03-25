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
#include "smt/smt_theory.h"
#include "smt/theory_lra.h"
#include "util/lp/lp_solver.h"

namespace smt {
    class theory_lra::imp {

    };

#if 0
        theory_lra(ast_manager& m);
        ~theory_lra();
        theory* mk_fresh(context* new_ctx);        
        void init(context * ctx);
        bool internalize_atom(app * atom, bool gate_ctx);                                                     
        bool internalize_term(app * term);
        void internalize_eq_eh(app * atom, bool_var v);
        void assign_eh(bool_var v, bool is_true);
        void new_eq_eh(theory_var v1, theory_var v2);
        bool use_diseqs() const;
        void new_diseq_eh(theory_var v1, theory_var v2);
        void push_scope_eh();
        void pop_scope_eh(unsigned num_scopes);
        void restart_eh();
        void relevant_eh(app* e);
        void init_search_eh();
        final_check_status final_check_eh();
        bool is_shared(theory_var v) const;
        bool can_propagate();
        void propagate();
        justification * why_is_diseq(theory_var v1, theory_var v2);
        // void flush_eh();
        void reset_eh();
        void init_model(model_generator & m);
        model_value_proc * mk_value(enode * n, model_generator & mg);
        bool validate_eq_in_model(theory_var v1, theory_var v2, bool is_true) const;
        void display(std::ostream & out) const;
        void collect_statistics(::statistics & st) const;
#endif

}

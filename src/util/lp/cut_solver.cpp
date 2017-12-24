/*
  Copyright (c) 2017 Microsoft Corporation
  Author: Nikolaj Bjorner, Lev Nachmanson
*/
#include "util/lp/cut_solver.h"
namespace lp {
    mpq polynomial::m_local_zero = zero_of_type<mpq>();

    size_t cut_solver::ccns_hash::operator() (const cut_solver::constraint* c) const { return c->id(); }
    
    bool cut_solver::ccns_equal::operator() (const cut_solver::constraint* a, const cut_solver::constraint * b) const { return a->id() == b->id(); }

    std::ostream& operator<<(std::ostream& out, pp_poly const& p) {
        p.s.print_polynomial(out, p.p);
        return out;
    }

}


/*++
  Copyright (c) 2017 Microsoft Corporation

  Module Name:

  <name>

  Abstract:

  <abstract>

  Author:
  Nikolaj Bjorner (nbjorner)
  Lev Nachmanson (levnach)

  Revision History:


  --*/
#pragma once
#include "util/lp/lp_types.h"
#include "util/lp/column_info.h"
#include "util/lp/explanation.h"

namespace nla {
typedef lp::constraint_index lpci;
typedef lp::lconstraint_kind llc;
typedef lp::constraint_index     lpci;
typedef lp::explanation          expl_set;
typedef lp::var_index            lpvar;

struct from_index_dummy{};

class signed_var {
    unsigned m_sv;
public:
    // constructor, sign = true means minus
    signed_var(lpvar v, bool sign): m_sv((v << 1) + (sign ? 1 : 0)) {}
    // constructor
    signed_var(unsigned sv, from_index_dummy): m_sv(sv) {}
    bool sign() const { return 0 != (m_sv & 0x1); }
    lpvar var() const { return m_sv >> 1; }
    unsigned index() const { return m_sv; }        
    void neg() { m_sv = m_sv ^ 1; }
    friend signed_var operator~(signed_var const& sv) {
        return signed_var(sv.var(), !sv.sign());
    }
    bool operator==(signed_var const& other) const {
        return m_sv == other.m_sv;
    }
    bool operator!=(signed_var const& other) const {
        return m_sv != other.m_sv;
    }
    rational rsign() const { return sign() ? rational::minus_one() : rational::one(); }

    std::ostream& display(std::ostream& out) const {
        return out << (sign()?"-":"") << var();
    }
};

inline std::ostream& operator<<(std::ostream& out, signed_var const& sv) { return sv.display(out); }

/*
 *  represents definition m_v = coeff* v1*v2*...*vn, 
 *  where m_vs = [v1, v2, .., vn]
 */
class monomial_coeff  {
    svector<lp::var_index> m_vs;
    rational m_coeff;
public:
    monomial_coeff(const svector<lp::var_index>& vs, rational const& coeff): m_vs(vs), m_coeff(coeff) {}
    rational const& coeff() const { return m_coeff; }
    const svector<lp::var_index> & vars() const { return m_vs; } 
};
}

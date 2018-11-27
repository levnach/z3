
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
#include "util/rational.h"
#include "util/lp/monomial.h"

namespace nla {

struct factorization_factory;
typedef unsigned lpvar;

class factorization {
    svector<lpvar>         m_vars;
public:
    bool is_empty() const { return m_vars.empty(); }
    svector<lpvar> & vars() { return m_vars; }
    const svector<lpvar> & vars() const { return m_vars; }
    unsigned operator[](unsigned k) const { return m_vars[k]; }
    size_t size() const { return m_vars.size(); }
    const lpvar* begin() const { return m_vars.begin(); }
    const lpvar* end() const { return m_vars.end(); }
};

struct const_iterator_mon {
    // fields
    svector<bool>                  m_mask;
    const factorization_factory *  m_ff;
    bool                           m_full_factorization_returned;

    //typedefs
    typedef const_iterator_mon self_type;
    typedef factorization value_type;
    typedef const factorization reference;
    typedef int difference_type;
    typedef std::forward_iterator_tag iterator_category;

    void init_vars_by_the_mask(unsigned_vector & k_vars, unsigned_vector & j_vars) const;
            
    bool get_factors(unsigned& k, unsigned& j, rational& sign) const;

    reference operator*() const;
    void advance_mask();
            
    self_type operator++();
    self_type operator++(int);

    const_iterator_mon(const svector<bool>& mask, const factorization_factory *f);
    
    bool operator==(const self_type &other) const;
    bool operator!=(const self_type &other) const;
            
    factorization create_binary_factorization(lpvar j, lpvar k) const;
    
    factorization create_full_factorization() const;
};

struct factorization_factory {
    const svector<lpvar>& m_vars;
    // returns true if found
    virtual bool find_monomial_of_vars(const svector<lpvar>& vars, monomial& m, rational & sign) const = 0;

    factorization_factory(const svector<lpvar>& vars) :
        m_vars(vars) {
    }

    const_iterator_mon begin() const {
        // we keep the last element always in the first factor to avoid
        // repeating a pair twice
        svector<bool> mask(m_vars.size() - 1, false);
        return const_iterator_mon(mask, this);
    }
    
    const_iterator_mon end() const {
        svector<bool> mask(m_vars.size() - 1, true);
        auto it = const_iterator_mon(mask, this);
        it.m_full_factorization_returned = true;
        return it;
    } 
   
};

}

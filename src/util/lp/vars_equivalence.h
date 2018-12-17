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

namespace nla {

typedef lp::constraint_index     lpci;
typedef lp::explanation          expl_set;
typedef lp::var_index            lpvar;
struct hash_svector {
    size_t operator()(const unsigned_vector & v) const {
        return svector_hash<unsigned_hash>()(v);
    }
};

struct index_with_sign {
    unsigned m_i; // the index
    rational m_sign; // the sign: -1 or 1
    index_with_sign(unsigned i, rational sign) : m_i(i), m_sign(sign) {}
    index_with_sign() {}
    bool operator==(const index_with_sign& b) {
        return m_i == b.m_i && m_sign == b.m_sign;
    }
    unsigned var() const { return m_i; }
    const rational& sign() const { return m_sign; }
};

struct rat_hash {
    typedef rational data;
    unsigned operator()(const rational& x) const { return x.hash(); }
};


struct hash_vector {
    size_t operator()(const vector<rational> & v) const {
        return vector_hash<rat_hash>()(v);
    }
};

struct vars_equivalence {
    
    std::unordered_map<rational, unsigned_vector> m_vars_by_abs_values;
    std::unordered_map<vector<rational>,
                       unsigned_vector,
                       hash_vector>               m_monomials_by_abs_vals;

    std::function<rational(lpvar)>                m_vvr;


    // constructor
    vars_equivalence(std::function<rational(lpvar)> vvr) : m_vvr(vvr) {}
    
    const std::unordered_map<vector<rational>,
                             unsigned_vector,
                             hash_vector>& monomials_by_abs_values() const {
        return m_monomials_by_abs_vals;
    }

    
    void clear() {
        m_vars_by_abs_values.clear();
        m_monomials_by_abs_vals.clear();
    }

    svector<lpvar> get_vars_with_the_same_abs_val(const rational& v) const {
        svector<unsigned> ret;
        auto it = m_vars_by_abs_values.find(abs(v));
        if (it == m_vars_by_abs_values.end())
            return ret;

        return it->second; 
    } 

    
    bool empty() const {
        return m_vars_by_abs_values.empty();
    }

    bool is_root(unsigned j) const {
        auto it = m_vars_by_abs_values.find(abs(m_vvr(j)));
        SASSERT(it != m_vars_by_abs_values.end());
        return *(it->second.begin()) == j;
    }

    
    // Finds the root var which is equivalent to j.
    // The sign is flipped if needed
    lpvar map_to_root(lpvar j, rational& sign) const {
        auto v = m_vvr(j);
        sign = v.is_pos()? rational(1) : -rational(1);
        auto it = m_vars_by_abs_values.find(abs(v));
        SASSERT(it != m_vars_by_abs_values.end());
        return *(it->second.begin());
    }

    // Finds the root var which is equivalent to j.
    lpvar map_to_root(lpvar j) const {
        auto it = m_vars_by_abs_values.find(abs(m_vvr(j)));
        SASSERT(it != m_vars_by_abs_values.end());
        return *(it->second.begin());
    }

    void register_var(unsigned j, const rational& val) {
        TRACE("nla_vars_eq", tout << "j = " << j;);
        rational v = abs(val);
        auto it = m_vars_by_abs_values.find(v);
        if (it == m_vars_by_abs_values.end()) {
            unsigned_vector uv;
            uv.push_back(j);
            m_vars_by_abs_values[v] = uv;
        } else {
            it->second.push_back(j);
        }
    }

    void deregister_monomial_from_abs_vals(const monomial & m, unsigned i){
        int sign;
        auto key = get_sorted_abs_vals_from_mon(m, sign);
        SASSERT(m_monomials_by_abs_vals.find(key)->second.back() == i);
        m_monomials_by_abs_vals.find(key)->second.pop_back();
    }

    vector<rational> get_sorted_abs_vals_from_mon(const monomial& m, int & sign) {
        sign = 1;
        vector<rational> abs_vals;
        for (lpvar j : m) {
            const rational v = m_vvr(j);
            abs_vals.push_back(abs(v));
            if (v.is_neg()) {
                sign = -sign;
            }
        }
        std::sort(abs_vals.begin(), abs_vals.end());
        return abs_vals;
    }
    
    void register_monomial_in_abs_vals(unsigned i, const monomial & m ) {
        int sign;
        vector<rational> abs_vals = get_sorted_abs_vals_from_mon(m, sign);
        auto it = m_monomials_by_abs_vals.find(abs_vals);
        if (it == m_monomials_by_abs_vals.end()) {
            unsigned_vector v;
            v.push_back(i);
            // v is a vector containing a single index_with_sign
            m_monomials_by_abs_vals.emplace(abs_vals, v);
        } 
        else {
            it->second.push_back(i);
        }
        
    }
}; // end of vars_equivalence
}

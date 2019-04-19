#include "util/vector.h"
#include "util/lp/factorization.h"
namespace nla {

void const_iterator_mon::init_vars_by_the_mask(unsigned_vector & k_vars, unsigned_vector & j_vars) const {
    // the last element for m_factorization.m_rooted_vars goes to k_vars
    SASSERT(m_mask.size() + 1  == m_ff->m_vars.size());
    k_vars.push_back(m_ff->m_vars.back()); 
    for (unsigned j = 0; j < m_mask.size(); j++) {
        if (m_mask[j]) {
            k_vars.push_back(m_ff->m_vars[j]);
        } else {
            j_vars.push_back(m_ff->m_vars[j]);
        }
    }
}
            
bool const_iterator_mon::get_factors(factor& k, factor& j, rational& sign) const {
    unsigned_vector k_vars;
    unsigned_vector j_vars;
    init_vars_by_the_mask(k_vars, j_vars);
    SASSERT(!k_vars.empty() && !j_vars.empty());
    std::sort(k_vars.begin(), k_vars.end());
    std::sort(j_vars.begin(), j_vars.end());

    if (k_vars.size() == 1) {
        k.set(k_vars[0], factor_type::VAR);
    } else {
        unsigned i;
        if (!m_ff->find_rm_monomial_of_vars(k_vars, i)) {
            return false;
        }
        k.set(i, factor_type::RM);
    }

    if (j_vars.size() == 1) {
        j.set(j_vars[0], factor_type::VAR);
    } else {
        unsigned i;
        if (!m_ff->find_rm_monomial_of_vars(j_vars, i)) {
            return false;
        }
        j.set(i, factor_type::RM);
    }
    return true;
}

factorization const_iterator_mon::operator*() const {
    if (m_full_factorization_returned == false)  {
        return create_full_factorization(m_ff->m_monomial);
    }
    factor j, k; rational sign;
    if (!get_factors(j, k, sign))
        return factorization(nullptr);
    return create_binary_factorization(j, k);
}
            
void const_iterator_mon::advance_mask() {
    if (!m_full_factorization_returned) {
        m_full_factorization_returned = true;
        return;
    }
    for (bool& m : m_mask) {
        if (m) {
            m = false;
        }
        else {
            m = true;
            break;
        }
    }
}

            
const_iterator_mon::self_type const_iterator_mon::operator++() {  self_type i = *this; operator++(1); return i;  }
const_iterator_mon::self_type const_iterator_mon::operator++(int) { advance_mask(); return *this; }

const_iterator_mon::const_iterator_mon(const svector<bool>& mask, const factorization_factory *f) : 
    m_mask(mask),
    m_ff(f) ,
    m_full_factorization_returned(false)
{}
            
bool const_iterator_mon::operator==(const const_iterator_mon::self_type &other) const {
    return
        m_full_factorization_returned == other.m_full_factorization_returned &&
        m_mask == other.m_mask;
}

bool const_iterator_mon::operator!=(const const_iterator_mon::self_type &other) const { return !(*this == other); }
            
factorization const_iterator_mon::create_binary_factorization(factor j, factor k) const {
    factorization f(nullptr);
    f.push_back(j);
    f.push_back(k);  
    return f;
}

factorization const_iterator_mon::create_full_factorization(const monomial* m) const {
    if (m != nullptr)
        return factorization(m);
    factorization f(nullptr);
    for (lpvar j : m_ff->m_vars) {
        f.push_back(factor(j, factor_type::VAR));
    }
    return f;
}

}


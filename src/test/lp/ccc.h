/*
  Copyright (c) 2017 Microsoft Corporation
  Author: Nikolaj Bjorner, Lev Nachmanson
*/
#pragma once
#include "util/vector.h"
#include "util/rational.h"
namespace lp {
struct ccc {
    vector<vector<std::pair<rational, unsigned>>> m_ineqs;
    vector<std::map<rational, char>> m_intdomains;
    void add_ineq(vector<std::pair<rational, var_index>> & lhs) {
        vector<std::pair<rational, var_index>>  local_lhs;
        for (auto & p : lhs) {
            local_lhs.push_back(std::make_pair(p.first, 0));
            add_var();
        }
        m_ineqs.push_back(local_lhs);
    }
        
    void add_var() {
        std::map<rational, char> v;
        m_intdomains.push_back(v);
    }
};
}

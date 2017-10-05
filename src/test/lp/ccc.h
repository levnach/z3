/*
  Copyright (c) 2017 Microsoft Corporation
  Author: Nikolaj Bjorner, Lev Nachmanson
*/
#pragma once
#include "util/vector.h"
#include "util/rational.h"
namespace lp {
struct ccc {
    vector<vector<rational>> m_ineqs;
    vector<std::map<rational, char>> m_intdomains;
    void add_ineq(vector<rational> & lhs) {
        m_ineqs.push_back(lhs);
    }
};
}

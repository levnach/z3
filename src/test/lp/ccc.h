/*
  Copyright (c) 2017 Microsoft Corporation
  Author: Nikolaj Bjorner, Lev Nachmanson
*/
#pragma once
#include "util/vector.h"
#include "util/rational.h"
namespace lp {
struct ccc {
    vector<vector<rational>> m_vecs;
    vector<std::map<rational, char>> m_maps;
    void add_ineq(vector<rational> & lhs) {
        m_vecs.push_back(lhs);
    }
};
}

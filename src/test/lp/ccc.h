/*
  Copyright (c) 2017 Microsoft Corporation
  Author: Nikolaj Bjorner, Lev Nachmanson
*/
#pragma once
#include "util/vector.h"
#include "util/rational.h"
namespace lp {
struct ccc {
    vector<std::map<rational, char>> m_maps;
};
}

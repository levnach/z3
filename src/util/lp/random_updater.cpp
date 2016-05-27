/*
  Copyright (c) 2013 Microsoft Corporation. All rights reserved.
  Released under Apache 2.0 license as described in the file LICENSE.

  Author: Lev Nachmanson
*/
#include "util/lp/random_updater.h"
namespace lean {
lean::random_updater::random_updater(lar_solver * lar_solver_par) : m_lar_solver(lar_solver_par) {
}

void random_updater::update() {
}

void random_updater::init(unsigned sz, var_index const * vars) {
}
}
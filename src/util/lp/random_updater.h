/*
Copyright (c) 2013 Microsoft Corporation. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.

Author: Lev Nachmanson
*/
#pragma once
#include <set>
#include <vector>
#include "util/lp/lp_settings.h"
// see http://research.microsoft.com/projects/z3/smt07.pdf
// The class searches for a feasible solution with as many different values of variables as it can find
namespace lean {
class lar_solver; // forward definition
class random_updater {
	std::set<mpq> m_initial_set_of_var_values;
	std::set<unsigned> m_candidate_set_for_random_update;
	lar_solver* m_lar_solver;
public:
	random_updater(lar_solver * lar_solver_par);
	void update();
	void init(unsigned sz, var_index const * vars);
	void init_column(unsigned j);
	
};
}

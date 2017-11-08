/*
  Copyright (c) 2017 Microsoft Corporation
  Author: Nikolaj Bjorner, Lev Nachmanson
*/
#include "util/lp/cut_solver_def.h"
namespace lp {
template <typename T>
T cut_solver<T>::m_local_zero = zero_of_type<T>();
template <> int cut_solver<int>::m_local_zero = 0;
template <> mpq cut_solver<mpq>::m_local_zero = zero_of_type<mpq>();
template lbool cut_solver<mpq>::check();
}


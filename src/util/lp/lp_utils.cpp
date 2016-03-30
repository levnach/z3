/*
  Copyright (c) 2013 Microsoft Corporation. All rights reserved.
  Released under Apache 2.0 license as described in the file LICENSE.

  Author: Lev Nachmanson
*/
#include "util/lp/lp_utils.h"
#ifdef lp_for_z3
namespace lp {
double numeric_traits<double>::g_zero = 0.0;
double numeric_traits<double>::g_one = 1.0;
}
#endif

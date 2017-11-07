/*
  Copyright (c) 2017 Microsoft Corporation
  Author: Lev Nachmanson
*/
#include "util/lp/scaler_def.h"
template bool lp::scaler<double, double>::scale();
template bool lp::scaler<lp::mpq, lp::mpq>::scale();

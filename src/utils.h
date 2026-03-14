#ifndef UTILS_H
#define UTILS_H
//#define DEBUG

#ifdef DEBUG
#include <iostream>
#endif

#ifndef USE_FC_LEN_T
#define USE_FC_LEN_T
#endif
#include <Rconfig.h>
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"
#ifndef FCONE
#define FCONE
#endif

#include "Rcpp.h"

#include <vector>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <algorithm>

/* Defining all R-accessible functions. */

extern "C" {

SEXP compute_apl (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP fit_levenberg (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP get_levenberg_start (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

}

/* Other utility functions and values */

const double low_value=std::pow(10.0, -10.0), log_low_value=std::log(low_value);

const double LNtwo=std::log(2), one_million=1000000, LNmillion=std::log(one_million);

#endif

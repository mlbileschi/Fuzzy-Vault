#ifndef __WELCH_BERLEKAMP_H__
#define __WELCH_BERLEKAMP_H__

#include "NTL/mat_ZZ.h"
#include "NTL/ZZ_pX.h"
NTL_CLIENT

static bool LinearSolve (vec_ZZ_p & solution, mat_ZZ_p & A);
bool WelchBerlekamp (ZZ_pX & result, const vec_ZZ_p & X, const vec_ZZ_p & Y,  int t);

#endif //__WELCH_BERLEKAMP_H__


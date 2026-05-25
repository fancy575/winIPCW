/*
 * Native helpers for winIPCW.
 *
 * These routines are intentionally small and package-facing. The statistical
 * API is written in R, while expensive pairwise kernels can be moved here
 * without changing user-facing functions.
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <cmath>

static inline double as_real_at(SEXP x, R_xlen_t i) {
  return REAL(x)[i];
}

static inline int as_int_at(SEXP x, R_xlen_t i) {
  return INTEGER(x)[i];
}

extern "C" SEXP winIPCW_cpp_version() {
  SEXP out = PROTECT(Rf_allocVector(STRSXP, 1));
  SET_STRING_ELT(out, 0, Rf_mkChar("winIPCW native core 0.0.1"));
  UNPROTECT(1);
  return out;
}

extern "C" SEXP winIPCW_pair_ipcw_sums(SEXP y_trt,
                                     SEXP d_trt,
                                     SEXP y_ctl,
                                     SEXP d_ctl,
                                     SEXP g_trt,
                                     SEXP g_ctl,
                                     SEXP tau,
                                     SEXP eps) {
  if (!Rf_isReal(y_trt) || !Rf_isInteger(d_trt) ||
      !Rf_isReal(y_ctl) || !Rf_isInteger(d_ctl) ||
      !Rf_isReal(g_trt) || !Rf_isReal(g_ctl)) {
    Rf_error("Invalid native input types.");
  }

  if (!Rf_isReal(tau) || XLENGTH(tau) != 1L ||
      !Rf_isReal(eps) || XLENGTH(eps) != 1L) {
    Rf_error("`tau` and `eps` must be numeric scalars.");
  }

  const R_xlen_t n1 = XLENGTH(y_trt);
  const R_xlen_t n0 = XLENGTH(y_ctl);
  if (XLENGTH(d_trt) != n1 || XLENGTH(g_trt) != n1 ||
      XLENGTH(d_ctl) != n0 || XLENGTH(g_ctl) != n0) {
    Rf_error("Native input vectors have incompatible lengths.");
  }
  if (n1 == 0 || n0 == 0) {
    Rf_error("Treatment and control groups must both be non-empty.");
  }

  const double lower = REAL(eps)[0];
  if (!std::isfinite(lower) || lower <= 0.0) {
    Rf_error("`eps` must be a positive finite scalar.");
  }

  SEXP out = PROTECT(Rf_allocVector(REALSXP, 2));
  REAL(out)[0] = 0.0;
  REAL(out)[1] = 0.0;

  for (R_xlen_t i = 0; i < n1; ++i) {
    const double yi = as_real_at(y_trt, i);
    for (R_xlen_t j = 0; j < n0; ++j) {
      const double yj = as_real_at(y_ctl, j);
      const double den = fmax(as_real_at(g_trt, i) * as_real_at(g_ctl, j), lower);
      if (yi > yj && as_int_at(d_ctl, j) == 1) {
        REAL(out)[0] += 1.0 / den;
      } else if (yj > yi && as_int_at(d_trt, i) == 1) {
        REAL(out)[1] += 1.0 / den;
      }
    }
  }

  REAL(out)[0] /= static_cast<double>(n1) * static_cast<double>(n0);
  REAL(out)[1] /= static_cast<double>(n1) * static_cast<double>(n0);

  UNPROTECT(1);
  return out;
}

static const R_CallMethodDef call_methods[] = {
  {"winIPCW_cpp_version", reinterpret_cast<DL_FUNC>(&winIPCW_cpp_version), 0},
  {"winIPCW_pair_ipcw_sums", reinterpret_cast<DL_FUNC>(&winIPCW_pair_ipcw_sums), 8},
  {NULL, NULL, 0}
};

extern "C" void R_init_winIPCW(DllInfo *dll) {
  R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}

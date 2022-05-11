#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */

extern SEXP imeiv_plauscontourMCMC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP imeiv_plauscontourMC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP imeiv_plauscontourMC2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP imeiv_sortmat(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"imeiv_plauscontourMCMC", (DL_FUNC) &imeiv_plauscontourMCMC, 11},
    {"imeiv_plauscontourMC", (DL_FUNC) &imeiv_plauscontourMC, 10},
    {"imeiv_plauscontourMC2", (DL_FUNC) &imeiv_plauscontourMC2, 10},
    {"imeiv_sortmat", (DL_FUNC) &imeiv_sortmat, 2},
    {NULL, NULL, 0}
};

void R_init_imeiv(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

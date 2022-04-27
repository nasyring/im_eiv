#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */

extern SEXP imeiv_plauscontourGF(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP imeiv_plauscontourIM(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP imeiv_plauscontourIMmarg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP imeiv_sortmat(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"imeiv_plauscontourGF", (DL_FUNC) &imeiv_plauscontourGF, 9},
    {"imeiv_plauscontourIM", (DL_FUNC) &imeiv_plauscontourIM, 9},
    {"imeiv_plauscontourIMmarg", (DL_FUNC) &imeiv_plauscontourIMmarg, 7},
    {"imeiv_sortmat", (DL_FUNC) &imeiv_sortmat, 2},
    {NULL, NULL, 0}
};

void R_init_imeiv(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

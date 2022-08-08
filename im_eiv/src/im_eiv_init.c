#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */

extern SEXP imeiv_plauscontourMCMC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP imeiv_plauscontourMCMCcond(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP imeiv_plauscontourSIR(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP imeiv_plauscontourMCMC2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP imeiv_plauscontourMC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP imeiv_plauscontourMC2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP imeiv_sortmat(SEXP, SEXP);
extern SEXP imeiv_grow(SEXP, SEXP);
extern SEXP imeiv_loglik(SEXP, SEXP, SEXP, SEXP);
extern SEXP imeiv_maxloglik(SEXP, SEXP, SEXP, SEXP);
extern SEXP imeiv_genIMplaus(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP imeiv_optimrcpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP imeiv_plaus_mc(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"imeiv_plauscontourMCMC", (DL_FUNC) &imeiv_plauscontourMCMC, 11},
    {"imeiv_plauscontourMCMCcond", (DL_FUNC) &imeiv_plauscontourMCMCcond, 15},
    {"imeiv_plauscontourSIR", (DL_FUNC) &imeiv_plauscontourSIR, 8},
    {"imeiv_plauscontourMCMC2", (DL_FUNC) &imeiv_plauscontourMCMC2, 11},
    {"imeiv_plauscontourMC", (DL_FUNC) &imeiv_plauscontourMC, 10},
    {"imeiv_plauscontourMC2", (DL_FUNC) &imeiv_plauscontourMC2, 10},
    {"imeiv_sortmat", (DL_FUNC) &imeiv_sortmat, 2},
    {"imeiv_grow", (DL_FUNC) &imeiv_grow, 1},
    {"imeiv_loglik", (DL_FUNC) &imeiv_loglik, 4},
    {"imeiv_maxloglik", (DL_FUNC) &imeiv_maxloglik, 4},
    {"imeiv_genIMplaus", (DL_FUNC) &imeiv_genIMplaus, 5},
    {"imeiv_optimrcpp", (DL_FUNC) &imeiv_genIMplaus, 4},
    {"imeiv_plaus_mc", (DL_FUNC) &imeiv_plaus_mc, 5},
    {NULL, NULL, 0}
};

void R_init_imeiv(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

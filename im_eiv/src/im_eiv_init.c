#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP impred_randsetsMCMC(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_randsetspred(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_sigmaSolve(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_zeroin(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP impred_root_function(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"impred_randsetsMCMC", (DL_FUNC) &impred_randsetsMCMC, 5},
    {"impred_randsetspred", (DL_FUNC) &impred_randsetspred, 8},
    {"impred_sigmaSolve", (DL_FUNC) &impred_sigmaSolve, 5}, 
    {"impred_zeroin", (DL_FUNC) &impred_zeroin, 8},
    {"impred_root_function", (DL_FUNC) &impred_root_function, 8}, 
    {NULL, NULL, 0}
};

void R_init_impred(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void odin_model_endectocide_initmod_desolve(void *);
extern void odin_model_endectocide_output_dde(void *);
extern void odin_model_endectocide_rhs_dde(void *);
extern void odin_model_endectocide_rhs_desolve(void *);

/* .Call calls */
extern SEXP odin_model_endectocide_contents(SEXP);
extern SEXP odin_model_endectocide_create(SEXP);
extern SEXP odin_model_endectocide_initial_conditions(SEXP, SEXP);
extern SEXP odin_model_endectocide_metadata(SEXP);
extern SEXP odin_model_endectocide_rhs_r(SEXP, SEXP, SEXP);
extern SEXP odin_model_endectocide_set_initial(SEXP, SEXP, SEXP, SEXP);
extern SEXP odin_model_endectocide_set_user(SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"odin_model_endectocide_initmod_desolve", (DL_FUNC) &odin_model_endectocide_initmod_desolve, 1},
    {"odin_model_endectocide_output_dde",      (DL_FUNC) &odin_model_endectocide_output_dde,      1},
    {"odin_model_endectocide_rhs_dde",         (DL_FUNC) &odin_model_endectocide_rhs_dde,         1},
    {"odin_model_endectocide_rhs_desolve",     (DL_FUNC) &odin_model_endectocide_rhs_desolve,     1},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"odin_model_endectocide_contents",           (DL_FUNC) &odin_model_endectocide_contents,           1},
    {"odin_model_endectocide_create",             (DL_FUNC) &odin_model_endectocide_create,             1},
    {"odin_model_endectocide_initial_conditions", (DL_FUNC) &odin_model_endectocide_initial_conditions, 2},
    {"odin_model_endectocide_metadata",           (DL_FUNC) &odin_model_endectocide_metadata,           1},
    {"odin_model_endectocide_rhs_r",              (DL_FUNC) &odin_model_endectocide_rhs_r,              3},
    {"odin_model_endectocide_set_initial",        (DL_FUNC) &odin_model_endectocide_set_initial,        4},
    {"odin_model_endectocide_set_user",           (DL_FUNC) &odin_model_endectocide_set_user,           2},
    {NULL, NULL, 0}
};

void R_init_ivRmectin(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

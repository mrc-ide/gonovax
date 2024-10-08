#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void model_initmod_desolve(void *);
extern void model_output_dde(void *);
extern void model_rhs_dde(void *);
extern void model_rhs_desolve(void *);
extern void model_trial_initmod_desolve(void *);
extern void model_trial_output_dde(void *);
extern void model_trial_rhs_dde(void *);
extern void model_trial_rhs_desolve(void *);
extern void model_trial_stochastic_rhs_dde(void *);
extern void model_withouthistory_initmod_desolve(void *);
extern void model_withouthistory_output_dde(void *);
extern void model_withouthistory_rhs_dde(void *);
extern void model_withouthistory_rhs_desolve(void *);
extern void model_withoutPN_initmod_desolve(void *);
extern void model_withoutPN_output_dde(void *);
extern void model_withoutPN_rhs_dde(void *);
extern void model_withoutPN_rhs_desolve(void *);

/* .Call calls */
extern SEXP model_contents(void *);
extern SEXP model_create(void *);
extern SEXP model_initial_conditions(void *, void *);
extern SEXP model_metadata(void *);
extern SEXP model_rhs_r(void *, void *, void *);
extern SEXP model_set_initial(void *, void *, void *, void *);
extern SEXP model_set_user(void *, void *);
extern SEXP model_trial_contents(void *);
extern SEXP model_trial_create(void *);
extern SEXP model_trial_initial_conditions(void *, void *);
extern SEXP model_trial_metadata(void *);
extern SEXP model_trial_rhs_r(void *, void *, void *);
extern SEXP model_trial_set_initial(void *, void *, void *, void *);
extern SEXP model_trial_set_user(void *, void *);
extern SEXP model_trial_stochastic_contents(void *);
extern SEXP model_trial_stochastic_create(void *);
extern SEXP model_trial_stochastic_initial_conditions(void *, void *);
extern SEXP model_trial_stochastic_metadata(void *);
extern SEXP model_trial_stochastic_rhs_r(void *, void *, void *);
extern SEXP model_trial_stochastic_set_initial(void *, void *, void *, void *);
extern SEXP model_trial_stochastic_set_user(void *, void *);
extern SEXP model_withouthistory_contents(void *);
extern SEXP model_withouthistory_create(void *);
extern SEXP model_withouthistory_initial_conditions(void *, void *);
extern SEXP model_withouthistory_metadata(void *);
extern SEXP model_withouthistory_rhs_r(void *, void *, void *);
extern SEXP model_withouthistory_set_initial(void *, void *, void *, void *);
extern SEXP model_withouthistory_set_user(void *, void *);
extern SEXP model_withoutPN_contents(void *);
extern SEXP model_withoutPN_create(void *);
extern SEXP model_withoutPN_initial_conditions(void *, void *);
extern SEXP model_withoutPN_metadata(void *);
extern SEXP model_withoutPN_rhs_r(void *, void *, void *);
extern SEXP model_withoutPN_set_initial(void *, void *, void *, void *);
extern SEXP model_withoutPN_set_user(void *, void *);

static const R_CMethodDef CEntries[] = {
    {"model_initmod_desolve",                (DL_FUNC) &model_initmod_desolve,                1},
    {"model_output_dde",                     (DL_FUNC) &model_output_dde,                     1},
    {"model_rhs_dde",                        (DL_FUNC) &model_rhs_dde,                        1},
    {"model_rhs_desolve",                    (DL_FUNC) &model_rhs_desolve,                    1},
    {"model_trial_initmod_desolve",          (DL_FUNC) &model_trial_initmod_desolve,          1},
    {"model_trial_output_dde",               (DL_FUNC) &model_trial_output_dde,               1},
    {"model_trial_rhs_dde",                  (DL_FUNC) &model_trial_rhs_dde,                  1},
    {"model_trial_rhs_desolve",              (DL_FUNC) &model_trial_rhs_desolve,              1},
    {"model_trial_stochastic_rhs_dde",       (DL_FUNC) &model_trial_stochastic_rhs_dde,       1},
    {"model_withouthistory_initmod_desolve", (DL_FUNC) &model_withouthistory_initmod_desolve, 1},
    {"model_withouthistory_output_dde",      (DL_FUNC) &model_withouthistory_output_dde,      1},
    {"model_withouthistory_rhs_dde",         (DL_FUNC) &model_withouthistory_rhs_dde,         1},
    {"model_withouthistory_rhs_desolve",     (DL_FUNC) &model_withouthistory_rhs_desolve,     1},
    {"model_withoutPN_initmod_desolve",      (DL_FUNC) &model_withoutPN_initmod_desolve,      1},
    {"model_withoutPN_output_dde",           (DL_FUNC) &model_withoutPN_output_dde,           1},
    {"model_withoutPN_rhs_dde",              (DL_FUNC) &model_withoutPN_rhs_dde,              1},
    {"model_withoutPN_rhs_desolve",          (DL_FUNC) &model_withoutPN_rhs_desolve,          1},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"model_contents",                            (DL_FUNC) &model_contents,                            1},
    {"model_create",                              (DL_FUNC) &model_create,                              1},
    {"model_initial_conditions",                  (DL_FUNC) &model_initial_conditions,                  2},
    {"model_metadata",                            (DL_FUNC) &model_metadata,                            1},
    {"model_rhs_r",                               (DL_FUNC) &model_rhs_r,                               3},
    {"model_set_initial",                         (DL_FUNC) &model_set_initial,                         4},
    {"model_set_user",                            (DL_FUNC) &model_set_user,                            2},
    {"model_trial_contents",                      (DL_FUNC) &model_trial_contents,                      1},
    {"model_trial_create",                        (DL_FUNC) &model_trial_create,                        1},
    {"model_trial_initial_conditions",            (DL_FUNC) &model_trial_initial_conditions,            2},
    {"model_trial_metadata",                      (DL_FUNC) &model_trial_metadata,                      1},
    {"model_trial_rhs_r",                         (DL_FUNC) &model_trial_rhs_r,                         3},
    {"model_trial_set_initial",                   (DL_FUNC) &model_trial_set_initial,                   4},
    {"model_trial_set_user",                      (DL_FUNC) &model_trial_set_user,                      2},
    {"model_trial_stochastic_contents",           (DL_FUNC) &model_trial_stochastic_contents,           1},
    {"model_trial_stochastic_create",             (DL_FUNC) &model_trial_stochastic_create,             1},
    {"model_trial_stochastic_initial_conditions", (DL_FUNC) &model_trial_stochastic_initial_conditions, 2},
    {"model_trial_stochastic_metadata",           (DL_FUNC) &model_trial_stochastic_metadata,           1},
    {"model_trial_stochastic_rhs_r",              (DL_FUNC) &model_trial_stochastic_rhs_r,              3},
    {"model_trial_stochastic_set_initial",        (DL_FUNC) &model_trial_stochastic_set_initial,        4},
    {"model_trial_stochastic_set_user",           (DL_FUNC) &model_trial_stochastic_set_user,           2},
    {"model_withouthistory_contents",             (DL_FUNC) &model_withouthistory_contents,             1},
    {"model_withouthistory_create",               (DL_FUNC) &model_withouthistory_create,               1},
    {"model_withouthistory_initial_conditions",   (DL_FUNC) &model_withouthistory_initial_conditions,   2},
    {"model_withouthistory_metadata",             (DL_FUNC) &model_withouthistory_metadata,             1},
    {"model_withouthistory_rhs_r",                (DL_FUNC) &model_withouthistory_rhs_r,                3},
    {"model_withouthistory_set_initial",          (DL_FUNC) &model_withouthistory_set_initial,          4},
    {"model_withouthistory_set_user",             (DL_FUNC) &model_withouthistory_set_user,             2},
    {"model_withoutPN_contents",                  (DL_FUNC) &model_withoutPN_contents,                  1},
    {"model_withoutPN_create",                    (DL_FUNC) &model_withoutPN_create,                    1},
    {"model_withoutPN_initial_conditions",        (DL_FUNC) &model_withoutPN_initial_conditions,        2},
    {"model_withoutPN_metadata",                  (DL_FUNC) &model_withoutPN_metadata,                  1},
    {"model_withoutPN_rhs_r",                     (DL_FUNC) &model_withoutPN_rhs_r,                     3},
    {"model_withoutPN_set_initial",               (DL_FUNC) &model_withoutPN_set_initial,               4},
    {"model_withoutPN_set_user",                  (DL_FUNC) &model_withoutPN_set_user,                  2},
    {NULL, NULL, 0}
};

void R_init_gonovax(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

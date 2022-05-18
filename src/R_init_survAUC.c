//
// 2022-05-18. Author F. Bertrand <frederic.bertrand@utt.fr>
//

#include <stdlib.h> // for NULL
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "surv_measures.h" // for NULL
#include "survAUC_UNO.h" // for NULL
#include "survAUC_Cham_Diao.h" // for NULL
#include "survAUC_Hung_Chiang.h" // for NULL
#include "survAUC_HZ.h" // for NULL
#include "survAUC_SongZhou.h" // for NULL
#include "utils.h" // for NULL

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_NativePrimitiveArgType C_partLCox_t[] = {
  REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, REALSXP 
};

static R_NativePrimitiveArgType C_partLCoxIndiv_t[] = {
  REALSXP, REALSXP, INTSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType C_begg_t[] = {
  REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP
};

static R_NativePrimitiveArgType C_GHCI_t[] = {
  REALSXP, INTSXP, REALSXP
};

static R_NativePrimitiveArgType C_XO_t[] = {
  REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType C_UnoC_t[] = {
  REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP
};

static R_NativePrimitiveArgType C_sens_uno_t[] = {
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

static R_NativePrimitiveArgType C_spec_uno_t[] = {
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP
};

static R_NativePrimitiveArgType C_auc_uno_t[] = {
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP
};

static R_NativePrimitiveArgType C_int_auc_t[] = {
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP
};

static R_NativePrimitiveArgType C_Hung_Chiang_t[] = {
  REALSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType C_km_weight_t[] = {
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP
};

static R_NativePrimitiveArgType C_km_Daim_t[] = {
  REALSXP, REALSXP, REALSXP, INTSXP
};

static R_NativePrimitiveArgType C_cox_weights_t[] = {
  REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, INTSXP
};

static R_NativePrimitiveArgType C_step_eval_R_t[] = {
  REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP 
};

static const R_CMethodDef cMethods[] = {
  {"C_partLCox", (DL_FUNC) &C_partLCox, 6, C_partLCox_t},
  {"C_partLCoxIndiv", (DL_FUNC) &C_partLCoxIndiv, 5, C_partLCoxIndiv_t},
  {"C_begg", (DL_FUNC) &C_begg, 11, C_begg_t},
  {"C_GHCI", (DL_FUNC) &C_GHCI, 3, C_GHCI_t},
  {"C_XO", (DL_FUNC) &C_XO, 6, C_XO_t},
  {"C_UnoC", (DL_FUNC) &C_UnoC, 10, C_UnoC_t},
  {"C_sens_uno", (DL_FUNC) &C_sens_uno, 12, C_sens_uno_t},
  {"C_spec_uno", (DL_FUNC) &C_spec_uno, 8, C_spec_uno_t},
  {"C_auc_uno", (DL_FUNC) &C_auc_uno, 15, C_auc_uno_t},
  {"C_int_auc", (DL_FUNC) &C_int_auc, 9, C_int_auc_t},
  {"C_Hung_Chiang", (DL_FUNC) &C_Hung_Chiang, 12, C_Hung_Chiang_t},
  {"C_km_weight", (DL_FUNC) &C_km_weight, 6, C_km_weight_t},
  {"C_km_Daim", (DL_FUNC) &C_km_Daim, 4, C_km_Daim_t},
  {"C_cox_weights", (DL_FUNC) &C_cox_weights, 6, C_cox_weights_t},
  {"C_step_eval_R", (DL_FUNC) &C_step_eval_R, 6, C_step_eval_R_t},
  {NULL, NULL, 0, NULL}
};

SEXP C_sens_SZ(SEXP THRESH, SEXP T, SEXP STIME, SEXP EVENT, SEXP N_TIME, SEXP LP, 
               SEXP N_LP, SEXP LPNEW, SEXP N_LPNEW, SEXP TYPE_SENS);
SEXP C_spec_SZ(SEXP THRESH, SEXP T, SEXP STIME, SEXP EVENT, SEXP N_TIME, SEXP LP, 
               SEXP N_LP, SEXP LPNEW, SEXP N_LPNEW);
SEXP C_auc_SZ(SEXP THRESH, SEXP T, SEXP STIME, SEXP EVENT, SEXP N_TIME, 
              SEXP STIME_NEW, SEXP EVENT_NEW, SEXP N_TIME_NEW, 
              SEXP LP, SEXP N_LP, SEXP LPNEW, SEXP N_LPNEW, SEXP TYPE_SENS);

static const R_CallMethodDef callMethods[]  = {
  {"C_predError", (DL_FUNC) &C_predError, 14},
  {"C_Cham_Diao", (DL_FUNC) &C_Cham_Diao, 11},
  {"C_sens_SZ", (DL_FUNC) &C_sens_SZ, 10},
  {"C_spec_SZ", (DL_FUNC) &C_spec_SZ, 9},
  {"C_auc_SZ", (DL_FUNC) &C_auc_SZ, 13},
  {"C_survfit_cox", (DL_FUNC) &C_survfit_cox, 7},
  {NULL, NULL, 0}
};

void R_init_survAUC(DllInfo * info)
{
  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}


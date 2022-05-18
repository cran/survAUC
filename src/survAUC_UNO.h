/*
 *  survAUC_UNO.h
 *  Daim
 *
 *  Created by Sergej Potapov on 01.06.10.
 *  Copyright 2010 __IMBE__. All rights reserved.
 *
 */


void C_spec_uno( double *spez, double *thres, double *t, double *marker, double *new_data, 
			  int *n_th, int *n_t, int *n_new_data );
void C_sens_uno( double *sens, double *surv_time, double *status, double *thres, double *t, 
			  double *marker, double *new_surv, double *new_event, int *n_th, int *n_t, 
			  int *n_new_data, int *n_surv );
void C_auc_uno( double *auc, double *i_auc, double *sens, double *spec, double *surv_time,
			 double *status, double *thres, double *t, double *marker, double *new_surv_t,
			 double *new_event, int *n_th, int *n_t, int *n_new_data, int *n_surv);
void C_int_auc( double *i_auc, double *auc, double *time, double *S, 
			 double *tmax, int *n_S, int *wChoice, int *maxI, double *maxS, int *Con_Inc);


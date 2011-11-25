/*
 *  surv_measures.h
 *  Daim
 *
 *  Created by Sergej Potapov on 11.10.10.
 *  Copyright 2010 __IMBE__. All rights reserved.
 *
 */


void GHCI(double *lp, int *n_lp, double *ans);
void cens_weights(double *times, int *n_times, double *stime, double *event, int *n_stime,
				  double *stime_new, double *event_new, int *n_stime_new, double *weights);
SEXP predError(SEXP TIME, SEXP EVENT, SEXP N_TIME,
			   SEXP TIME_NEW, SEXP EVENT_NEW, SEXP N_TIME_NEW, 
			   SEXP TH_TIME, SEXP N_TH_TIME, SEXP LP, SEXP N_LP,
			   SEXP LPNEW, SEXP N_LPNEW, SEXP TYPE, SEXP INTEG);
void partLCoxIndiv(double *stime, double *time, int *n_stime, double *lp, double *LL);
void XO(double *stime, double *event, int *n_stime, double *lp, double *lp0, double *XO);
void UnoC(double *stime, double *event, int *n_stime, double *new_stime,
		  double *new_event, int *new_n_stime, double *lp, double *tau,
		  int *n_tau, double *CStat);
void c_begg(double *new_stime, double *new_event, int *new_n_stime,
			  double *times, int *n_times, double *lp, double *lpnew, 
			  double *surv_prob, double *surv_times, int *n_surv_times, double *CStat)

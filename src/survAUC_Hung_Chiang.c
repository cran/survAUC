/*
 *  survAUC_Hung_Chiang.c
 *  Daim
 *
 *  Created by Sergej Potapov on 10.10.10.
 *  Copyright 2010 __IMBE__. All rights reserved.
 *
 */



#include "utils.h"




void C_Hung_Chiang(double *time, int *n_time, double *stime, double *event, int *n_stime,
				 double *stime_new, double *event_new, int *n_stime_new,
				 double *lpnew, int *n_lpnew, double *ans, double *i_auc)
{
	int i, j, k;

	double *SX, *SXa, *St, *Sta, *Sc, *Sca; 
	double *tmp_status, *Sc_stime, *Sc_event;
	double tmp_nen = 0.0;
	Sc_stime = Calloc(*n_stime, double);
	Sc_event = Calloc(*n_stime, double);
	SX = Calloc(*n_stime, double);
	St = Calloc(*n_stime, double);
	Sc = Calloc(*n_stime, double);
	tmp_status = Calloc(*n_stime, double);
	SXa = Calloc(*n_time, double);
	Sta = Calloc(*n_time, double);
	Sca = Calloc(*n_stime_new, double);
	
	for(i=0; i<*n_stime; i++){
		tmp_status[i] = 1.0;
		Sc_stime[i] = stime[i];
		Sc_event[i] = 1.0 - event[i];
	}
	
	C_km_Daim(St, stime, event, n_stime);
	step_eval2(Sta, time, St, stime, *n_time, *n_stime);
	
	C_km_Daim(SX, stime, tmp_status, n_stime);
	step_eval2(SXa, time, SX, stime, *n_time, *n_stime);
	
	C_km_Daim(Sc, Sc_stime, Sc_event, n_stime);
	step_eval2(Sca, stime_new, Sc, Sc_stime, *n_stime_new, *n_stime);

	/* Calculation of AUC */
	for(k=0; k<*n_time; k++){
		for(i=0; i<*n_lpnew; i++){
			for(j=0; j<*n_lpnew; j++){
				if(i != j && ((event_new[i] && (lpnew[i] > lpnew[j]) && (stime_new[i] <= time[k] && stime_new[j] > time[k])) && (Sca[i] > FLT_EPSILON))){
					ans[k] += 1.0 / (Sca[i]);
				}
			}
		}
		tmp_nen = SXa[k]*(1.0-Sta[k])*(*n_lpnew)*(*n_lpnew-1);
		if(tmp_nen > FLT_EPSILON)
			ans[k] /= tmp_nen;
		else
			ans[k] = 0.0;
	}
	Free(SX);Free(SXa);Free(Sc);Free(Sca);Free(St);
	Free(Sta);Free(Sc_event);Free(Sc_stime);Free(tmp_status);
	/* Calculation of iAUC */
	
	double *f, *S, *S_new;
	f = Calloc(*n_time, double);
	S_new = Calloc(*n_stime_new, double);
	S = Calloc(*n_time, double);
	
	C_km_Daim(S_new, stime_new, event_new, n_stime_new);
	step_eval2(S, time, S_new, stime_new, *n_time,  *n_stime_new);
	
	f[0] = 1.0 - S[0];
	for(i=1; i<*n_time; i++){
		f[i] = S[i-1] - S[i];
	}
	double wT = 0.0;
	for(i=0; i < *n_time; i++){
		if(f[i] > FLT_EPSILON){
			wT += f[i];
		}
	}
	for(i=0; i < *n_time; i++){
		if(wT != 0.0){
			/* cumulative case*/
			if(f[i] > FLT_EPSILON && R_finite(ans[i]) ){
				*i_auc += ans[i] * (f[i]) / wT;
			}
		}
	}
	Free(f);Free(S);Free(S_new);
}





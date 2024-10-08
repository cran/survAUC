/*
 *  surv_measures.c
 *  Daim
 *
 *  Created by Sergej Potapov on 11.10.10.
 *  Copyright 2010 __IMBE__. All rights reserved.
 *  2022-05-18. Updated by F. Bertrand <frederic.bertrand@utt.fr>
 *
 */

#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>


#include "utils.h"

/*  Calculation of Gonen and Hellers Concordance Index for Cox models */
void C_GHCI(double *lp, int *n_lp, double *ans){
	int i,j;
	double temp_m = 0.;
	double K = 0.;
	
	for(i = 0; i < *n_lp-1; i++){
		for(j = i+1; j < *n_lp; j++){
			temp_m = lp[i] -lp[j];
			if(temp_m > 0.0)
				K += 1. / (1. + exp(-1*temp_m));
			if(temp_m < 0.0)
				K += 1. / (1. + exp(temp_m));
		}
	}
	*ans = 2.0*K / (double) *n_lp / (double) (*n_lp - 1); 
}





/* Calculation of censoring weights Cox-Model*/
void C_cens_weights(double *times, int *n_times, double *stime, double *event, int *n_stime,
				  double *stime_new, double *event_new, int *n_stime_new, double *weights)
{
	int i, j;
	
	/* Calculation of Kaplan-Maier */
	double *surv;
	surv = R_Calloc(*n_stime, double);
	C_km_Daim(surv, stime, event, n_stime);

	/* Calculation of survival for test data time points */
	double *W_1_2;
	W_1_2 = R_Calloc(*n_stime_new, double);
	step_eval2(W_1_2, stime_new, surv, stime, *n_stime_new, *n_stime);
	
	/* Calculation of survival for new time points */
	double *W_2_2;
	W_2_2 = R_Calloc(*n_times, double);
	step_eval2(W_2_2, times, surv, stime, *n_times, *n_stime);

	
	for(i=0; i < *n_times; i++){
		for(j=0; j < *n_stime_new; j++){
			if(stime_new[j] <= times[i]){
				weights[i+*n_times*j] = event_new[j] / W_1_2[j];
			}else{
				weights[i+*n_times*j] = 1. / W_2_2[i];
			}
		}
	}
	R_Free(surv);R_Free(W_1_2);R_Free(W_2_2);
}




/* Calculation of 'Prediction Error Curve': brier & robust */
SEXP C_predError(SEXP TIME, SEXP EVENT, SEXP N_TIME, 
			   SEXP TIME_NEW, SEXP EVENT_NEW, SEXP N_TIME_NEW, 
			   SEXP TH_TIME, SEXP N_TH_TIME, SEXP LP, SEXP N_LP, 
			   SEXP LPNEW, SEXP N_LPNEW, SEXP TYPE, SEXP INTEG)
{
	SEXP xdims, S1a, ERR;
	int i, j, nrx, ncx;
	
	PROTECT(S1a = C_survfit_cox(LP, TIME, EVENT, N_TIME, N_LP, LPNEW, N_LPNEW));
	xdims = getAttrib(VECTOR_ELT(S1a,0), R_DimSymbol);
	nrx = INTEGER(xdims)[0];
	ncx = INTEGER(xdims)[1];
	
	int n_time_new = INTEGER(N_TIME_NEW)[0];
	int N_th_times = LENGTH(TH_TIME);
	double *surv_new;
	surv_new = R_Calloc(N_th_times*ncx,double);
	
	double *Y;
	Y = R_Calloc(N_th_times*ncx,double);	
	step_eval3(surv_new, REAL(TH_TIME), REAL(VECTOR_ELT(S1a,0)), REAL(VECTOR_ELT(S1a,1)), N_th_times, ncx, nrx);
	
	double *weights;
	weights = R_Calloc(N_th_times*ncx,double);
	
	double *event;
	event = R_Calloc(INTEGER(N_TIME)[0],double);
	for(i=0; i < INTEGER(N_TIME)[0]; i++){
		event[i] = 1. - REAL(EVENT)[i];
	}
	
	C_cens_weights(REAL(TH_TIME), INTEGER(N_TH_TIME), REAL(TIME), event, INTEGER(N_TIME),
				 REAL(TIME_NEW), REAL(EVENT_NEW), INTEGER(N_TIME_NEW), weights);
	R_Free(event);
	
	/* INTEGER(TYPE)[0] = 1 - robust predicion Error */
	/* INTEGER(TYPE)[0] = 0 - brier predicion Error */

	if(INTEGER(TYPE)[0]){
		for(i=0; i < n_time_new; i++){
			for(j=0; j < N_th_times; j++){
				Y[j+i*N_th_times] = fabs((REAL(TIME_NEW)[i] > REAL(TH_TIME)[j]) - surv_new[j+i*N_th_times]) * weights[j+i*N_th_times];
			}
		}
	}
	else{
		for(i=0; i < n_time_new; i++){
			for(j=0; j < N_th_times; j++){
				Y[j+i*N_th_times] = pow(fabs((REAL(TIME_NEW)[i] > REAL(TH_TIME)[j]) - surv_new[j+i*N_th_times]),2) * weights[j+i*N_th_times];
			}
		}
	}
	double temp_err;
	PROTECT(ERR = allocVector(REALSXP,N_th_times));
	for (i = 0; i < N_th_times; i++){
		temp_err = 0.;
		for(j=0; j < n_time_new; j++){
			temp_err += Y[j*N_th_times + i];
		}
		REAL(ERR)[i] = temp_err / (double) n_time_new;
	}
	R_Free(weights);R_Free(Y);R_Free(surv_new);
	
	SEXP IERR, result, names_result;
	PROTECT(IERR = allocVector(REALSXP,1));
	/*INTEGER(INTEG)[0] = 1  - integration without weighting*/
	/*INTEGER(INTEG)[0] = 0  - integration with weighting*/
	REAL(IERR)[0] = 0.;
	if(INTEGER(INTEG)[0]){
		for (i = 0; i < N_th_times-1; i++){
			REAL(IERR)[0] += ((REAL(TH_TIME)[i+1] + REAL(TH_TIME)[i])/2.0) * fabs(REAL(ERR)[i+1] - REAL(ERR)[i]);
		}
		REAL(IERR)[0] /= dmax(REAL(TH_TIME), N_th_times) - dmin(REAL(TH_TIME), N_th_times);
	}else{
		/* Calculation of weighted iERR */
		int n_new_data = INTEGER(N_TIME_NEW)[0];
		double *f, *S, *S_new;
		f = R_Calloc(N_th_times, double);
		S_new = R_Calloc(n_new_data, double);
		S = R_Calloc(N_th_times, double);
		C_km_Daim(S_new, REAL(TIME_NEW), REAL(EVENT_NEW), INTEGER(N_TIME_NEW));
		step_eval2(S, REAL(TH_TIME), S_new, REAL(TIME_NEW), N_th_times, n_new_data);
		
		f[0] = 1.0 - S[0];
		for(i=1; i < N_th_times; i++){
			f[i] = S[i-1] - S[i];
		}
		double wT = 0.0;
		for(i=0; i < N_th_times; i++){
			wT += f[i];
		}
		double i_err = 0.;
		for(i=0; i < N_th_times; i++){
			if(wT != 0.0){
				if(f[i] > FLT_EPSILON)
					i_err += REAL(ERR)[i] * f[i] / wT;
			}
		}
		R_Free(f);R_Free(S);R_Free(S_new);
		REAL(IERR)[0] = i_err;
	}
	
	PROTECT(result = allocVector(VECSXP,3));
	PROTECT(names_result = allocVector(STRSXP, 3));
	SET_STRING_ELT(names_result, 0, mkChar("error"));
	SET_STRING_ELT(names_result, 1, mkChar("times"));
	SET_STRING_ELT(names_result, 2, mkChar("ierror"));
	setAttrib(result, R_NamesSymbol, names_result);
	UNPROTECT(1);
	SET_VECTOR_ELT(result, 0, ERR);
	SET_VECTOR_ELT(result, 1, TH_TIME);
	SET_VECTOR_ELT(result, 2, IERR);
	UNPROTECT(4);
	return(result);
}	
	



/* Calculation of 'Partial Likelihood Cox-Model */
void C_partLCox(double *time, double *event, int *n_time, double *lp, int *n_lp, double *LL)
{
	int i,j;
	double tmp_risk;
	
	double *risk;
	risk = R_Calloc(*n_time, double);
	for(i=0; i < *n_time; i++){
		tmp_risk=0.;
		for(j=0; j < *n_time; j++){
			if(time[j] >= time[i])
				tmp_risk += exp(lp[j]);
		}
		risk[i] = tmp_risk;
	}
	double *f;
	f = R_Calloc(*n_time, double);
	for(i=0; i < *n_time; i++){
		f[i] = lp[i]*event[i];
	}
	rsort_xyzv(time, event, risk, f, *n_time);
	
	double time0 = time[0];
	for (i = 1, j = 0; i < *n_time; i++){
		if(fabs(time0 - time[i]) > DBL_EPSILON){
			time0 = time[i];
			j++;
			event[j] = event[i];
			f[j] = f[i];
		}else{
			time[j] = time[i];
			event[j] += event[i];
			f[j] += f[i];
			risk[j] = risk[i];
		}
	}
	for(i=0; i < j+1; i++){
		*LL += f[i] - event[i]*log(risk[i]);
	}
	R_Free(risk);R_Free(f);
}




/* Calculation of 'Partial Likelihood Indiv */
void C_partLCoxIndiv(double *stime, double *time, int *n_stime, double *lp, double *LL){
	int i;
	double nen = 0.0;
	for(i=0; i < *n_stime; i++){
		if(stime[i] >= *time){
			nen += exp(lp[i]);
		}
	}
	for(i=0; i < *n_stime; i++){
		if(stime[i] >= *time){
			LL[i] = exp(lp[i]) / nen;
		}else{
			LL[i] = 0.;
		}
	}
}



/* Calculation of 'XO' measure */
void C_XO(double *stime, double *event, int *n_stime, double *lp, double *lp0, double *XO)
{
	int i,j;
	double *sF;
	sF = R_Calloc(*n_stime, double);
	
	double *pijbeta;
	pijbeta = R_Calloc(*n_stime, double);
	
	double *pij0;
	pij0 = R_Calloc(*n_stime, double);
	
	double tmp_sum;
	for(i=0; i<*n_stime; i++){
		C_partLCoxIndiv(stime, &stime[i], n_stime, lp, pijbeta);
		C_partLCoxIndiv(stime, &stime[i], n_stime, lp0, pij0);
		tmp_sum=0.;
		for(j=0; j<*n_stime; j++){
			if(pij0[j] > 0.){
				tmp_sum += pijbeta[j] * log(pijbeta[j] / pij0[j]);
			}
		}
		sF[i] = tmp_sum;
	}
	R_Free(pijbeta);R_Free(pij0);
	double *surv;
	surv = R_Calloc(*n_stime, double);
	C_km_Daim(surv, stime, event, n_stime);
	
	for(i=*n_stime-1; i>0; i--){
		surv[i] = surv[i-1] - surv[i];
	}
	surv[0] = 1. - surv[0];
	
	double GammaHat=0.0;
	for(i=0; i<*n_stime; i++){
		GammaHat += surv[i] * sF[i];
	}
	*XO = 1. - exp(-2*GammaHat);
	R_Free(sF);R_Free(surv);
}






/* Calculation of 'C-Statistic' suggest by Uno */

void C_UnoC(double *stime, double *event, int *n_stime, double *new_stime, 
            double *new_event, int *new_n_stime, double *lp, double *tau, 
            int *n_tau, double *CStat)
{
	double *surv;
	surv = R_Calloc(*n_stime, double);
	C_km_Daim(surv, stime, event, n_stime);
	
	double *G;
	G = R_Calloc(*new_n_stime, double);
	step_eval2(G, new_stime, surv, stime, *new_n_stime, *n_stime);
	
	int i,j;
	if(*n_tau < 2){
		double denom=0.0, num=0.0;
		for(i=0; i<*new_n_stime; i++){
			for(j=0; j<*new_n_stime; j++){
				if(new_stime[i] < new_stime[j] && G[i] > 0.0){
					denom += (1.0/G[i]/G[i])*(new_event[i])*(new_stime[i] < *tau);
					num += (1.0/G[i]/G[i])*(new_event[i])*(new_stime[i] < *tau)*(lp[i] > lp[j]);// + 0.5*(lp[i] == lp[j]) - 0.5*(i == j));
				}
			}
		}
		*CStat = num / denom;
	}
	else{
		double *denom;
		denom = R_Calloc(*n_tau, double);
		double *num;
		num = R_Calloc(*n_tau, double);
		int k=0;
		for(k=0; k<*n_tau; k++){
			denom[k]=0.0, num[k]=0.0;
			for(i=0; i<*new_n_stime; i++){
				for(j=0; j<*new_n_stime; j++){
					if(new_stime[i] < new_stime[j] && G[i] > 0.0){
						denom[k] += (1.0/G[i]/G[i])*(new_event[i])*(new_stime[i] < *tau);
						num[k] +=  (1.0/G[i]/G[i])*(new_event[i])*(new_stime[i] < *tau)*(lp[i] > lp[j]);
					}
				}
			}
			if(denom[k] <= 0){
				CStat[k]=0.0;
			}else{
				CStat[k] = num[k]/denom[k];
			}
		}
		R_Free(denom);R_Free(num);
	}
	R_Free(surv);R_Free(G);
}




/* Calculation of 'C-Statistic' suggest by Begg */

void C_begg(double *new_stime, double *new_event, int *new_n_stime,
			  double *times, int *n_times, double *lp, double *lpnew, 
			  double *surv_prob, double *surv_times, int *n_surv_times, double *CStat)
{
	int i, j;
	double *surv_new;
	surv_new = R_Calloc((*n_times)*(*new_n_stime),double);
	step_eval3(surv_new, times, surv_prob, surv_times, *n_times, *new_n_stime, *n_surv_times);
	
	double tempC=0.0, Cindex=0.0;
	for(i=0; i < *new_n_stime; i++){
		for(j=0; j < *new_n_stime; j++){
			tempC=0.0;
			if(fabs(lpnew[i] - lpnew[j]) <= FLT_EPSILON){
				tempC = 0.5;
			}else{
				if(lpnew[i] > lpnew[j]){
					if((new_event[i] == 1.0) && (new_event[j] == 1.0) && (new_stime[i] < new_stime[j])){
						tempC = 1.0;
					}
					if((new_event[i] == 0.0) && (new_event[j] == 1.0) && (new_stime[i] < new_stime[j])){
						if(surv_new[i+*new_n_stime*i] > FLT_EPSILON){
							tempC = (surv_new[i+*new_n_stime*i] - surv_new[j+*new_n_stime*i])/surv_new[i+*new_n_stime*i];
						}else{
							tempC = 0.0;
						}
					}
					if((new_event[i] == 1.0) && (new_event[j] == 0.0) && (new_stime[i] < new_stime[j])){
						tempC = 1.0;
					}
					if((new_event[i] == 1.0) && (new_event[j] == 0.0) && (new_stime[i] > new_stime[j])){
						if(surv_new[j+*new_n_stime*j] > FLT_EPSILON){
							tempC = surv_new[i+*new_n_stime*j]/surv_new[j+*new_n_stime*j];
						}else{
							tempC = 0.0;
						}
					}
					if((new_event[i] == 0.0) && (new_event[j] == 0.0) && (new_stime[i] < new_stime[j])){
						if(surv_new[i+*new_n_stime*i] > FLT_EPSILON){
							tempC = (surv_new[i+*new_n_stime*i] - surv_new[j+*new_n_stime*i]/2.0)/surv_new[i+*new_n_stime*i];
						}else{
							tempC = 0.0;
						}
					}
					if((new_event[i] == 0.0) && (new_event[j] == 0.0) && (new_stime[i] > new_stime[j])){
						if(surv_new[j+*new_n_stime*j] > FLT_EPSILON){
							tempC = surv_new[i+*new_n_stime*j]/surv_new[j+*new_n_stime*j]/2.0;
						}else{
							tempC = 0.0;
						}
					}
				}else{
					if(lpnew[i] < lpnew[j]){
						tempC = 0.0;
					}
				}
			}
			if((fabs(lpnew[i] - lpnew[j]) <= FLT_EPSILON) && (i == j)){
				tempC = tempC/2.0;
			}
			if(i == j){
				tempC = 0.0;
			}
			Cindex = Cindex + tempC;
		}
	}
	*CStat = Cindex / (*new_n_stime * (*new_n_stime - 1.0) / 2.0);
	R_Free(surv_new);
}














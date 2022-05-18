/*
 *  survAUC_HZ.c
 *  Daim
 *
 *  Created by Sergej Potapov on 22.06.10.
 *  Copyright 2010 __IMBE__. All rights reserved.
 *
 */

#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>


#include "utils.h"

void C_cox_weights (double *eta, double *time, int *status, double *target, double *AUC, int *N_st)
{
	int i, j=0, n=0, n_dead=0;
	LDOUBLE sum_e_eta=0.0;
	
	int *ID_dead;
	int *ID_at_risk;
	ID_dead = Calloc(*N_st, int);
	ID_at_risk = Calloc(*N_st, int);
	
	
	for (i=0, j=0; i<*N_st ; i++) {
		ID_at_risk[i] = time[i] >= *target;
		ID_dead[i] = time[i] == *target && status[i] == 1;
				
		if(ID_at_risk[i]){
			sum_e_eta += exp(eta[i]);
			eta[j] = eta[i];
			ID_dead[j] = ID_dead[i];
			j++;
		}
		n += ID_at_risk[i];
		n_dead += ID_dead[i];
	}
	
	LDOUBLE *P0;
	LDOUBLE *P1;
//	double *P0;
//	double *P1;
	int *indx;
	P0 = Calloc(n, LDOUBLE);
	P1 = Calloc(n, LDOUBLE);
//	P0 = Calloc(n, double);
//	P1 = Calloc(n, double);
	indx = Calloc(n, int);
	
	for(i=0; i < n; i++){
		indx[i] = i;
		if(ID_dead[i]){
			P0[i] = 0.;
		}else{
			P0[i] = (LDOUBLE) (1.0 / (double) (n - n_dead));
		}
		P1[i] = (LDOUBLE) exp(eta[i]) / sum_e_eta;
	}
	Free(ID_dead); Free(ID_at_risk);

	double *FP;
	double *TP;

	FP = Calloc(n+2, double);
	TP = Calloc(n+2, double);
	
	FP[0] = 1.0; TP[0] = 1.0;
	rsort_index(eta, indx, n);
	LDOUBLE csum=0., csum1=0.;
//	double csum=0.0, csum1=0.0;
	
	for(i=1; i < n; i++){
		csum += P1[indx[i]];
		TP[i] = 1.0 - csum;
		csum1 += P0[indx[i]];
		FP[i] = 1.0 - csum1;
	}
	Free(indx);Free(P0);Free(P1);
	
	for(i=0; i < (n+1); i++){
		*AUC += 0.5*(TP[i+1] + TP[i]) * fabs(FP[i+1] - FP[i]);
	}
}

#' @title AUC estimator proposed by Uno et al.
#' 
#' @rdname survAUC_Uno
#' @name survAUC_Uno
#' 
#' @description Uno's estimator of cumulative/dynamic AUC for right-censored time-to-event
#' data
#' 
#' @details The \code{sens.uno} and \code{spec.uno} functions implement the estimators
#' of time-dependent true and false positive rates proposed in Section 5.1 of
#' Uno et al. (2007).
#' 
#' The \code{AUC.uno} function implements the estimator of cumulative/dynamic
#' AUC that is based on the TPR and FPR estimators proposed by Uno et al.
#' (2007). It is given by the area(s) under the time-dependent ROC curve(s)
#' estimated by \code{sens.uno} and \code{spec.uno}. The \code{iauc} summary
#' measure is given by the integral of AUC on [0, max(\code{times})] (weighted
#' by the estimated probability density of the time-to-event outcome).
#' 
#' Uno's estimators are based on inverse-probability-of-censoring weights and
#' do not assume a specific working model for deriving the predictor
#' \code{lpnew}. It is assumed, however, that there is a one-to-one
#' relationship between the predictor and the expected survival times
#' conditional on the predictor. Note that the estimators implemented in
#' \code{sens.uno}, \code{spec.uno} and \code{AUC.uno} are restricted to
#' situations where the random censoring assumption holds.
#' 
#' @aliases AUC.uno spec.uno sens.uno
#' @param Surv.rsp A \code{Surv(.,.)} object containing to the outcome of the
#' training data.
#' @param Surv.rsp.new A \code{Surv(.,.)} object containing the outcome of the
#' test data.
#' @param lpnew The vector of predictors obtained from the test data.
#' @param times A vector of time points at which to evaluate AUC.
#' @param savesensspec A logical specifying whether sensitivities and
#' specificities should be saved.
#' @return \code{AUC.uno} returns an object of class \code{survAUC}.
#' Specifically, \code{AUC.uno} returns a list with the following components:
#' \item{auc}{The cumulative/dynamic AUC estimates (evaluated at
#' \code{times}).} \item{times}{The vector of time points at which AUC is
#' evaluated.} \item{iauc}{The summary measure of AUC.} \code{sens.uno} and
#' \code{spec.uno} return matrices of dimensions \code{times} x \code{(lpnew +
#' 1)}. The elements of these matrices are the sensitivity and specificity
#' estimates for each threshold of \code{lpnew} and for each time point
#' specified in \code{times}.
#' @seealso \code{\link{AUC.cd}}, \code{\link{AUC.sh}}, \code{\link{AUC.hc}},
#' \code{\link{IntAUC}}
#' @references
#' 
#' Uno, H., T. Cai, L. Tian, and L. J. Wei (2007).\cr Evaluating prediction
#' rules for t-year survivors with censored regression models.\cr \emph{Journal
#' of the American Statistical Association} \bold{102}, 527--537.\cr
#' @keywords classif
#' @examples
#' 
#' data(cancer,package="survival")
#' TR <- ovarian[1:16,]
#' TE <- ovarian[17:26,]
#' train.fit  <- survival::coxph(survival::Surv(futime, fustat) ~ age,
#'                     x=TRUE, y=TRUE, method="breslow", data=TR)
#' 
#' lpnew <- predict(train.fit, newdata=TE)
#' Surv.rsp <- survival::Surv(TR$futime, TR$fustat)
#' Surv.rsp.new <- survival::Surv(TE$futime, TE$fustat)
#' times <- seq(10, 1000, 10)                  
#' 
#' AUC_Uno <- AUC.uno(Surv.rsp, Surv.rsp.new, lpnew, times)
#' names(AUC_Uno)
#' AUC_Uno$iauc
NULL

##############################################################
## Uno sensetifity 
##############################################################
# Surv.rsp = Zielvariable train, Surv-Objekt (time, status)
# Surv.rsp.new = Zielvariable test, Surv-Objekt (time, status)
# lpnew = lin. Praedikt. aus Cox-Modell, z.B. predict(train.fit, newdata=test.data)
# times = Vektor der Zeitpunkte, an denen ausgewertet werden soll

#' @rdname survAUC_Uno
#' @export
sens.uno <- function(Surv.rsp, Surv.rsp.new, lpnew, times){
	thresh <- my.sort(unique(lpnew))
	n_th <- length(thresh)
	n_t <- length(times)
	ERG <- .C(`C_sens_uno`,
			  as.numeric(rep(1, n_t*(n_th+1))),
			  as.numeric(Surv.rsp[,1]),
			  as.numeric((1-Surv.rsp[,2])),
			  as.numeric(thresh),
			  as.numeric(times),
			  as.numeric(lpnew),
			  as.numeric(Surv.rsp.new[,1]),
			  as.numeric(Surv.rsp.new[,2]),
			  as.integer(n_th),
			  as.integer(n_t),
			  as.integer(dim(Surv.rsp.new)[1]),
			  as.integer(dim(Surv.rsp)[1]))
	#No longer needed since the symbol is registered in the NAMESPACE
	#          ,PACKAGE="survAUC")
	matrix(ERG[[1]], n_t, n_th+1)
}


##############################################################
## Uno specificity
##############################################################
# Surv.rsp.new = Zielvariable test, Surv-Objekt (time, status)
# lpnew = lin. Praedikt. aus Cox-Modell, z.B. predict(train.fit, newdata=test.data)
# times = Vektor der Zeitpunkte, an denen ausgewertet werden soll

#' @rdname survAUC_Uno
#' @export
spec.uno <- function(Surv.rsp.new, lpnew, times){
	thresh <- my.sort(unique(lpnew))
	n_th <- length(thresh)
	n_t <- length(times)
	ERG <- .C(`C_spec_uno`, 
			  as.numeric(rep(0, n_t*(n_th+1))), 
			  as.numeric(thresh), 
			  as.numeric(times),
			  as.numeric(lpnew), 
			  as.numeric(Surv.rsp.new[,1]), 
			  as.integer(n_th),
			  as.integer(n_t), 
			  as.integer(dim(Surv.rsp.new)[1]))
	#No longer needed since the symbol is registered in the NAMESPACE
	#          ,PACKAGE="survAUC")
	matrix(ERG[[1]], n_t, n_th+1)
}



##############################################################
## Uno AUC
##############################################################
# Surv.rsp = Zielvariable train, Surv-Objekt (time, status)
# Surv.rsp.new = Zielvariable test, Surv-Objekt (time, status)
# lpnew = lin. Praediktoren aus Cox-Modell, z.B. predict(train.fit, newdata=test.data)
# times = Vektor der Zeitpunkte, an denen ausgewertet werden soll
# weight = Welche Gewichtung der Integrated AUC?; rescale oder conditional.

#' @rdname survAUC_Uno
#' @export
AUC.uno <- function(Surv.rsp, Surv.rsp.new, lpnew, times, savesensspec=FALSE){

	thresh <- my.sort(unique(lpnew))
	n_th <- length(thresh)
	n_t <- length(times)
	
	#### Sensetivity, Specificity and AUC.
	auc.uno <- .C(`C_auc_uno`,
				  as.numeric(vector("numeric",length=n_t)),
				  as.numeric(0),
				  as.numeric(vector("numeric",length=n_t*(n_th+1))+1),
				  as.numeric(vector("numeric",length=n_t*(n_th+1))),
				  as.numeric(Surv.rsp[,1]),
				  as.numeric(1-Surv.rsp[,2]),
				  as.numeric(thresh), 
				  as.numeric(times),
				  as.numeric(lpnew), 
				  as.numeric(Surv.rsp.new[,1]),
				  as.numeric(Surv.rsp.new[,2]),
				  as.integer(n_th),
				  as.integer(n_t), 
				  as.integer(dim(Surv.rsp.new)[1]),
				  as.integer(dim(Surv.rsp)[1]))
	#No longer needed since the symbol is registered in the NAMESPACE
	#          ,PACKAGE="survAUC")
	if(!savesensspec){
		erg <- list(auc=auc.uno[[1]], times=auc.uno[[8]], iauc=auc.uno[[2]])
	}else{
		erg <- list(auc=auc.uno[[1]], times=auc.uno[[8]], iauc=auc.uno[[2]],
			 sens=matrix(auc.uno[[3]], n_t, n_th+1), 
			 spec=matrix(auc.uno[[4]], n_t, n_th+1),
			 thresh=auc.uno[[7]])
	}
	class(erg) <- "survAUC"
	erg
}



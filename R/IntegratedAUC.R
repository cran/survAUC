##############################################################
## Integrated AUC
##############################################################
# AUC	- vector of AUC's
# times	- vector of times
# S		- vector of survival probability
# tmax	- maximum timepoint



#' @title Integration of time-dependent AUC curves
#' 
#' @description Compute summary measures of a time-dependent AUC curve
#' 
#' @details This function calculates the integral under a time-dependent AUC curve
#' (\dQuote{IAUC} measure) using the integration limits [0, \code{tmax}]. The
#' values of the AUC curve are specified via the \code{AUC} argument.
#' 
#' In case \code{auc.type = "cumulative"} (cumulative/dynamic IAUC), the values
#' of \code{AUC} are weighted by the estimated probability density of the
#' time-to-event outcome. In case \code{auc.type = "incident"}
#' (incident/dynamic IAUC), the values of \code{AUC} are weighted by 2 times
#' the product of the estimated probability density and the (estimated)
#' survival function of the time-to-event outcome.  The survival function has
#' to be specified via the \code{S} argument.
#' 
#' As shown by Heagerty and Zheng (2005), the incident/dynamic version of IAUC
#' can be interpreted as a global concordance index measuring the probability
#' that observations with a large predictor value have a shorter survival time
#' than observations with a small predictor value. The incident/dynamic version
#' of IAUC has the same interpretation as Harrell's C for survival data.
#' 
#' @param AUC A vector of AUCs.
#' @param times The vector of time points corresponding to \code{AUC}.
#' @param S A vector of survival probabilities corresponding to \code{times}.
#' @param tmax A number specifying the upper limit of the time range for which
#' to compute the summary measure.
#' @param auc.type A string defining the type of AUC. 'cumulative' refers to
#' cumulative/dynamic AUC, 'incident' refers to incident/dynamic AUC.
#' @return A scalar number corresponding to the summary measure of interest.
#' @seealso \code{\link{AUC.cd}}, \code{\link{AUC.sh}}, \code{\link{AUC.uno}},
#' \code{\link{AUC.hc}}
#' @references
#' 
#' Harrell, F. E., R. M. Califf, D. B. Pryor, K. L. Lee and R. A. Rosati
#' (1982). \cr Evaluating the yield of medical tests.\cr \emph{Journal of the
#' American Medical Association} \bold{247}, 2543--2546.\cr
#' 
#' Harrell, F. E., K. L. Lee, R. M. Califf, D. B. Pryor and R. A. Rosati
#' (1984). \cr Regression modeling strategies for improved prognostic
#' prediction.\cr \emph{Statistics in Medicine} \bold{3}, 143--152.\cr
#' 
#' Heagerty, P. J. and Y. Zheng (2005). \cr Survival model predictive accuracy
#' and ROC curves.\cr \emph{Biometrics} \bold{61}, 92--105.\cr
#' @keywords classif manip
#' @examples
#' 
#' data(cancer,package="survival")
#' TR <- ovarian[1:16,]
#' TE <- ovarian[17:26,]
#' train.fit  <- survival::coxph(survival::Surv(futime, fustat) ~ age,
#'                               x=TRUE, y=TRUE, method="breslow", data=TR)
#' 
#' lp <- predict(train.fit)
#' lpnew <- predict(train.fit, newdata=TE)
#' Surv.rsp <- survival::Surv(TR$futime, TR$fustat)
#' Surv.rsp.new <- survival::Surv(TE$futime, TE$fustat)
#' times <- seq(10, 1000, 10)                  
#' 
#' 
#' AUC_CD <- AUC.cd(Surv.rsp, Surv.rsp.new, lp, lpnew, times)
#' IntAUC(AUC_CD$auc, AUC_CD$times, runif(length(times),0,1), median(times), auc.type="cumulative")
#' 
#' @export IntAUC
IntAUC <- function(AUC, times, S, tmax, auc.type="cumulative")
{
	n_S <- length(S)
	n_AUC <- length(AUC)
	n_times <- length(times)
	if(!((n_S == n_AUC) && (n_AUC == n_times)))
		stop("AUC, times and S must be the same length!")
	auc.type <- charmatch( auc.type, c("cumulative","incident") )
	if (is.na(auc.type))
		stop("auc.type must be one of 'cumulative' or 'incident'")
	maxI <- sum( times <= tmax )
	ind_S <- S[min(maxI+1,length(S))]
	iAUC <- .C(`C_int_auc`,
			   as.numeric(AUC),
			   as.numeric(times),
			   as.numeric(S),
			   as.numeric(tmax),
			   as.integer(n_S),
			   as.integer(maxI),
			   as.numeric(ind_S),
			   as.integer(auc.type-1),
			   as.numeric(0))
	#No longer needed since the symbol is registered in the NAMESPACE
	#          ,PACKAGE="survAUC")
	iAUC[[9]]
}


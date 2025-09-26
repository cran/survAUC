##############################################################
## C-Statistic suggest by Begg
##############################################################
## Surv.rsp		- the Surv(.,.) Outcome of training data
## Surv.rsp.new	- the Surv(.,.) Outcome of test data
## lp			- the vector of linear predictors of training data
## lpnew		- the vector of linear predictors of test data






#' @title C-statistic by Begg et al.
#' 
#' @description C-statistic by Begg et al.
#' 
#' @details This function implements the C-statistic proposed by Begg et al. (2000). It
#' has the same interpretation as Harrell's C for survival data (implemented in
#' the \code{rcorr.cens} function of the \bold{Hmisc} package).  \code{BeggC}
#' is restricted to Cox regression.  Specifically, it is assumed that \code{lp}
#' and \code{lpnew} are the predictors of a Cox proportional hazards model.
#' Estimates obtained from \code{BeggC} are valid as long as the Cox model is
#' specified correctly.
#' 
#' @param Surv.rsp A \code{Surv(.,.)} object containing to the outcome of the
#' training data.
#' @param Surv.rsp.new A \code{Surv(.,.)} object containing the outcome of the
#' test data.
#' @param lp The vector of predictors estimated from the training data.
#' @param lpnew The vector of predictors obtained from the test data.
#' @return The estimated C-statistic.
#' @seealso \code{\link{UnoC}}, \code{\link{GHCI}}, \code{\link{AUC.sh}},
#' \code{\link{IntAUC}}
#' @references
#' 
#' Begg, B. C., L. D. Craemer, E. S. Venkatraman and J. Rosai (2000). \cr
#' Comparing tumor staging and grading systems: a case study and a review of
#' the issues, using thymoma as a model.\cr \emph{Statistics in Medicine}
#' \bold{19}, 1997--2014.\cr
#' @keywords classif
#' @examples
#' 
#' data(cancer,package="survival")
#' TR <- ovarian[1:16,]
#' TE <- ovarian[17:26,]
#' train.fit  <- survival::coxph(survival::Surv(futime, fustat) ~ age,
#'                     x=TRUE, y=TRUE, method="breslow", data=TR)
#' 
#' lp <- predict(train.fit)
#' lpnew <- predict(train.fit, newdata=TE)
#' Surv.rsp <- survival::Surv(TR$futime, TR$fustat)
#' Surv.rsp.new <- survival::Surv(TE$futime, TE$fustat)              
#' 
#' Cstat <- BeggC(Surv.rsp, Surv.rsp.new, lp, lpnew)
#' Cstat
#' 
#' 
#' @export BeggC
BeggC <- function(Surv.rsp, Surv.rsp.new, lp, lpnew){

## Surv-train
	stime <- Surv.rsp[,1]
	event <- Surv.rsp[,2]
	
## Surv-test
	stime.new <- Surv.rsp.new[,1]
	event.new <- Surv.rsp.new[,2]
## Times
	times <- stime.new
	n_times <- length(times)
	n_stime <- length(stime)
	n_lp <- length(lp)
	n_stime_new <- length(stime.new)
	n_lpnew <- length(lpnew)
	if(n_stime != n_lp)
		stop(" 'Surv.rsp' and 'linear predictors' must have the same length!\n")
	if(n_stime_new != n_lpnew)
		stop(" 'Surv.rsp.new' and 'linear predictors new' must have the same length!\n")

	#### Cox survival function estimates for lpnew	
	surv.cox <- .Call(`C_survfit_cox`,
					  as.numeric(lp), 
					  as.numeric(stime),
					  as.numeric(event), 
					  as.integer(n_stime),
					  as.integer(n_lp),
					  as.numeric(lpnew),
					  as.integer(n_lpnew))
	#No longer needed since the symbol is registered in the NAMESPACE
	#          ,PACKAGE="survAUC")
	#### C-Statistic
	c.begg <- .C(`C_begg`,
				   as.numeric(stime.new), 
				   as.numeric(event.new),
				   as.integer(n_stime_new),
				   as.numeric(times),
				   as.integer(n_times),
				   as.numeric(lp),
				   as.numeric(lpnew),
				   as.numeric(surv.cox[[1]]),
				   as.numeric(surv.cox[[2]]),
				   as.integer(length(surv.cox[[2]])),
				   as.numeric(vector("numeric",length=1)))
				   #No longer needed since the symbol is registered in the NAMESPACE
				   #          ,PACKAGE="survAUC")
	c.begg[[11]]
}



################################################################
###						Hung & Chiang						####
################################################################
## Surv.rsp		- the Surv(.,.) Outcome of training data
## Surv.rsp.new	- the Surv(.,.) Outcome of test data
## lpnew		- the vector of linear predictors of test data
## times		- the vector of times



#' @title AUC estimator proposed by Hung and Chiang
#' 
#' @rdname survAUC_HungChiang
#' @name survAUC_HungChiang
#' 
#' @description Hung and Chiang's estimator of cumulative/dynamic AUC for right-censored
#' time-to-event data
#' 
#' @details This function implements the estimator of cumulative/dynamic AUC proposed by
#' Hung and Chiang (2010). The estimator is based on
#' inverse-probability-of-censoring weights and does not assume a specific
#' working model for deriving the predictor \code{lpnew}. It is assumed,
#' however, that there is a one-to-one relationship between the predictor and
#' the expected survival times conditional on the predictor. The \code{iauc}
#' summary measure is given by the integral of AUC on [0, max(\code{times})]
#' (weighted by the estimated probability density of the time-to-event
#' outcome).
#' 
#' Note that the estimator implemented in \code{AUC.hc} is restricted to
#' situations where the random censoring assumption holds (formula (4) in Hung
#' and Chiang 2010).
#' 
#' @param Surv.rsp A \code{Surv(.,.)} object containing to the outcome of the
#' training data.
#' @param Surv.rsp.new A \code{Surv(.,.)} object containing the outcome of the
#' test data.
#' @param lpnew The vector of predictors obtained from the test data.
#' @param times A vector of time points at which to evaluate AUC.
#' @return \code{AUC.hc} returns an object of class \code{survAUC}.
#' Specifically, \code{AUC.hc} returns a list with the following components:
#' \item{auc}{The cumulative/dynamic AUC estimates (evaluated at
#' \code{times}).} \item{times}{The vector of time points at which AUC is
#' evaluated.} \item{iauc}{The summary measure of AUC.}
#' @seealso \code{\link{AUC.uno}}, \code{\link{AUC.sh}}, \code{\link{AUC.cd}},
#' \code{\link{IntAUC}}
#' @references
#' 
#' Hung, H. and C.-T. Chiang (2010). \cr Estimation methods for time-dependent
#' AUC models with survival data.\cr \emph{Canadian Journal of Statistics}
#' \bold{38}, 8--26.\cr
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
#' AUC_hc <- AUC.hc(Surv.rsp, Surv.rsp.new, lpnew, times)
#' AUC_hc
#' 
#' 
#' @export AUC.hc
AUC.hc <- function(Surv.rsp, Surv.rsp.new, lpnew, times)
{
## Surv-train
	stime <- Surv.rsp[,1]
	event <- Surv.rsp[,2]
	
## Surv-test
	stime.new <- Surv.rsp.new[,1]
	event.new <- Surv.rsp.new[,2]

	n_time <- length(times)
	n_stime <- length(stime)
	n_stime_new <- length(stime.new)
	n_lpnew <- length(lpnew)
	auc <- vector("numeric",n_time)
	
	ans <- .C(`C_Hung_Chiang`,
			  as.numeric(times),
			  as.integer(n_time),
			  as.numeric(stime),
			  as.numeric(event),
			  as.integer(n_stime),
			  as.numeric(stime.new),
			  as.numeric(event.new),
			  as.integer(n_stime_new),
			  as.numeric(lpnew),
			  as.integer(n_lpnew),
			  as.numeric(auc),
			  as.numeric(0))
	#No longer needed since the symbol is registered in the NAMESPACE
	#          ,PACKAGE="survAUC")
	erg <- list(auc=ans[[11]], times=times, iauc=ans[[12]])
	class(erg) <- "survAUC"
	erg
}



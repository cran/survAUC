################################################################
###						Gonen and Hellers					 ###
###              Concordance Index for Cox models			 ###
################################################################
## lpnew			- the vector of linear predictors of data




#' @title Gonen and Heller's Concordance Index for Cox models
#' 
#' @description Gonen and Heller's Concordance Index for Cox proportional hazards models
#' 
#' @details This function implements the concordance probability estimator proposed by
#' Gonen and Heller (2005). It has the same interpretation as Harrell's C for
#' survival data (implemented in the \code{rcorr.cens} function of the
#' \bold{Hmisc} package).
#' 
#' The results obtained from \code{GHCI} are valid as long as \code{lpnew} is
#' the predictor of a correctly specified Cox proportional hazards model. In
#' this case, the estimator remains valid even if the censoring times depend on
#' the values of the predictor.
#' 
#' Note that the smoothed version of \code{GHCI}, which is proposed in Section
#' 3 of Gonen and Heller (2005), is not implemented in R package
#' \bold{survAUC}.
#' 
#' @param lpnew The vector of predictors obtained from the test data.
#' @return A length-one numeric vector containing the concordance probability
#' estimate.
#' @seealso \code{\link{AUC.sh}}, \code{\link{IntAUC}}
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
#' Gonen, M. and G. Heller (2005). \cr Concordance probability and
#' discriminatory power in proportional hazards regression.\cr
#' \emph{Biometrika} \bold{92}, 965--970.\cr
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
#'                  
#' GHCI(lpnew)
#' 
#' 
#' @export GHCI
GHCI <- function(lpnew){
	ans <- .C(`C_GHCI`,
			  as.numeric(lpnew),
			  as.integer(length(lpnew)),
			  as.numeric(0.0))
	#No longer needed since the symbol is registered in the NAMESPACE
	#          ,PACKAGE="survAUC")
	ans[[3]]
}






################################################################
###						Prediction Error					 ###
###						robust & brier						 ###
################################################################
## Surv.rsp		- Surv(.,.) Outcome of training data
## Surv.rsp.new	- Surv(.,.) Outcome of test data
## lp			- vector of linear predictors of training data
## lpnew		- vector of linear predictors of test data
## times		- vector of times
## type			- kind of prediction error curve: 'brier' or 'robust'



#' @title Distance-based estimators of survival predictive accuracy
#' 
#' @description Inverse-probability-of-censoring weighted estimators of absolute and squared
#' deviations between survival functions
#' 
#' @details This function implements two types of prediction error curves for
#' right-censored time-to-event data: The Brier Score (\code{type = "brier"},
#' Gerds and Schumacher 2006) estimates the \emph{squared} deviation between
#' predicted and observed survival whereas the method proposed by Schmid et al.
#' (2011) estimates the \emph{absolute} deviation between predicted and
#' observed survival (\code{type = "robust"}).
#' 
#' Both methods are based on inverse-probability-of-censoring weights and do
#' not assume a specific working model for survival prediction.  Note, however,
#' that the estimators implemented in \code{predErr}, are restricted to
#' situations where the random censoring assumption holds.
#' 
#' Time-independent summary measures of prediction error are given by the the
#' areas under the prediction error curves. If \code{int.type = "weighted"},
#' prediction errors are weighted by the estimated probability density of the
#' time-to-event outcome.
#' 
#' @param Surv.rsp A \code{Surv(.,.)} object containing to the outcome of the
#' training data.
#' @param Surv.rsp.new A \code{Surv(.,.)} object containing the outcome of the
#' test data.
#' @param lp The vector of predictors estimated from the training data.
#' @param lpnew The vector of predictors obtained from the test data.
#' @param times A vector of time points at which to evaluate the prediction
#' error curve.
#' @param type A string specifying the type of prediction error curve: 'brier'
#' refers to the squared deviation between predicted and observed survival
#' (Brier score), 'robust' refers to the absolute deviation between predicted
#' and observed survival.
#' @param int.type A string specifying the type of integration method for the
#' prediction error curves. Either 'unweighted' or 'weighted'.
#' @return \code{predErr} returns an object of class \code{survErr}.
#' Specifically, \code{predErr} returns a list containing the following
#' components: \item{error}{The prediction error estimates (evaluated at
#' \code{times}).} \item{times}{The vector of time points at which prediction
#' errors are evaluated.} \item{ierror}{The integrated prediction error.}
#' @seealso \code{\link{IntAUC}}, \code{\link{OXS}}, \code{\link{schemper}}
#' @references
#' 
#' Gerds, T. A. and M. Schumacher (2006).\cr Consistent estimation of the
#' expected Brier score in general survival models with right-censored event
#' times.\cr \emph{Biometrical Journal} \bold{48}, 1029--1040.\cr
#' 
#' Schmid, M., T. Hielscher, T. Augustin, and O. Gefeller (2011).\cr A robust
#' alter- native to the Schemper-Henderson estimator of prediction error.\cr
#' \emph{Biometrics} \bold{67}, 524--535.\cr
#' @keywords classif
#' @examples
#' 
#' data(cancer,package="survival")
#' TR <- ovarian[1:16,]
#' TE <- ovarian[17:26,]
#' train.fit  <- survival::coxph(survival::Surv(futime, fustat) ~ age, x=TRUE, y=TRUE, 
#'                     method="breslow", data=TR)
#' 
#' lp <- predict(train.fit)
#' lpnew <- predict(train.fit, newdata=TE)
#' Surv.rsp <- survival::Surv(TR$futime, TR$fustat)
#' Surv.rsp.new <- survival::Surv(TE$futime, TE$fustat)
#' times <- 1:500                  
#' 
#' predErr(Surv.rsp, Surv.rsp.new, lp, lpnew, times, 
#'         type = "brier", int.type = "unweighted")
#' 
#' predErr(Surv.rsp, Surv.rsp.new, lp, lpnew, times, 
#'         type = "robust", int.type = "unweighted")
#' 
#' predErr(Surv.rsp, Surv.rsp.new, lp, lpnew, times, 
#'         type = "brier", int.type = "weighted")
#' 
#' 
#' @export predErr
predErr <- function(Surv.rsp, Surv.rsp.new, lp, lpnew, times, 
					type = "brier", int.type = "unweighted")
{
	type <- charmatch( type, c("brier","robust") )
	if (is.na(type))
		stop("'type' must be one of 'brier' or 'robust'")
	int.type <- charmatch( int.type, c("weighted","unweighted") )
	if (is.na(int.type))
		stop("'int.type' must be one of 'weighted' or 'unweighted'")

	## Surv-train
	stime <- Surv.rsp[,1]
	event <- Surv.rsp[,2]
	
	## Surv-test
	stime.new <- Surv.rsp.new[,1]
	event.new <- Surv.rsp.new[,2]

	n.times <- length(times)
	n.stime <- length(stime)
	n.stime.new <- length(stime.new)
	n.lp <- length(lp)
	n.lpnew <- length(lpnew)
	
	erg <- .Call(`C_predError`,
				 as.numeric(stime),
				 as.numeric(event),
				 as.integer(n.stime),
				 as.numeric(stime.new),
				 as.numeric(event.new),
				 as.integer(n.stime.new),
				 as.numeric(times),
				 as.integer(length(times)),
				 as.numeric(lp),
				 as.integer(n.lp),
				 as.numeric(lpnew),
				 as.integer(n.lpnew),
				 as.integer(type-1),
				 as.integer(int.type-1))
	#No longer needed since the symbol is registered in the NAMESPACE
	#          ,PACKAGE="survAUC")
	class(erg) <- "survErr"
	erg
}






#' @title R2-type Coefficients for Cox proportional hazards models
#' 
#' @description R2-type Coefficients for Cox proportional hazards models
#' 
#' @name R2_type_Coef
#' @rdname R2_type_Coef
#' @details The \code{OXS}, \code{Nagelk} and \code{XO} functions implement three types
#' of R2 coefficients for right-censored time-to-event data: (a) The
#' coefficient proposed by O'Quigley et al. (2005) (\code{OXS}), (b) the
#' coefficient proposed by Nagelkerke (1991) (\code{Nagelk}) and (c) the
#' coefficient proposed by Xu and O'Quigley (1999) (\code{XO}).
#' 
#' Because the \code{OXS}, \code{Nagelk} and \code{XO} functions assume that
#' \code{lp} and \code{lpnew} were derived from a correctly specified Cox
#' proportional hazards model, estimates obtained from these functions are only
#' valid if the Cox model holds.
#' 
#' @aliases OXS Nagelk XO
#' @param Surv.rsp A \code{Surv(.,.)} object containing to the outcome of the
#' test data.
#' @param lp The vector of predictors.
#' @param lp0 The vector of predictors obtained from the covariate-free null
#' model.
#' @return The estimated R2 coefficient.
#' @seealso \code{\link{predErr}}, \code{\link{schemper}}, \code{\link{GHCI}}
#' @references
#' 
#' Nagelkerke, N. J. D. (1991).\cr A note on a general definition of the
#' coefficient of determination.\cr \emph{Biometrika} \bold{78}, 691--692.\cr
#' 
#' O'Quigley, J., R. Xu, and J. Stare (2005).\cr Explained randomness in
#' proportional hazards models.\cr \emph{Statistics in Medicine} \bold{24},
#' 479--489.\cr
#' 
#' Xu, R. and J. O'Quigley (1999).\cr A measure of dependence for proportional
#' hazards models.\cr \emph{Journal of Nonparametric Statistics} \bold{12},
#' 83--107.\cr
#' @keywords classif
#' @examples
#' 
#' data(cancer,package="survival")
#' TR <- ovarian[1:16,]
#' TE <- ovarian[17:26,]
#' train.fit  <- survival::coxph(survival::Surv(futime, fustat) ~ age,
#'                     x=TRUE, y=TRUE, method="breslow", data=TR)
#' 
#' model0 <- survival::coxph(survival::Surv(futime, fustat)~1, data=TR)
#' model1 <- survival::coxph(survival::Surv(futime, fustat)~age, data=TR)
#' f0 <- rep(0,nrow(TE))
#' f1 <- predict(model1, newdata=TE)               
#' Surv.res <- survival::Surv(TE$futime, TE$fustat)
#' 
#' OXS(Surv.res, f1, f0)
#' Nagelk(Surv.res, f1, f0)
#' XO(Surv.res, f1, f0)
#' 
NULL

################################################################
###				measure by O''Quigley et al. (2005)			 ###
###						  R^2_{OXS}  						 ###
################################################################
## Surv.rsp		- Surv(.,.) Outcome of training data
## lp			- vector of linear predictors
## lp0			- vector of linear predictors of null-model

#' @rdname R2_type_Coef
#' @export
OXS <- function(Surv.rsp, lp, lp0)
{
	
	L <- PartLCox(Surv.rsp, lp)
	L0 <- PartLCox(Surv.rsp, lp0)	
	1 - exp( - 2 * (L-L0) / sum(Surv.rsp[,2]))
}






################################################################
###				measure by Nagelkerke						 ###
###						  R^2_{N}							 ###
################################################################
## Surv.rsp		- Surv(.,.) Outcome of training data
## lp			- vector of linear predictors
## lp0			- vector of linear predictors of null-model

#' @rdname R2_type_Coef
#' @export
Nagelk <- function(Surv.rsp, lp, lp0)
{	
	L <- PartLCox(Surv.rsp, lp)
	L0 <- PartLCox(Surv.rsp, lp0)
	n <- length(lp)
	(1 - exp( - 2 * (L-L0) / n)) / (1 - exp( 2 * L0 / n))
}





################################################################
###				 measure by Xu & O''Quigley					 ###
###						  R^2_{XO}							 ###
################################################################
## Surv.rsp		- Surv(.,.) Outcome of training data
## lp			- vector of linear predictors
## lp0			- vector of linear predictors of null-model

#' @rdname R2_type_Coef
#' @export
XO <- function(Surv.rsp, lp, lp0){

	time <- Surv.rsp[,1]
	event <- Surv.rsp[,2]
	n <- length(time)
	n_lp <- length(lp)
	n_lp0 <- length(lp0)
	if(n != n_lp || n_lp != n_lp0 || n != n_lp0)
		stop(" 'Surv.rsp', 'linear predictors' and 'linear predictors of null-model' must have the same length!\n")
	ans <- .C(`C_XO`,
			  as.numeric(time),
			  as.numeric(event),
			  as.integer(n), 
			  as.numeric(lp),
			  as.numeric(lp0),
			  as.numeric(0))
	#No longer needed since the symbol is registered in the NAMESPACE
	#          ,PACKAGE="survAUC")
	ans[[6]]
}





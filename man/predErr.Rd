\name{predErr}
\alias{predErr}

\title{Distance-based estimators of survival predictive accuracy}
\description{
Inverse-probability-of-censoring weighted estimators of absolute and squared
deviations between survival functions
}
\usage{
predErr(Surv.rsp, Surv.rsp.new, lp, lpnew, times, 
        type = "brier", int.type = "unweighted")
}
\arguments{
  \item{Surv.rsp}{A \code{Surv(.,.)} object containing to the outcome of the training data.}
  \item{Surv.rsp.new}{A \code{Surv(.,.)} object containing the outcome of the test data.}
  \item{lp}{The vector of predictors estimated from the training data.}
  \item{lpnew}{The vector of predictors obtained from the test data.}
  \item{times}{A vector of time points at which to evaluate the prediction error curve.}
  \item{type}{A string specifying the type of prediction error curve: 'brier' refers
              to the squared deviation between predicted and observed survival (Brier score), 'robust'
              refers to the absolute deviation between predicted and observed survival.}
  \item{int.type}{A string specifying the type of integration method for the prediction
                  error curves. Either 'unweighted' or 'weighted'.}  
}
\details{
  This function implements two types of prediction error curves for right-censored
  time-to-event data: The Brier Score (\code{type = "brier"}, Gerds and Schumacher 2006) estimates the \emph{squared}
  deviation between predicted and observed survival whereas the method proposed
  by Schmid et al. (2011) estimates the \emph{absolute} deviation between predicted and
  observed survival (\code{type = "robust"}).
  
  Both methods are based on inverse-probability-of-censoring
  weights and do not assume a specific working model for survival prediction.
  Note, however, that the estimators implemented in \code{predErr},
  are restricted to situations where the random censoring assumption holds.
  
  Time-independent summary measures of prediction error are given by the
  the areas under the prediction error curves. If \code{int.type = "weighted"},
  prediction errors are weighted by the estimated probability density of
  the time-to-event outcome.

}
\references{

Gerds, T. A. and M. Schumacher (2006).\cr Consistent estimation of the expected
Brier score in general survival models with right-censored event times.\cr
\emph{Biometrical Journal} \bold{48}, 1029--1040.\cr

Schmid, M., T. Hielscher, T. Augustin, and O. Gefeller (2011).\cr A robust alter-
native to the Schemper-Henderson estimator of prediction error.\cr \emph{Biometrics}
\bold{67}, 524--535.\cr

}
\value{
  \code{predErr} returns an object of class \code{survErr}. Specifically,
  \code{predErr} returns a list containing the following components:
  \item{error}{The prediction error estimates (evaluated at \code{times}).}
  \item{times}{The vector of time points at which prediction errors are evaluated.}
  \item{ierror}{The integrated prediction error.}
}

\seealso{
  \code{\link{IntAUC}}, \code{\link{OXS}}, \code{\link{schemper}}
}
\examples{

TR <- ovarian[1:16,]
TE <- ovarian[17:26,]
train.fit  <- coxph(Surv(futime, fustat) ~ age, x=TRUE, y=TRUE, 
                    method="breslow", data=TR)

lp <- predict(train.fit)
lpnew <- predict(train.fit, newdata=TE)
Surv.rsp <- Surv(TR$futime, TR$fustat)
Surv.rsp.new <- Surv(TE$futime, TE$fustat)
times <- 1:500                  

predErr(Surv.rsp, Surv.rsp.new, lp, lpnew, times, 
        type = "brier", int.type = "unweighted")

predErr(Surv.rsp, Surv.rsp.new, lp, lpnew, times, 
        type = "robust", int.type = "unweighted")

predErr(Surv.rsp, Surv.rsp.new, lp, lpnew, times, 
        type = "brier", int.type = "weighted")

}

\keyword{classif}


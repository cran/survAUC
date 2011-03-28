\name{predErr}
\alias{predErr}

\title{Distance-based Estimators of Survival Predictive Accuracy}
\description{
Inverse-probability-of-censoring Weighted Estimators of Absolute and Squared
Deviations between Survival Functions (Prediction Error Curves)
}
\usage{
predErr(Surv.rsp, Surv.rsp.new, lp, lpnew, times, 
					type = "brier", int.type = "unweighted")
}
\arguments{
  \item{Surv.rsp}{A \code{Surv(.,.)} object containing to the outcome of the training data.}
  \item{Surv.rsp.new}{A \code{Surv(.,.)} object containing the outcome of the test data.}
  \item{lp}{The vector of linear predictors estimated from the training data.}
  \item{lpnew}{The vector of linear predictors obtained from the test data.}
  \item{times}{A vector of time points at which to compute the prediction error curve.}
  \item{type}{A string specifying the type of prediction error curve: 'brier' refers
              to the squared deviation between predicted and observed survival (Brier score), 'robust'
              refers to the absolute deviation between predicted and observed survival.}
  \item{int.type}{A string specifying the type of the integration method for prediction
                  error curves. Either 'unweighted' or 'weighted'.}  
}
\references{

Gerds, T. A. and M. Schumacher (2006).\cr Consistent estimation of the expected
Brier score in general survival models with right-censored event times.\cr
\emph{Biometrical Journal} \bold{48}, 1029--1040.\cr

Schmid, M., T. Hielscher, T. Augustin, and O. Gefeller (2011).\cr A robust alter-
native to the Schemper-Henderson estimator of prediction error.\cr \emph{Biometrics}.
To appear.\cr

}

\seealso{
  \code{\link{IntAUC}}, \code{\link{OXS}}, \code{\link{schemper}}
}
\examples{

TR <- ovarian[1:16,]
TE <- ovarian[17:26,]
train.fit  <- coxph(Surv(futime, fustat) ~ age,
                    x=TRUE, y=TRUE, method="breslow", data=TR)

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


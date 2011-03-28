\name{AUC.cd}
\alias{AUC.cd}

\title{AUC estimation proposed by Chambless and Diao}
\description{
  Chambless and Diao's estimator of cumulative/dynamic AUC for right-censored
  time-to-event data
}
\usage{
AUC.cd(Surv.rsp, Surv.rsp.new = NULL, lp, lpnew, times)
}
\arguments{
  \item{Surv.rsp}{A \code{Surv(.,.)} object containing to the outcome of the training data.}
  \item{Surv.rsp.new}{A \code{Surv(.,.)} object containing the outcome of the test data.}
  \item{lp}{The vector of linear predictors estimated from the training data.}
  \item{lpnew}{The vector of linear predictors obtained from the test data.}
  \item{times}{A vector of time points at which to compute the AUC.}
}
\references{

Chambless, L. E. and G. Diao (2006). \cr Estimation of time-dependent area under
the ROC curve for long-term risk prediction.\cr \emph{Statistics in Medicine}
\bold{25}, 3474--3486.\cr

}

\seealso{
  \code{\link{AUC.uno}}, \code{\link{AUC.sh}}, \code{\link{AUC.hc}},
  \code{\link{IntAUC}}
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
times <- seq(10, 1000, 10)                  

AUC_CD <- AUC.cd(Surv.rsp, Surv.rsp.new, lp, lpnew, times)

}

\keyword{classif}
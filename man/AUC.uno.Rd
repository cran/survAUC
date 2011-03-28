\name{AUC.uno}
\alias{AUC.uno}
\alias{spec.uno}
\alias{sens.uno}
\title{AUC estimation proposed by Uno et al.}
\description{
Uno's estimator of cumulative/dynamic AUC for right-censored time-to-event data
}
\usage{
AUC.uno(Surv.rsp, Surv.rsp.new, lpnew, times, savesensspec=FALSE)
sens.uno(Surv.rsp, Surv.rsp.new, lpnew, times)
spec.uno(Surv.rsp.new, lpnew, times)
}
\arguments{

  \item{Surv.rsp}{A \code{Surv(.,.)} object containing to the outcome of the training data.}
  \item{Surv.rsp.new}{A \code{Surv(.,.)} object containing the outcome of the test data.}
  \item{lpnew}{The vector of linear predictors obtained from the test data.}
  \item{times}{A vector of time points at which to compute the AUC.}
  \item{savesensspec}{A logical specifying whether sensitivities and specificities
                      should be saved.}
}
\references{

Uno, H., T. Cai, L. Tian, and L. J. Wei (2007).\cr Evaluating prediction rules for
t-year survivors with censored regression models.\cr \emph{Journal of the American
Statistical Association} \bold{102}, 527--537.\cr

}

\seealso{
  \code{\link{AUC.cd}}, \code{\link{AUC.sh}},  \code{\link{AUC.hc}},
  \code{\link{IntAUC}}
}
\examples{

TR <- ovarian[1:16,]
TE <- ovarian[17:26,]
train.fit  <- coxph(Surv(futime, fustat) ~ age,
                    x=TRUE, y=TRUE, method="breslow", data=TR)

lpnew <- predict(train.fit, newdata=TE)
Surv.rsp <- Surv(TR$futime, TR$fustat)
Surv.rsp.new <- Surv(TE$futime, TE$fustat)
times <- seq(10, 1000, 10)                  

AUC_Uno <- AUC.uno(Surv.rsp, Surv.rsp.new, lpnew, times)
names(AUC_Uno)
AUC_Uno$iauc

}

\keyword{classif}

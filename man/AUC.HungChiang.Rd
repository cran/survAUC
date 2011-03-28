\name{AUC.hc}
\alias{AUC.hc}

\title{AUC estimation proposed by Hung and Chiang}
\description{
Hung and Chiang's estimator of cumulative/dynamic AUC for right-censored
time-to-event data
}
\usage{
AUC.hc(Surv.rsp, Surv.rsp.new, lpnew, times)
}
\arguments{
  \item{Surv.rsp}{A \code{Surv(.,.)} object containing to the outcome of the training data.}
  \item{Surv.rsp.new}{A \code{Surv(.,.)} object containing the outcome of the test data.}
  \item{lpnew}{The vector of linear predictors obtained from the test data.}
  \item{times}{A vector of time points at which to compute the AUC.}
}
\references{

Hung, H. and C.-T. Chiang (2010). \cr Estimation methods for time-dependent
AUC models with survival data.\cr \emph{Canadian Journal of Statistics}
\bold{38}, 8--26.\cr

}

\seealso{
  \code{\link{AUC.uno}}, \code{\link{AUC.sh}}, \code{\link{AUC.cd}},
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

AUC_hc <- AUC.hc(Surv.rsp, Surv.rsp.new, lpnew, times)
AUC_hc

}

\keyword{classif}
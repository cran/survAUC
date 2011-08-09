\name{AUC.hc}
\alias{AUC.hc}

\title{AUC estimator proposed by Hung and Chiang}
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
  \item{lpnew}{The vector of predictors obtained from the test data.}
  \item{times}{A vector of time points at which to evaluate AUC.}
}
\details{
  This function implements the estimator of cumulative/dynamic AUC proposed by
  Hung and Chiang (2010). The estimator is based on inverse-probability-of-censoring
  weights and does not assume a specific working model for deriving the predictor
  \code{lpnew}. It is assumed, however, that there is a one-to-one
  relationship between the predictor and the expected survival times conditional
  on the predictor. The \code{iauc} summary measure is given by the integral of AUC on
  [0, max(\code{times})] (weighted by the estimated probability density of
  the time-to-event outcome).

  Note that the estimator implemented in \code{AUC.hc} is restricted to situations
  where the random censoring assumption holds (formula (4) in Hung and Chiang 2010).
}
\value{
  \code{AUC.hc} returns an object of class \code{survAUC}. Specifically,
  \code{AUC.hc} returns a list with the following components:
  \item{auc}{The cumulative/dynamic AUC estimates (evaluated at \code{times}).}
  \item{times}{The vector of time points at which AUC is evaluated.}
  \item{iauc}{The summary measure of AUC.}
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
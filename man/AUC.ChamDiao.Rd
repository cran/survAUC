\name{AUC.cd}
\alias{AUC.cd}

\title{AUC estimator proposed by Chambless and Diao}
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
  \item{lp}{The vector of predictors estimated from the training data.}
  \item{lpnew}{The vector of predictors obtained from the test data.}
  \item{times}{A vector of time points at which to evaluate AUC.}
}
\details{
  This function implements the estimator of cumulative/dynamic AUC proposed in
  Section 3.3 of Chambless and Diao (2006). In contrast to the general form of
  Chambless and Diao's estimator, \code{AUC.cd} is restricted to Cox regression.
  Specifically, it is assumed that \code{lp} and \code{lpnew} are the predictors of
  a Cox proportional hazards model. Estimates obtained from
  \code{AUC.cd} are valid as long as the Cox model is specified correctly.
  The \code{iauc} summary measure is given by the integral of AUC on
  [0, max(\code{times})] (weighted by the estimated probability density of
  the time-to-event outcome).
  
  Note that the recursive estimators proposed in Sections 3.1 and 3.2 of 
  Chambless and Diao (2006) are not implemented in the \bold{survAUC} package.
}
\value{
  \code{AUC.cd} returns an object of class \code{survAUC}. Specifically,
  \code{AUC.cd} returns a list with the following components:
  \item{auc}{The cumulative/dynamic AUC estimates (evaluated at \code{times}).}
  \item{times}{The vector of time points at which AUC is evaluated.}
  \item{iauc}{The summary measure of AUC.}
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
data(cancer,package="survival")
TR <- ovarian[1:16,]
TE <- ovarian[17:26,]
train.fit  <- survival::coxph(survival::Surv(futime, fustat) ~ age,
                    x=TRUE, y=TRUE, method="breslow", data=TR)

lp <- predict(train.fit)
lpnew <- predict(train.fit, newdata=TE)
Surv.rsp <- survival::Surv(TR$futime, TR$fustat)
Surv.rsp.new <- survival::Surv(TE$futime, TE$fustat)
times <- seq(10, 1000, 10)                  

AUC_CD <- AUC.cd(Surv.rsp, Surv.rsp.new, lp, lpnew, times)

}

\keyword{classif}

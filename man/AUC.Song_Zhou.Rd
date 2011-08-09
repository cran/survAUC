\name{AUC.sh}
\alias{AUC.sh}
\alias{spec.sh}
\alias{sens.sh}

\title{AUC estimator proposed by Song and Zhou}
\description{
  Song and Zhou's estimators of AUC for right-censored time-to-event data
}
\usage{
AUC.sh(Surv.rsp, Surv.rsp.new=NULL, lp, lpnew, times, 
		type="incident", savesensspec=FALSE)
sens.sh(Surv.rsp, lp, lpnew, times, type="incident")
spec.sh(Surv.rsp, lp, lpnew, times)
}
\arguments{
  \item{Surv.rsp}{A \code{Surv(.,.)} object containing to the outcome of the training data.}
  \item{Surv.rsp.new}{A \code{Surv(.,.)} object containing the outcome of the test data.}
  \item{lp}{The vector of predictors estimated from the training data.}
  \item{lpnew}{The vector of predictors obtained from the test data.}
  \item{times}{A vector of time points at which to evaluate AUC.}
  \item{type}{A string defining the type of true positive rate (TPR):
              \code{"incident"} refers to incident TPR ,
              \code{"cumulative"} refers to cumulative TPR.}
  \item{savesensspec}{A logical specifying whether sensitivities and specificities
                      should be saved.}
}
\details{
  The \code{sens.sh} and \code{spec.sh} functions implement the estimators of
  time-dependent true and false positive rates proposed by Song and Zhou (2008).
  
  The \code{AUC.sh} function implements the estimators of cumulative/dynamic and
  incident/dynamic AUC proposed by Song and Zhou (2008). These estimators are given
  by the areas under the time-dependent ROC curves estimated by
  \code{sens.sh} and \code{spec.sh}. In case of cumulative/dynamic
  AUC, the \code{iauc} summary measure is given by the integral of AUC on
  [0, max(\code{times})] (weighted by the estimated probability density of
  the time-to-event outcome). In case of incident/dynamic AUC, \code{iauc} is
  given by the integral of AUC on [0, max(\code{times})] (weighted by 2 times
  the product of the estimated probability density and the estimated survival
  function of the time-to-event outcome).
  
  The results obtained from \code{spec.sh}, \code{spec.sh} and \code{AUC.sh}
  are valid as long as \code{lp} and \code{lpnew} are the predictors of
  a correctly specified Cox proportional hazards model. In this case, the
  estimators remain valid even if the censoring times depend on the values of
  the predictors.

}
\value{
  \code{AUC.sh} returns an object of class \code{survAUC}. Specifically,
  \code{AUC.sh} returns a list with the following components:
  \item{auc}{The cumulative/dynamic or incident/dynamic AUC estimates
  (evaluated at \code{times}).}
  \item{times}{The vector of time points at which AUC is evaluated.}
  \item{iauc}{The summary measure of AUC.}
  
  \code{sens.sh} and \code{spec.sh} return matrices of dimensions \code{times} x
  \code{lpnew + 1}. The elements of these matrices are the sensitivity and
  specificity estimates for each threshold of \code{lpnew} and for each time point
  specified in \code{times}.
}
\references{

Song, X. and X.-H. Zhou (2008).\cr A semiparametric approach for the covariate
specific ROC curve with survival outcome.\cr \emph{Statistica Sinica}
\bold{18}, 947--965.\cr

}

\seealso{
  \code{\link{AUC.uno}}, \code{\link{AUC.cd}}, \code{\link{AUC.hc}},
  \code{\link{GHCI}}, \code{\link{IntAUC}}
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

AUC_sh <- AUC.sh(Surv.rsp, Surv.rsp.new, lp, lpnew, times)
names(AUC_sh)
AUC_sh$iauc

}

\keyword{classif}


\name{AUC.sh}
\alias{AUC.sh}
\alias{spec.sh}
\alias{sens.sh}

\title{AUC estimation proposed by Shong and Zhou}
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
  \item{lp}{The vector of linear predictors estimated from the training data.}
  \item{lpnew}{The vector of linear predictors obtained from the test data.}
  \item{times}{A vector of time points at which to compute the AUC.}
  \item{type}{A string defining the type of true positive rate (TPR):
              \code{"incident"} refers to incident TPR ,
              \code{"cumulative"} refers to cumulative TPR.}
  \item{savesensspec}{A logical specifying whether sensitivities and specificities
                      should be saved.}
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


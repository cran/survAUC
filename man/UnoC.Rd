\name{UnoC}
\alias{UnoC}

\title{C-statistic by Uno et al.}
\description{
Censoring-adjusted C-statistic by Uno et al.
}
\usage{
UnoC(Surv.rsp, Surv.rsp.new, lpnew, time = NULL)
}
\arguments{
  \item{Surv.rsp}{A \code{Surv(.,.)} object containing to the outcome of the training data.}
  \item{Surv.rsp.new}{A \code{Surv(.,.)} object containing the outcome of the test data.}
  \item{lpnew}{The vector of predictors obtained from the test data.}
  \item{time}{A positive number restricting the upper limit of the time range under consideration.}
}
\details{
  This function implements the censoring-adjusted C-statistic proposed by
  Uno et al. (2011). It has the same interpretation as Harrell's C for
  survival data (implemented in the \code{rcorr.cens} function of the \bold{Hmisc}
  package).
  
  Uno's estimator is based on inverse-probability-of-censoring
  weights and does not assume a specific working model for deriving the predictor
  \code{lpnew}. It is assumed, however, that there is a one-to-one
  relationship between the predictor and the expected survival times conditional
  on the predictor. Note that the estimator implemented in \code{UnoC}
  is restricted to situations where the random censoring assumption holds.
  
}
\value{
  The estimated C-statistic.
}
\references{

Harrell, F. E., R. M. Califf, D. B. Pryor, K. L. Lee and R. A. Rosati (1982). \cr Evaluating
the yield of medical tests.\cr \emph{Journal of the American Medical Association}
\bold{247}, 2543--2546.\cr

Harrell, F. E., K. L. Lee, R. M. Califf, D. B. Pryor and R. A. Rosati (1984). \cr Regression
modeling strategies for improved prognostic prediction.\cr \emph{Statistics in Medicine}
\bold{3}, 143--152.\cr

Uno, H., T. Cai T, M. J. Pencina, R. B. D'Agostino and W. L. Wei (2011). \cr On the
C-statistics for evaluating overall adequacy of risk prediction procedures with censored
survival data.\cr \emph{Statistics in Medicine} \bold{30}, 1105--1117.\cr

}

\seealso{
  \code{\link{GHCI}}, \code{\link{AUC.sh}}, \code{\link{IntAUC}}
}
\examples{

TR <- ovarian[1:16,]
TE <- ovarian[17:26,]
train.fit  <- coxph(Surv(futime, fustat) ~ age,
                    x=TRUE, y=TRUE, method="breslow", data=TR)

lpnew <- predict(train.fit, newdata=TE)
Surv.rsp <- Surv(TR$futime, TR$fustat)
Surv.rsp.new <- Surv(TE$futime, TE$fustat)              

Cstat <- UnoC(Surv.rsp, Surv.rsp.new, lpnew)
Cstat

}

\keyword{classif}
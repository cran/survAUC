\name{IntAUC}
\alias{IntAUC}
\title{Integration of time-dependent AUC curves}
\description{
Compute summary measures of a time-dependent AUC curve
}
\usage{
IntAUC(AUC, times, S, tmax, auc.type="cumulative")
}
\arguments{
  \item{AUC}{A vector of AUCs.}
  \item{times}{The vector of time points corresponding to \code{AUC}.}
  \item{S}{A vector of survival probabilities corresponding to \code{times}.}
  \item{tmax}{A number specifying the upper limit of the time
              range for which to compute the summary measure.}
  \item{auc.type}{A string defining the type of AUC. 'cumulative' refers
                  to cumulative/dynamic AUC, 'incident' refers to incident/dynamic
                  AUC.}
}
\details{
  This function calculates the integral under a time-dependent AUC curve (\dQuote{IAUC}
  measure) using the integration limits [0, \code{tmax}]. The values of the AUC curve are
  specified via the \code{AUC} argument.
  
  In case \code{auc.type = "cumulative"} (cumulative/dynamic IAUC), the values of
  \code{AUC} are weighted by the estimated probability density of
  the time-to-event outcome. In case \code{auc.type = "incident"} (incident/dynamic
  IAUC), the values of \code{AUC} are weighted by 2 times the product of the estimated
  probability density and the (estimated) survival function of the time-to-event outcome.
  The survival function has to be specified via the \code{S} argument.
  
  As shown by Heagerty and Zheng (2005), the incident/dynamic version of IAUC
  can be interpreted as a global concordance index measuring the probability
  that observations with a large predictor value have a shorter survival time than
  observations with a small predictor value. The incident/dynamic version of IAUC
  has the same interpretation as Harrell's C for survival data.

}
\value{
A scalar number corresponding to the summary measure of interest.
}
\references{

Harrell, F. E., R. M. Califf, D. B. Pryor, K. L. Lee and R. A. Rosati (1982). \cr Evaluating
the yield of medical tests.\cr \emph{Journal of the American Medical Association}
\bold{247}, 2543--2546.\cr

Harrell, F. E., K. L. Lee, R. M. Califf, D. B. Pryor and R. A. Rosati (1984). \cr Regression
modeling strategies for improved prognostic prediction.\cr \emph{Statistics in Medicine}
\bold{3}, 143--152.\cr

Heagerty, P. J. and Y. Zheng (2005). \cr Survival model predictive accuracy and
ROC curves.\cr \emph{Biometrics} \bold{61}, 92--105.\cr

}

\seealso{
  \code{\link{AUC.cd}}, \code{\link{AUC.sh}}, \code{\link{AUC.uno}},
  \code{\link{AUC.hc}}
}

\keyword{classif} \keyword{manip}

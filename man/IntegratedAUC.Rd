\name{IntAUC}
\alias{IntAUC}
\title{IntAUC}
\description{
Compute summary measures (discrimination indices) from a time-dependent
AUC curve
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
\value{
A scalar number corresponding to the summary measure of interest.
}
\seealso{
  \code{\link{AUC.cd}}, \code{\link{AUC.sh}}, \code{\link{AUC.uno}},
  \code{\link{AUC.hc}}
}

\keyword{classif} \keyword{manip}

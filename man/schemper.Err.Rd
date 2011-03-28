\name{schemper}
\alias{schemper}
\title{Distance-based Estimator of Survival Predictive Accuracy proposed
by Schemper and Henderson}
\description{
Schemper and Henderson's Estimator of the Absolute
Deviation between Survival Functions (Prediction Error Curve)
}
\usage{
schemper(train.fit, traindata, newdata)
}
\arguments{
  \item{train.fit}{A \code{cph} object corresponding to the fit of a Cox model.}
  \item{traindata}{A data frame corresponding to the set of training data.}
  \item{newdata}{A data frame corresponding to the set of test data.}

}
\details{
This code has been adapted from Lusa et al. (2007).
}

\references{

Lusa, L., R. Miceli and L. Mariani (2007).\cr Estimation of predictive accuracy in 
survival analysis using R and S-PLUS.\cr 
\emph{Computer Methods and Programms in Biomedicine} \bold{87}, 132--137.\cr


Schemper, M. and R. Henderson (2000).\cr Predictive accuracy and explained
variation in Cox regression.\cr \emph{Biometrics} \bold{56}, 249--255.\cr

}

\seealso{
  \code{\link{IntAUC}}, \code{\link{predErr}}, \code{\link{OXS}}
}

\keyword{classif}


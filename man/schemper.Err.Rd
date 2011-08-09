\name{schemper}
\alias{schemper}
\title{Distance-based estimator of survival predictive accuracy proposed
by Schemper and Henderson}
\description{
Schemper and Henderson's estimator of the absolute
deviation between survival functions
}
\usage{
schemper(train.fit, traindata, newdata)
}
\arguments{
  \item{train.fit}{A \code{cph} object containing the fit of a Cox proportional hazards model.}
  \item{traindata}{A data frame containing the set of training data.}
  \item{newdata}{A data frame containing the set of test data.}

}
\details{
This code has been adapted from Lusa et al. (2007). Schemper and Henderson's estimator
(as implemented by Lusa et al. 2007) assumes that predictions of the time-to-event outcome
were obtained from a Cox proportional hazards model. The estimator is valid as long as the Cox
model is specified correctly.

Technical details:
\itemize{
  \item The Cox model has to be estimated via the \code{cph} function of the
        \bold{Design} package.
  \item The survival times and the censoring indicators have to be
        labelled \dQuote{time} and \dQuote{status}, respectively (see example below).
  \item In contrast to the other estimators implemented in the \bold{survAUC} package,
        \code{schemper} does not estimate the survival function of the censoring distribution
        from the training data but from the test data.
}
For details on the estimator and its implementation, we refer to
Schemper and Henderson (2000) and Lusa et al. (2007).

}
\value{
  \code{schemper} returns a list with the following components:
  \item{Model}{The call to \code{cph}.}
  \item{D}{The estimator of predictive accuracy obtained from the
           covariate-free null model.}
  \item{Dx}{The estimator of predictive accuracy obtained from the
            Cox model.}
  \item{V}{The estimator of relative gains in predictive accuracy.}
  \item{Mhat}{The absolute distance estimator obtained from the
              Cox model (evaluated at the event times of the test data).}
  \item{Mhat.0}{The absolute distance estimator obtained from the
              covariate-free null model (evaluated at the event times of
              the test data).}
  \item{timep}{The event times of the test data.}
}

\references{

Schemper, M. and R. Henderson (2000).\cr Predictive accuracy and explained
variation in Cox regression.\cr \emph{Biometrics} \bold{56}, 249--255.\cr

Lusa, L., R. Miceli and L. Mariani (2007).\cr Estimation of predictive accuracy in 
survival analysis using R and S-PLUS.\cr 
\emph{Computer Methods and Programms in Biomedicine} \bold{87}, 132--137.\cr

}

\seealso{
  \code{\link{IntAUC}}, \code{\link{predErr}}, \code{\link{OXS}}
}

\examples{

library(Design)
ovarian$time <- ovarian$futime
ovarian$status <- ovarian$fustat
set.seed(2011)
trobs <- sample(1:26,16)
TR <- ovarian[trobs,]
TE <- ovarian[-trobs,]
train.fit  <- cph(Surv(time, status) ~ age,
                  x=TRUE, y=TRUE, method="breslow", data=TR)

schemper(train.fit, TR, TE)

}

\keyword{classif}


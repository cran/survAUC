\name{GHCI}
\alias{GHCI}

\title{Gonen and Heller's Concordance Index for Cox models}
\description{
Gonen and Heller's Concordance Index for Cox proportional hazards models
}
\usage{
GHCI(lpnew)
}
\arguments{
  \item{lpnew}{The vector of predictors obtained from the test data.}
}
\details{
  This function implements the concordance probability estimator proposed by
  Gonen and Heller (2005). It has the same interpretation as Harrell's C for
  survival data (implemented in the \code{rcorr.cens} function of the \bold{Hmisc}
  package).
  
  The results obtained from \code{GHCI} are valid as long as \code{lpnew} is the
  predictor of a correctly specified Cox proportional hazards model. In this case, the
  estimator remains valid even if the censoring times depend on the values of
  the predictor.
  
  Note that the smoothed version of \code{GHCI}, which is proposed in Section 3 of
  Gonen and Heller (2005), is not implemented in R package \bold{survAUC}.
  
}
\value{
  A length-one numeric vector containing the concordance probability estimate.
}
\references{

Harrell, F. E., R. M. Califf, D. B. Pryor, K. L. Lee and R. A. Rosati (1982). \cr Evaluating
the yield of medical tests.\cr \emph{Journal of the American Medical Association}
\bold{247}, 2543--2546.\cr

Harrell, F. E., K. L. Lee, R. M. Califf, D. B. Pryor and R. A. Rosati (1984). \cr Regression
modeling strategies for improved prognostic prediction.\cr \emph{Statistics in Medicine}
\bold{3}, 143--152.\cr

Gonen, M. and G. Heller (2005). \cr Concordance probability and discriminatory
power in proportional hazards regression.\cr \emph{Biometrika}
\bold{92}, 965--970.\cr

}

\seealso{
  \code{\link{AUC.sh}}, \code{\link{IntAUC}}
}
\examples{
data(cancer,package="survival")
TR <- ovarian[1:16,]
TE <- ovarian[17:26,]
train.fit  <- survival::coxph(survival::Surv(futime, fustat) ~ age,
                    x=TRUE, y=TRUE, method="breslow", data=TR)

lpnew <- predict(train.fit, newdata=TE)
                 
GHCI(lpnew)

}

\keyword{classif}


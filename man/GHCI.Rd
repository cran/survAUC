\name{GHCI}
\alias{GHCI}

\title{Gonen and Heller's Concordance Index for Cox models}
\description{
Gonen and Heller's Concordance Index for Cox models
}
\usage{
GHCI(lpnew)
}
\arguments{
  \item{lpnew}{The vector of linear predictors obtained from the test data.}
}
\references{

Gonen, M. and G. Heller (2005). \cr Concordance probability and discriminatory
power in proportional hazards regression.\cr \emph{Biometrika}
\bold{92}, 965--970.\cr

}

\seealso{
  \code{\link{AUC.sh}}, \code{\link{IntAUC}}
}
\examples{

TR <- ovarian[1:16,]
TE <- ovarian[17:26,]
train.fit  <- coxph(Surv(futime, fustat) ~ age,
                    x=TRUE, y=TRUE, method="breslow", data=TR)

lpnew <- predict(train.fit, newdata=TE)
                 
GHCI(lpnew)

}

\keyword{classif}


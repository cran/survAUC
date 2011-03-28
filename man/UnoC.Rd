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
  \item{lpnew}{The vector of linear predictors obtained from the test data.}
  \item{time}{A positive number restricting the upper limit of the time range under consideration.}
}
\references{

Uno, H., T. Cai T, M. J. Pencina, R. B. D'Agostino and W. L. Wei (2011). \cr On the
C-statistics for evaluating overall adequacy of risk prediction procedures with censored
survival data.\cr \emph{Statistics in Medicine}, to appear.\cr

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
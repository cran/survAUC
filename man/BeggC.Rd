\name{BeggC}
\alias{BeggC}

\title{C-statistic by Begg et al.}
\description{
C-statistic by Begg et al.
}
\usage{
BeggC(Surv.rsp, Surv.rsp.new, lp, lpnew)
}
\arguments{
  \item{Surv.rsp}{A \code{Surv(.,.)} object containing to the outcome of the training data.}
  \item{Surv.rsp.new}{A \code{Surv(.,.)} object containing the outcome of the test data.}
  \item{lp}{The vector of predictors estimated from the training data.}
  \item{lpnew}{The vector of predictors obtained from the test data.}
}
\details{
  This function implements the C-statistic proposed by
  Begg et al. (2000). It has the same interpretation as Harrell's C for
  survival data (implemented in the \code{rcorr.cens} function of the \bold{Hmisc}
  package).  \code{BeggC} is restricted to Cox regression.
  Specifically, it is assumed that \code{lp} and \code{lpnew} are the predictors of
  a Cox proportional hazards model. Estimates obtained from
  \code{BeggC} are valid as long as the Cox model is specified correctly.
}
\value{
  The estimated C-statistic.
}
\references{

Begg, B. C., L. D. Craemer, E. S. Venkatraman and J. Rosai (2000). \cr Comparing tumor
staging and grading systems: a case study and a review of the issues, using thymoma as a model.\cr
\emph{Statistics in Medicine} \bold{19}, 1997--2014.\cr

}

\seealso{
  \code{\link{UnoC}}, \code{\link{GHCI}}, \code{\link{AUC.sh}}, \code{\link{IntAUC}}
}
\examples{
data(cancer,package="survival")
TR <- ovarian[1:16,]
TE <- ovarian[17:26,]
train.fit  <- survival::coxph(survival::Surv(futime, fustat) ~ age,
                    x=TRUE, y=TRUE, method="breslow", data=TR)

lp <- predict(train.fit)
lpnew <- predict(train.fit, newdata=TE)
Surv.rsp <- survival::Surv(TR$futime, TR$fustat)
Surv.rsp.new <- survival::Surv(TE$futime, TE$fustat)              

Cstat <- BeggC(Surv.rsp, Surv.rsp.new, lp, lpnew)
Cstat

}

\keyword{classif}

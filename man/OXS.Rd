\name{OXS}
\alias{OXS}
\alias{Nagelk}
\alias{XO}
\title{R2-type coefficients for Cox proportional hazards models}
\description{
R2-type Coefficients for Cox proportional hazards models
}
\usage{
OXS(Surv.rsp, lp, lp0)
Nagelk(Surv.rsp, lp, lp0)
XO(Surv.rsp, lp, lp0)
}
\arguments{
  \item{Surv.rsp}{A \code{Surv(.,.)} object containing to the outcome of the test data.}
  \item{lp}{The vector of predictors.}
  \item{lp0}{The vector of predictors obtained from the covariate-free null model.}
}
\details{
  The \code{OXS}, \code{Nagelk} and \code{XO} functions implement three types of R2 coefficients for right-censored
  time-to-event data: (a) The coefficient proposed by O'Quigley et al. (2005) (\code{OXS}), (b)
  the coefficient proposed by Nagelkerke (1991) (\code{Nagelk}) and (c) the
  coefficient proposed by Xu and O'Quigley (1999) (\code{XO}).

  Because the \code{OXS}, \code{Nagelk} and \code{XO} functions assume that \code{lp} and \code{lpnew} were
  derived from a correctly specified Cox proportional hazards model, estimates obtained from these functions are only
  valid if the Cox model holds.

}
\references{

Nagelkerke, N. J. D. (1991).\cr A note on a general definition of the coefficient of
determination.\cr \emph{Biometrika} \bold{78}, 691--692.\cr

O'Quigley, J., R. Xu, and J. Stare (2005).\cr Explained randomness in proportional
hazards models.\cr \emph{Statistics in Medicine} \bold{24}, 479--489.\cr

Xu, R. and J. O'Quigley (1999).\cr A measure of dependence for proportional
hazards models.\cr \emph{Journal of Nonparametric Statistics} \bold{12}, 83--107.\cr

}
\value{
  The estimated R2 coefficient.
}

\seealso{
  \code{\link{predErr}}, \code{\link{schemper}}, \code{\link{GHCI}}
}
\examples{

TR <- ovarian[1:16,]
TE <- ovarian[17:26,]
train.fit  <- coxph(Surv(futime, fustat) ~ age,
                    x=TRUE, y=TRUE, method="breslow", data=TR)

model0 <- coxph(Surv(futime, fustat)~1, data=TR)
model1 <- coxph(Surv(futime, fustat)~age, data=TR)
f0 <- rep(0,nrow(TE))
f1 <- predict(model1, newdata=TE)               
Surv.res <- Surv(TE$futime, TE$fustat)

OXS(Surv.res, f1, f0)
Nagelk(Surv.res, f1, f0)
XO(Surv.res, f1, f0)

}

\keyword{classif}


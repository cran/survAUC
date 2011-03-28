\name{OXS}
\alias{OXS}
\alias{Nagelk}
\alias{XO}
\title{R2-type Coefficients for Time-to-event Data}
\description{
R2-type Coefficients for Time-to-event Data
}
\usage{
OXS(Surv.rsp, lp, lp0)
Nagelk(Surv.rsp, lp, lp0)
XO(Surv.rsp, lp, lp0)
}
\arguments{
  \item{Surv.rsp}{A \code{Surv(.,.)} object containing to the outcome of the test data.}
  \item{lp}{The vector of linear predictors.}
  \item{lp0}{The vector of linear predictors obtained from the null model.}
}
\references{

Nagelkerke, N. J. D. (1991).\cr A note on a general definition of the coefficient of
determination.\cr \emph{Biometrika} \bold{78}, 691--692.\cr

O'Quigley, J., R. Xu, and J. Stare (2005).\cr Explained randomness in proportional
hazards models.\cr \emph{Statistics in Medicine} \bold{24}, 479--489.\cr

Xu, R. and J. O'Quigley (1999).\cr A measure of dependence for proportional
hazards models.\cr \emph{Journal of Nonparametric Statistics} \bold{12}, 83--107.\cr

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


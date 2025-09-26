###########################################################
###             Schemper-Henderson estimator
###########################################################

### Programm fuer Cox-Modell (aus Lusa et al. 2007)
### liefert Objekt mit Zeitpunkten wie survfit ("timep")
### und prediction error curve ("Mhat")
### -> mit seval, etc. weiterverwertbar

### train.fit = fit von cph()



#' @title Distance-based estimator of survival predictive accuracy proposed by
#' Schemper and Henderson
#' 
#' @description Schemper and Henderson's estimator of the absolute deviation between
#' survival functions
#' 
#' @details This code has been adapted from Lusa et al. (2007). Schemper and Henderson's
#' estimator (as implemented by Lusa et al. 2007) assumes that predictions of
#' the time-to-event outcome were obtained from a Cox proportional hazards
#' model. The estimator is valid as long as the Cox model is specified
#' correctly.
#' 
#' Technical details: \itemize{ \item The Cox model has to be estimated via the
#' \code{cph} function of the \bold{Design} package.  \item The survival times
#' and the censoring indicators have to be labelled \dQuote{time} and
#' \dQuote{status}, respectively (see example below).  \item In contrast to the
#' other estimators implemented in the \bold{survAUC} package, \code{schemper}
#' does not estimate the survival function of the censoring distribution from
#' the training data but from the test data. } For details on the estimator and
#' its implementation, we refer to Schemper and Henderson (2000) and Lusa et
#' al. (2007).
#' 
#' @param train.fit A \code{cph} object containing the fit of a Cox
#' proportional hazards model.
#' @param traindata A data frame containing the set of training data.
#' @param newdata A data frame containing the set of test data.
#' @return \code{schemper} returns a list with the following components:
#' \item{Model}{The call to \code{cph}.} \item{D}{The estimator of predictive
#' accuracy obtained from the covariate-free null model.} \item{Dx}{The
#' estimator of predictive accuracy obtained from the Cox model.} \item{V}{The
#' estimator of relative gains in predictive accuracy.} \item{Mhat}{The
#' absolute distance estimator obtained from the Cox model (evaluated at the
#' event times of the test data).} \item{Mhat.0}{The absolute distance
#' estimator obtained from the covariate-free null model (evaluated at the
#' event times of the test data).} \item{timep}{The event times of the test
#' data.}
#' @seealso \code{\link{IntAUC}}, \code{\link{predErr}}, \code{\link{OXS}}
#' @references
#' 
#' Schemper, M. and R. Henderson (2000).\cr Predictive accuracy and explained
#' variation in Cox regression.\cr \emph{Biometrics} \bold{56}, 249--255.\cr
#' 
#' Lusa, L., R. Miceli and L. Mariani (2007).\cr Estimation of predictive
#' accuracy in survival analysis using R and S-PLUS.\cr \emph{Computer Methods
#' and Programms in Biomedicine} \bold{87}, 132--137.\cr
#' @keywords classif
#' @examples
#' 
#' data(cancer,package="survival")
#' ovarian$time <- ovarian$futime
#' ovarian$status <- ovarian$fustat
#' set.seed(2011)
#' trobs <- sample(1:26,16)
#' TR <- ovarian[trobs,]
#' TE <- ovarian[-trobs,]
#' train.fit  <- rms::cph(survival::Surv(time, status) ~ age,
#'                   x=TRUE, y=TRUE, method="breslow", data=TR)
#' 
#' schemper(train.fit, TR, TE)
#' 
#' 
#' @export schemper
schemper <- function(train.fit, traindata, newdata)
{
	if(!inherits(train.fit,"rms"))
		stop("\nThe Cox model has to be estimated via the cph function of the rms package.\n")
    f.Mt <- function(tempo, tutti.tempi, stima.surv, tempi.evento,
					 Stj, ind.censura, num.sogg)
		{
			Stj1 <- unique(Stj[tempi.evento == tempo])
			primo <- rep(1 - Stj1, num.sogg)
			primo[tutti.tempi <= tempo] <- 0
			secondo <- Stj1 * (1 - ind.censura)
			secondo[tutti.tempi > tempo] <- 0
			terzo <- ind.censura * (((1 - Stj1) * Stj1)/stima.surv + Stj1 * (1 - Stj1/stima.surv))
			terzo[tutti.tempi > tempo] <- 0
			terzo[is.na(terzo)] <- 0
			ris <- primo + secondo + terzo
			return(sum(ris)/num.sogg)
		}
    f.Mt.cox <- function(tempo, tutti.tempi, stima.surv, tempi.evento,
						 Stj0, ind.censura, num.sogg, lin.pred)
		{
			Stj00 <- unique(Stj0[tempi.evento == tempo])
			Stj1 <- Stj00^exp(lin.pred)
			primo <- 1 - Stj1
			primo[tutti.tempi <= tempo] <- 0
			secondo <- Stj1 * (1 - ind.censura)
			secondo[tutti.tempi > tempo] <- 0
			terzo <- ind.censura * (((1 - Stj1) * Stj1)/stima.surv + Stj1 * (1 - Stj1/stima.surv))
			terzo[tutti.tempi > tempo] <- 0
			terzo[is.na(terzo)] <- 0
			ris <- primo + secondo + terzo
			return(sum(ris)/num.sogg)
		}
    f.assegna.surv <- function(tempo, tempi.eventi)
		{
			if (any(tempo == tempi.eventi)) {
				pos <- (c(1:length(tempi.eventi)) * as.numeric(tempo == tempi.eventi))
				pos <- pos[pos != 0]
			}
			else {
				tmp <- (tempo - tempi.eventi)
				if (all(tmp < 0))
					pos <- NA
				else {
					tmp[tmp < 0] <- Inf
					pos <- order(tmp)[1]
				}
			}
			return(pos)
		}
    tsurv <- as.numeric(newdata$time)
    surv  <- as.numeric(newdata$status)
    lin.pred <- rms::predictrms(train.fit, newdata,"lp")
    num.sogg <- length(tsurv)
    km <- survival::survfit(survival::Surv(tsurv, surv) ~ 1)
    km.fit <- survival::survfit(survival::Surv(time, status) ~ 1, data=traindata)
    tempi.eventi <- km$time[km$n.event != 0]
    pos.surv <- apply(as.matrix(tsurv), 1, f.assegna.surv, tempi.eventi)
    surv.tj <- stats::approx(km.fit$time,km.fit$surv,xout =tempi.eventi , method = "constant", f = 0, yleft=1,yright=min(km.fit$surv, na.rm=T))$y
    surv.tot.km <- (surv.tj)[pos.surv]
    ind.censura <- as.numeric(!as.logical(surv))
    Mt <- apply(as.matrix(tempi.eventi), 1, f.Mt, tsurv, surv.tot.km,
				tempi.eventi, surv.tj, ind.censura, num.sogg)
    numero.eventi <- km$n.event[km$n.event != 0]
    surv0.tj.cox <- stats::approx(train.fit$time,train.fit$surv,xout =tempi.eventi , method = "constant", f = 0, yleft=1,yright=min(km.fit$surv, na.rm=T))$y
    surv0.tot.cox <- (surv0.tj.cox)[pos.surv]
    surv.tot.cox <- surv0.tot.cox^exp(lin.pred)
    Mtx <- apply(as.matrix(tempi.eventi), 1, f.Mt.cox, tsurv,
				 surv.tot.cox, tempi.eventi, surv0.tj.cox, ind.censura,
				 num.sogg, lin.pred)
    Gkm <- survival::survfit(survival::Surv(tsurv, ind.censura) ~ 1)
    tempi.censure <- Gkm$time[Gkm$n.event != 0]
    if (!length(tempi.censure))
	cens.tot.km <- rep(1, length(tempi.eventi))
    else {
        pos.surv.censure <- apply(as.matrix(tempi.eventi), 1,
								  f.assegna.surv, tempi.censure)
        cens.tot.km <- (Gkm$surv[Gkm$n.event != 0])[pos.surv.censure]
        cens.tot.km[tempi.eventi < min(Gkm$time[Gkm$n.event !=
									   0])] <- 1
    }
    pesi <- numero.eventi/cens.tot.km
    peso.tot <- sum(pesi)
    D <- sum(Mt * pesi)/peso.tot
    Dx <- sum(Mtx * pesi)/peso.tot
    V <- (D - Dx)/D
    return(list(Model = train.fit$call, D = D, Dx = Dx, V = V,  Mhat=Mtx, Mhat.0=Mt, timep=tempi.eventi))
}



##############################################################
## C-Statistic suggest by Uno
##############################################################
## Surv.rsp		- Surv(.,.) Outcome of training or test data
## lpnew			- vector of linear predictors of training or test data
## time		- time point



UnoC <- function(Surv.rsp, Surv.rsp.new, lpnew, time = NULL)
{
	if(is.null(time)){
		tau <- max(Surv.rsp.new[,1])
	}else{
		tau <- time
	}
	time <- Surv.rsp[,1]
	event <- 1-Surv.rsp[,2]

	time.new <- Surv.rsp.new[,1]
	event.new <- Surv.rsp.new[,2]
	
	n <- length(time)
	n.new <- length(time.new)
	n_lp <- length(lpnew)
	n_tau <- length(tau)
	if(n.new != n_lp)
		stop(" 'Surv.rsp' and 'linear predictors' must have the same length!\n")
	if(n_tau > 1){
		UnoC <- vector("numeric",length=n_tau)
	}else{
		UnoC <- 0
	}
	ans <- .C("C_UnoC",
			  as.numeric(time),
			  as.numeric(event),
			  as.integer(n),
			  as.numeric(time.new),
			  as.numeric(event.new),
			  as.integer(n.new),
			  as.numeric(lpnew),
			  as.numeric(tau),
			  as.integer(n_tau),
			  as.numeric(UnoC))
	#No longer needed since the symbol is registered in the NAMESPACE
	#          ,PACKAGE="survAUC")
	ans[[10]]
}






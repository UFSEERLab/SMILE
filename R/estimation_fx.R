#' LIZ negative log likelihood
#'
#' to be used with the estimation function
#'
#' @param pars list of parameters to estimate
#' @param period periodicity of infection in weeks
#' @param years number of years
#' @param SMILE.obs observation data
#'
#' @returns description
#' @export
#'
#' @examples
#'
LIZ_negll	<-	function(pars=c(theta,tau,b0,b1),period,years,SMILE.obs){

  theta	<-	exp(pars[1])
  tau		<-	exp(pars[2])
  b0		<-	pars[3]
  b1		<-	pars[4]
  # Likelihood of infection seasonality as a function of climate

  #wt.bar.hat	<-	mean(clim)
  #clim.pred	<-	infect2clim(wt.bar=wt.bar.hat,b0=b0,b1=b1,kw=kw,t=1:length(clim),period=period)
  #clim.ssq	<-	sum((clim-clim.pred)^2)

  SMILE.pred	<-	smile_main(b0,b1,period,theta,tau,years, sigmaa = 0.92^(1/52))
  ts2keep		<-	names(SMILE.pred)%in%names(SMILE.obs)
  SMILE.pred	<-	SMILE.pred[which(ts2keep==TRUE)]
  loglik.ls	<-	vector("list",length(SMILE.obs))

  for(i in 1:length(SMILE.obs)){

    loglik.ls[[i]]	<-	stats::dpois(SMILE.obs[[i]],SMILE.pred[[i]],log=TRUE)

  }

  loglik.ls	<-	lapply(loglik.ls,sum)

  loglik	<-	sum(unlist(loglik.ls))

  if(!is.finite(loglik)){loglik=-.Machine$double.xmax}
  negll	<-	-loglik #+ clim.ssq
  return(negll)

}

#' Parameter estimation from data
#'
#' @param b0 starting value
#' @param b1 starting value
#' @param theta starting value
#' @param tau starting value
#' @param period starting value
#' @param SMILE.obs data
#' @param method set to BFGS
#'
#' @returns description
#' @export
#'
#' @examples
#'
SMILE_param_estim	<-	function(b0,b1,theta,tau,period,SMILE.obs,method="BFGS"){

  years		<- 	length(SMILE.obs[[1]])/52
  pars		<-	c(theta,tau,b0,b1)

  optim.res	<-	stats::optim(par=pars,fn=LIZ_negll,method=method,SMILE.obs=SMILE.obs,
                            years=years,period=period)

  theta.hat	<-	exp(optim.res$par[1])
  tau.hat		<-	exp(optim.res$par[2])
  b0.hat		<-	optim.res$par[3]
  b1.hat		<-	optim.res$par[4]


  neg.ll		<-	optim.res$value

  results		<-	c(tau=tau.hat,theta=theta.hat,b0=b0.hat,b1=b1.hat,loglik=-neg.ll)
  return(results)
}

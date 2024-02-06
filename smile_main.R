## Code adapted from the original smile model
## building a main function that allows for the different smile scenarios
## author: javirudolph
##  date: Feb 2, 2024


# I'm putting all the parameters as options so we can fully modify it

smile_fx <- function(b0,b1,period,theta,tau,years,
                     # fixed parameters
                     # Fixed parameters
                     alpha =	1/52,
                     zeta_novax	=	0.88,
                     gamma = 0.9868,
                     sigmaa = 1,
                     psi = 1,
                     N1 = 5000,
                     K = 5000,
                     # adding a parameter to add vaccine data if we want
                     vax = NULL,
                     beta_0 = NULL,
                     beta_1 = NULL,
                     # To remove seasonality, give b a value
                     b_fixed = NULL,
                     # To remove host age structuring, change to false
                     age_struc = TRUE,
                     stochastic = FALSE,
                     LIZ_init = 1, 
                     output_df = FALSE) {
  
  if(is.null(vax) == FALSE) {
    if(is.null(beta_0) == TRUE){
      print("stop! you need beta_0 and beta_1 to calculate vaccine coverage")
    } else {
      # zeta_0 <- -0.25 * (vax - 1)^2 + 0.999
      zeta_0 <- 1 / (1 + exp(-(beta_0 + beta_1 * vax)))
      zeta <- rep(zeta_0,each=52)
      zeta <- c(zeta[1],zeta)
    }
  } else {
    zeta <- rep(zeta_novax, years*52 + 1)
  }
  

  
  n.weeks		<-	years*52 + 1
  S<-	M<-	I<-	L<-E <-N<-lambda<- array(0,dim=c(n.weeks),dimnames=list(1:(n.weeks)))
  N[1]		<-	N1
  S[1]		<-	N[1]	
  L[1]	<-	LIZ_init; E[1]	<-	L[1]*psi	
  
  for(t in 2:n.weeks){
    
    tm1	<-	t-1
    # b	<-	b.season(b0,b1,period,t)
    
    if(is.null(b_fixed) == TRUE) {
      b <- b.season(b0,b1,period,t)
    } else {
      b <- b_fixed
    }
    
    lambda[tm1]	<-	lambda.t(theta=theta,tau=tau,b=b,E=E[tm1])
    
    if(age_struc == FALSE) {
      births.happen = 0
      rep.prob = 0
    } else {
      births.happen	<-	as.numeric(t%%52==0)
      rep.prob	<-	rho_n(N[tm1], K)
    }
    
    if(stochastic == TRUE) {
      
      I[t]	<-	rbinom(1,S[tm1],lambda[tm1])
      births	<-	rbinom(1,N[tm1],rep.prob)*births.happen
      M.surv	<-	rbinom(1,M[tm1],sigmaa)
      M.recov	<-	rbinom(1,M.surv,alpha)
      M.new	<-	rbinom(1,I[tm1],zeta)
      S[t]	<-	rbinom(1,S[tm1]-I[t],sigmaa) + M.recov + births
      M[t]	<-	M.new + (M.surv - M.recov)
      L[t]	<-	(I[tm1]-M.new)
      E[t]	<-	rpois(1,psi*L[tm1]) + rbinom(1,E[tm1],gamma)
      N[t]	<-	S[t]+M[t]
      
    } else {
      
      S[t]	<-	(S[tm1]*(1-lambda[tm1]))*sigmaa + M[tm1]*sigmaa*1/52 + rep.prob*(N[tm1])*births.happen
      I[t]	<-	S[tm1]*lambda[tm1]
      M[t]	<-	I[tm1]*zeta[tm1] + M[tm1]*sigmaa*(1-(1/52))
      L[t]	<-	I[tm1]*(1-zeta[tm1])
      E[t]	<-	psi*L[tm1] + E[tm1]*gamma
      N[t]	<-	S[t]+M[t]
    }

  }
  
  
  results		<-	list(Suceptibles=S[-1]
                   ,Immune=M[-1]
                   ,Infected=I[-1]
                   ,LIZ=L[-1]
                   ,Environment=E[-1])
  
  if(output_df == TRUE) {
    results <- data.frame(week = 1:n.weeks, S, M, I, L, E)
  }
  
  if(output_df == TRUE) {
    if(is.null(vax) == FALSE) {
      results <- data.frame(week = 1:n.weeks, S, M, I, L, E, Z = zeta)
    } else {
      results <- data.frame(week = 1:n.weeks, S, M, I, L, E)
    }
  }
  
  return(results)
}


rho_n	<-	function(N, K){
  
  0.41/(1+(N/K)^(10))
  
}

# No change to these ones -------------------------------------

# Infection probability based on Ponciano and Capistran 2011.
lambda.t	<-	function(theta,tau,b,E){
  
  1-(theta/(theta+b*E))^tau
  
}

# Introducing seasonality in the infection probability through b.

b.season	<-	function(b0,b1,period,t){
  
  exp(b0*(1+b1*cos((2*pi*t)/period)))
  
}

# Adapting the estimation functions -----------------------------------

LIZ_negll	<-	function(pars=c(theta,tau,b0,b1),period,years,SMILE.obs){
  
  theta	<-	exp(pars[1])
  tau		<-	exp(pars[2])
  b0		<-	pars[3]
  b1		<-	pars[4]
  # Likelihood of infection seasonality as a function of climate
  
  #wt.bar.hat	<-	mean(clim)
  #clim.pred	<-	infect2clim(wt.bar=wt.bar.hat,b0=b0,b1=b1,kw=kw,t=1:length(clim),period=period)
  #clim.ssq	<-	sum((clim-clim.pred)^2)
  
  SMILE.pred	<-	smile_fx(b0,b1,period,theta,tau,years, sigmaa = 0.92^(1/52))
  ts2keep		<-	names(SMILE.pred)%in%names(SMILE.obs)
  SMILE.pred	<-	SMILE.pred[which(ts2keep==TRUE)]
  loglik.ls	<-	vector("list",length(SMILE.obs))
  
  for(i in 1:length(SMILE.obs)){
    
    loglik.ls[[i]]	<-	dpois(SMILE.obs[[i]],SMILE.pred[[i]],log=TRUE)
    
  }
  
  loglik.ls	<-	lapply(loglik.ls,sum)
  
  loglik	<-	sum(unlist(loglik.ls))
  
  if(!is.finite(loglik)){loglik=-.Machine$double.xmax}
  negll	<-	-loglik #+ clim.ssq
  return(negll)
  
}

SMILE_param_estim	<-	function(b0,b1,theta,tau,period,SMILE.obs,method="BFGS"){
  
  years		<- 	length(SMILE.obs[[1]])/52
  pars		<-	c(theta,tau,b0,b1)
  
  optim.res	<-	optim(par=pars,fn=LIZ_negll,method=method,SMILE.obs=SMILE.obs
                     ,years=years,period=period)
  
  theta.hat	<-	exp(optim.res$par[1])
  tau.hat		<-	exp(optim.res$par[2])
  b0.hat		<-	optim.res$par[3]
  b1.hat		<-	optim.res$par[4]
  
  
  neg.ll		<-	optim.res$value
  
  results		<-	c(tau=tau.hat,theta=theta.hat,b0=b0.hat,b1=b1.hat,loglik=-neg.ll)
  return(results)
}




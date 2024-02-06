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
                     sigmaa = 0.92^(1/52),
                     psi = 1,
                     N1 = 5000,
                     K = 5000,
                     # adding a parameter to add vaccine data if we want
                     vax = NULL,
                     beta_0 = NULL,
                     beta_1 = NULL,
                     # To remove seasonality, give b a value
                     b_fixed = NULL,
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
    births.happen	<-	as.numeric(t%%52==0)
    rep.prob	<-	rho_n(N[tm1], K)
    
    S[t]	<-	(S[tm1]*(1-lambda[tm1]))*sigmaa + M[tm1]*sigmaa*1/52 + rep.prob*(N[tm1])*births.happen
    I[t]	<-	S[tm1]*lambda[tm1]
    M[t]	<-	I[tm1]*zeta[tm1] + M[tm1]*sigmaa*(1-(1/52))
    L[t]	<-	I[tm1]*(1-zeta[tm1])
    E[t]	<-	psi*L[tm1] + E[tm1]*gamma
    N[t]	<-	S[t]+M[t]
    
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


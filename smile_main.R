## Code adapted from the original smile model
## building a main function that allows for the different smile scenarios
## author: javirudolph
##  date: Feb 2, 2024


# I'm putting all the parameters as options so we can fully modify it

smile_fx <- function(b0,b1,period,theta,tau,years,
                     # fixed parameters
                     # Fixed parameters
                     alpha =	1/52,
                     zeta	=	0.88,
                     gamma = 0.9868,
                     sigmaa = 0.92^(1/52),
                     psi = 1,
                     N1 = 5000) {
  
  n.weeks		<-	years*52 + 1
  S<-	M<-	I<-	L<-E <-N<-lambda<- array(0,dim=c(n.weeks),dimnames=list(1:(n.weeks)))
  N[1]		<-	N1
  S[1]		<-	N[1]	
  L[1]	<-	1; E[1]	<-	L[1]*psi	
  
  for(t in 2:n.weeks){
    
    tm1	<-	t-1
    b	<-	b.season(b0,b1,period,t)
    lambda[tm1]	<-	lambda.t(theta=theta,tau=tau,b=b,E=E[tm1])
    births.happen	<-	as.numeric(t%%52==0)
    rep.prob	<-	rho.n(N[tm1])
    
    S[t]	<-	(S[tm1]*(1-lambda[tm1]))*sigmaa + M[tm1]*sigmaa*1/52 + rep.prob*(N[tm1])*births.happen
    I[t]	<-	S[tm1]*lambda[tm1]
    M[t]	<-	I[tm1]*zeta + M[tm1]*sigmaa*(1-(1/52))
    L[t]	<-	I[tm1]*(1-zeta)
    E[t]	<-	psi*L[tm1] + E[tm1]*gamma
    N[t]	<-	S[t]+M[t]
    
  }
  
  
  results		<-	list(Suceptibles=S[-1]
                   ,Immune=M[-1]
                   ,Infected=I[-1]
                   ,LIZ=L[-1]
                   ,Environment=E[-1])
  
  return(results)
}
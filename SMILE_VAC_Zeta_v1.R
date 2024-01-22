# This code is the source of functions that define the SMILE model.
# It is a compilation of functions required to evalluate the model with different 
# assumptions starting from the determinisitc simple model withou population dynamics
# or seasonal forcing and ending with the stochastic version of the model with 
# population dynamics and seasonal forcing. It also includes some useful functions 
# evaluated inside the functions to simulate the time series.
# Finally I provide the functions required to estimate the parameters of the model
# using maximum likelihood optimization through the optim function.

# Author: Juan Pablo Gomez
# Version 1.1.
# Last revised: March 7 2018.

#################################
# Infection probability based on Ponciano and Capistran 2011.
lambda.t	<-	function(theta,tau,b,E){
	
	1-(theta/(theta+b*E))^tau
	
}

# Introducing seasonality in the infection probability through b.

b.season	<-	function(b0,b1,period,t){
	
	exp(b0*(1+b1*cos((2*pi*t)/period)))
	
}

# Density dependent reproduction

rho.n	<-	function(N){
	
	0.41/(1+(N/5000)^(10))
	
}

# First case of R0

r0	<-	function(b,E,theta,tau){
	
	(b*E*tau)/theta
	
}

# Function to simulate climatic variable
# a = upper asymptote
# c = lower asymptote
# b = time required to reach half of the climatic maxima
# d = time at which the mid point between max and minimum climate is reached.
clim.func	<-	function(a,c,d,d2,b,b2,t,sd){
	
	if(missing(b2)){b2 <- b}
	if(missing(c)){c<-0}
	
	det.clim	<-	 ((a-c)/(1+exp((d-t)/b)) - (a-c)/(1+exp((d2-t)/b2))) + c
	rand.clim	<-	rnorm(length(t),det.clim,sd=sd)
	return(rand.clim)
	
}

# Function used to estimate seasonality in infection probability using climatic variables
# wt.bar = mean of climatic covariable
# b0 = as in b.season
# b1 = as in b.season
# kw = scaling coefficient
# t = time
# period = as in b.season

infect2clim	<-	function(wt.bar,b0,b1,kw,t,period){
	
	wt.bar + ((b0*b1)/kw)*cos((2*pi*t)/period)
	
}

#### modeling zeta (recovery prob.) as a function of the fraction of vaccinated
#### individuals (annually)
vacc.prop <- seq(0,1,by=0.01)
beta0 <- -3.5
beta1 <- 7.4
zeta.v <- 1/(1+exp(-(beta0+beta1*vacc.prop)))
plot(vacc.prop,zeta.v, type="l", col="red", lwd=2)


# Function of SMILE without any age structuring
# This is the most basic SMILE function in which I have assumed no population dynamics
# and no deaths other than disease related. 
# Fixed parameters are: 
# alpha= 1/52 Probability of an immune individual to become Suceptible 
# zeta = 0.88 Probability that an Infected individual becomes immune
# gamma= 0.9868 Spore decay probability
# psi = 1 Number of spores in one LIZ

# Variable parameters
# b: Number of infections caused by one LIZ if there was no dispersion effort
# tau and theta: Shape and Rate parameters of the gamma distribution defining dispersion effort.

smile.vacc1	<-	function(b,theta,tau,years,beta0,beta1,vacc.vec,N1){
  
  # Fixed parameters
  alpha	<-	1/52
  zeta0 <- 1/(1+exp(-(beta0+beta1*vacc.vec)))
  gamma	<-	0.9868
  psi		<-	1	 
  
  n.weeks		<-	years*52 + 1
  zeta <- rep(zeta0,each=52)
  zeta <- c(zeta[1],zeta)
  
  S<-	M<-	I<-	L<-E <-N<-lambda<- array(0,dim=c(n.weeks),dimnames=list(1:(n.weeks)))
  N[1]		<-	N1
  S[1]		<-	N[1]
  
  L[1]	<-	1; E[1]	<-	L[1]*psi	
  
  for(t in 2:n.weeks){
    
    tm1	<-	t-1
    lambda[tm1]	<-	lambda.t(theta=theta,tau=tau,b=b,E=E[tm1])
    
    S[t]	<-	(S[tm1]*(1-lambda[tm1])) + M[tm1]*1/52
    I[t]	<-	S[tm1]*lambda[tm1]
    M[t]	<-	I[tm1]*zeta[tm1] + M[tm1]*(1-(1/52))
    L[t]	<-	I[tm1]*(1-zeta[tm1])
    E[t]	<-	psi*L[tm1] + E[tm1]*gamma
    N[t]	<-	S[t]+M[t]
    
  }
  
  
  results		<-	cbind(Suceptibles=S[-1]
                    ,Immune=M[-1]
                    ,Infected=I[-1]
                    ,LIZ=L[-1]
                    ,Environment=E[-1]
                    ,lambda=lambda[-1])
  
  return(results)
}







#' SMILE function
#'
#' This function encompasses the original functions and assumptions. It outputs a time series for each of the compartments given a set of parameters.
#'
#' @param b0 strength of seasonality parameter
#' @param b1 strength of seasonality parameter
#' @param period periodicity of the infection in weeks
#' @param theta dispersion effort parameter
#' @param tau dispersion effort parameter
#' @param years number of years to run the simulation
#' @param alpha probability of transitioning from immune to susceptible in one week timestep
#' @param zeta_novax probability of becoming immune after infection (fixed when no vaccination)
#' @param gamma spore persistence rate
#' @param sigmaa rate of death by other causes than disease
#' @param psi number of spores per carcass
#' @param N1 starting population size
#' @param K carrying capacity for population dynamics
#' @param vax vaccination rates for each year
#' @param beta_0 survival parameter
#' @param beta_1 vaccine survival parameter
#' @param b_fixed Number of infections caused by one LIZ when dispersion effort is zero and assume no seasonality
#' @param age_struc default is TRUE, if set to false, no population dynamics are considered
#' @param stochastic default is FALSE, so the simulations are deterministic
#' @param LIZ_init Initial number of carcasses, set to default 1
#' @param rho_pop average reproduction rate to be used when simulating with population dynamics
#' @param output_df if set to TRUE, will return a dataframe instead of a list
#'
#' @returns description
#' @export
#'
#' @examples
#' SMILE::smile_fx(-30, 0.85, 3*52, 10, 10, 10)

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
                     rho_pop = NULL,
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
      b <- b_season(b0,b1,period,t)
    } else {
      b <- b_fixed
    }

    lambda[tm1]	<-	lambda_t(theta=theta,tau=tau,b=b,E=E[tm1])

    if(is.null(rho_pop) == TRUE) {
      rho_pop = 0.41
    }

    if(age_struc == FALSE) {
      births.happen = 0
      rep.prob = 0
    } else {
      births.happen	<-	as.numeric(t%%52==0)
      rep.prob	<-	rho_n(N[tm1], K, rho_pop)
    }

    if(stochastic == TRUE) {

      I[t]	<-	stats::rbinom(1,S[tm1],lambda[tm1])
      births	<-	stats::rbinom(1,N[tm1],rep.prob)*births.happen
      M.surv	<-	stats::rbinom(1,M[tm1],sigmaa)
      M.recov	<-	stats::rbinom(1,M.surv,alpha)
      M.new	<-	stats::rbinom(1,I[tm1],zeta)
      S[t]	<-	stats::rbinom(1,S[tm1]-I[t],sigmaa) + M.recov + births
      M[t]	<-	M.new + (M.surv - M.recov)
      L[t]	<-	(I[tm1]-M.new)
      E[t]	<-	stats::rpois(1,psi*L[tm1]) + stats::rbinom(1,E[tm1],gamma)
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

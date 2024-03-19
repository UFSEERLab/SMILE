

#' Calculate average reproduction rate at time t
#'
#' @param N population size at time t
#' @param K carrying capacity
#' @param rho average reproduction rate
#'
#' @returns description
#' @export
#'
#' @examples
#' rho_n(2000, 5000)
rho_n	<-	function(N, K, rho=0.41){

  rho/(1+(N/K)^(10))

}


#' Calculate infection probability for time t
#'
#' Infection probability based on Ponciano and Capistran 2011.
#'
#' @param theta disease dispersion parameter
#' @param tau disease dispersion parameter
#' @param b number of infections caused by one LIZ
#' @param E Number of spores in the environment
#'
#' @returns description
#' @export
#'
#' @examples
#' lambda_t(10, 10, 0.001, 1)
lambda_t	<-	function(theta,tau,b,E){

  1-(theta/(theta+b*E))^tau

}


#' Estimating the number of infections caused by 1 LIZ under seasonality and no dispersion effort
#'
#' Introducing seasonality in the infection probability through b.
#'
#' @param b0 strength of seasonality parameter
#' @param b1 strength of seasonality parameter
#' @param period periodicity of infection in weeks
#' @param t vector of time
#'
#' @returns description
#' @export
#'
#' @examples
#' b_season(-30, 0.85, 3*52, 1:10)
b_season	<-	function(b0,b1,period,t){

  exp(b0*(1+b1*cos((2*pi*t)/period)))

}


# Adapting the estimation functions -----------------------------------

#' Calculating local R0
#'
#' @param tau dispersion parameter
#' @param theta dispersion parameter
#' @param b number of infections caused by one LIZ
#' @param E Number of spores in the environment
#'
#' @returns description
#' @export
#'
#' @examples
#' calc_local_R0(10, 10, 0.001, 1)
calc_local_R0 <- function(tau, theta, b, E) {
  (tau*b*E)/theta
}



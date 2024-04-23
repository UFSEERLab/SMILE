#' Calculate survival probability based on vaccination rate
#'
#' This function allows us to visualize and also calculate the survival probability
#' based on a set of parameters
#'
#' @param beta_0 vector or single value. The starting point for natural immunity. A higher value means greater natural immunity
#' @param beta_1 vector or single value. The effect of the vaccine. Higher value means the vaccine has a greater positive effect on survival
#' @param vax_rate value between 0 and 1, to represent the percentage of the population that is vaccinated
#'
#' @return description
#' @export
#'
#' @examples
#' calc_survival_prob(beta_0 = -1, beta_1 = 10, vax_rate = seq(0, 1, 0.01))
#'
calc_survival_prob <- function(beta_0, beta_1, vax_rate) {
  (1 / (1 + exp(-(beta_0 + beta_1 * vax_rate))))
}



library(ggplot2)


exponential_pdf <- function(x = a, rate = lambda) {
  
  lambda = rate
  a = x
  lambda * exp(-lambda * a)
}


lomax_pdf <- function(x = a, shape = tau, scale = theta) {
  # shape parameter is the tau
  # scale parameter is theta
  theta = scale
  tau = shape
  a = x
  (tau * theta^tau) / (a + theta)^(tau + 1)
}

mean_lomax <- function(theta, tau) {
  theta / (tau - 1)
}

gamma_pdf <- function(x = a, shape = tau, rate = theta) {
  # shape parameter here is tau
  # rate parameter is theta
  # scale is 1/theta
  theta = rate
  tau = shape
  a = x
  ((theta^tau) / gamma(tau)) * a^(tau - 1) * exp(-theta*a)
  

}

mean_gamma <- function(theta, tau) {
  tau / theta
}

mean_gamma(theta = 100, tau = 10) # for theta = 100, tau = 10 -- it's 0.1, low mean dispersion effort
mean_gamma(10,10)
mean_gamma(1, 10) # High dispersion, mean dispersion effort is 10

# But the lomax is backwards... increasing theta increases dispersion effort
mean_lomax(theta = 100, tau = 10)
mean_lomax(10, 10)
mean_lomax(1, 10)

a <- seq(from = 0, to = 10, length.out = 200)

ggplot() +
  stat_function(fun = dgamma, args = list(shape = 10, rate = 10)) +
  stat_function(fun = gamma_pdf, args = list(shape = 10, rate = 10), color = "blue") +
  lims(x = c(0,5))

# As we increase the value of the scale in the lomax distribution, we increase the variance, and get a heavier tail

library(cowplot)

plot_grid(
  ggplot() +
    stat_function(fun = lomax_pdf, args = list(shape = 10, scale = 5), aes(color = "5")) +
    stat_function(fun = lomax_pdf, args = list(shape = 10, scale = 10), aes(color = "10")) +
    stat_function(fun = lomax_pdf, args = list(shape = 10, scale = 15), aes(color = "15")) +
    stat_function(fun = lomax_pdf, args = list(shape = 10, scale = 50), aes(color = "50")) +
    stat_function(fun = lomax_pdf, args = list(shape = 10, scale = 100), aes(color = "100")) +
    lims(x = c(0,50), y = c(0, 2)) +
    theme_classic() +
    labs(title = "Lomax distribution with fixed shape = 10",
         color = "Scale"),
  ggplot() +
    # stat_function(fun = gamma_pdf, args = list(shape = 10, rate = 5), aes(color = "5")) +
    # stat_function(fun = gamma_pdf, args = list(shape = 10, rate = 10), aes(color = "10")) +
    # stat_function(fun = gamma_pdf, args = list(shape = 10, rate = 15), aes(color = "15")) +
    # stat_function(fun = gamma_pdf, args = list(shape = 10, rate = 50), aes(color = "50")) +
    stat_function(fun = gamma_pdf, args = list(shape = 10, rate = 100), aes(color = "100")) +
    lims(x = c(0,50), y = c(0,2)) +
    theme_classic() +
    labs(title = "Gamma distribution with fixed shape = 10",
         color = "Rate")
)

ggplot() +
  stat_function(fun = gamma_pdf, args = list(shape = 10, rate = 100), aes(color = "100"))



ggplot() +
  stat_function(fun = lomax_pdf, args = list(shape = 2, scale = 5), aes(color = "5")) +
  stat_function(fun = lomax_pdf, args = list(shape = 2, scale = 10), aes(color = "10")) +
  stat_function(fun = lomax_pdf, args = list(shape = 2, scale = 15), aes(color = "15")) +
  stat_function(fun = lomax_pdf, args = list(shape = 2, scale = 50), aes(color = "50")) +
  stat_function(fun = lomax_pdf, args = list(shape = 2, scale = 100), aes(color = "100")) +
  lims(x = c(0,100), y = c(0, 0.1)) +
  theme_classic() +
  labs(title = "Lomax distribution with fixed shape = 2",
       color = "Scale")


ggplot() +
  stat_function(fun = dexp, args = list(rate = 1), aes(color = "1")) +
  stat_function(fun = dexp, args = list(rate = 10), aes(color = "10")) +
  stat_function(fun = dexp, args = list(rate = 100), aes(color = "100")) +
  theme_classic()



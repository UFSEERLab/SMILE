

library(ggplot2)



lomax_pdf <- function(a, theta, tau) {
  # shape parameter is the tau
  # scale parameter is theta
  (tau * theta^tau) / (a + theta)^(tau + 1)
}

mean_lomax <- function(theta, tau) {
  theta / (tau - 1)
}

gamma_pdf <- function(a, theta, tau) {
  # shape parameter here is tau
  # scale parameter is theta
  (theta^tau / gamma(tau)) * a^(tau - 1) * exp(-theta*a)
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

a <- seq(from = 0, to = 10, length.out = 100)
tau <- 10
theta <- 1
my_vals <- data.frame(a = a, lomax = lomax_pdf(a, theta, tau), gamma = gamma_pdf(a, theta, tau))


data.frame(a = seq(from = 0, to = 10, length.out = 100)) |>
ggplot() +
  # geom_line(aes(x = a, y = gamma_pdf(a, 1, 10)), color = "black") +
  # geom_line(aes(x = a, y = gamma_pdf(a, 10, 10)), color = "black", linetype = "dashed") +
  # geom_line(aes(x = a, y = gamma_pdf(a, 100, 10)), color = "black", linetype = "dotted") +
  geom_line(aes(x = a, y = lomax_pdf(a, 10, 2)), color = "red") +
  geom_line(aes(x = a, y = lomax_pdf(a, 10, 3)), color = "red", linetype = "dashed") +
  geom_line(aes(x = a, y = lomax_pdf(a, 10, 4)), color = "red", linetype = "dotted") +
  labs(y = "f(a)", x = "a") +
  theme_classic()




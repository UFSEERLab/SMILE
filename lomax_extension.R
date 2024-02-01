


lomax_pdf <- function(a, theta, tau) {
  (tau * theta^tau) / (a + theta)^(tau + 1)
}

gamma_pdf <- function(a, theta, tau) {
  (theta^tau / gamma(tau)) * a^(tau - 1) * exp(-theta*a)
}




a <- seq(from = 0, to = 20, length.out = 100)
tau <- 10
theta <- 10


library(ggplot2)
ggplot(data = data.frame(a = a, lomax = lomax_pdf(a, theta, tau), gamma = gamma_pdf(a, theta, tau))) +
  geom_line(aes(x = a, y = lomax), color = "red") +
  geom_line(aes(x = a, y = gamma), color = "black") +
  labs(title = paste("theta =", theta, "tau =", tau),
       subtitle = "Lomax in red") +
  theme_classic()






my_beta_1 = 10
my_betas_0 = seq(-11, 5, length.out = 30)

vax_prcnt <- seq(0, 1, 0.01)

data.frame(vax_prcnt, betas_0 = rep(my_betas_0, each = length(vax_prcnt)), betas_1 = my_beta_1)  %>% 
  mutate(surv_prob = (1 / (1 + exp(-(betas_0 + betas_1 * vax_prcnt))))) %>% 
  ggplot(aes(x = vax_prcnt, y = surv_prob, color = betas_0, group = betas_0)) +
  geom_line(linewidth = 1) +
  scale_color_viridis_c(option = "A") +
  scale_fill_viridis_c(option = "A") +
  theme_bw() +
  scale_x_continuous(labels = scales::label_percent(), breaks = seq(0, 1, 0.1)) +
  labs(x = "Vaccination rate", y = "Survival Probability", color = expression(beta[0]),
       title = expression("Survival probability curves with fixed"~beta[1]~"=10")) -> B



my_beta_1 = 5
my_betas_0 = seq(-11, 5, length.out = 30)

vax_prcnt <- seq(0, 1, 0.01)

data.frame(vax_prcnt, betas_0 = rep(my_betas_0, each = length(vax_prcnt)), betas_1 = my_beta_1)  %>% 
  mutate(surv_prob = (1 / (1 + exp(-(betas_0 + betas_1 * vax_prcnt))))) %>% 
  ggplot(aes(x = vax_prcnt, y = surv_prob, color = betas_0, group = betas_0)) +
  geom_line(linewidth = 1) +
  scale_color_viridis_c(option = "A") +
  scale_fill_viridis_c(option = "A") +
  theme_bw() +
  scale_x_continuous(labels = scales::label_percent(), breaks = seq(0, 1, 0.1)) +
  labs(x = "Vaccination rate", y = "Survival Probability", color = expression(beta[0]),
       title = expression("Survival probability curves with fixed"~beta[1]~"=5")) -> A


my_beta_1 = 20
my_betas_0 = seq(-11, 5, length.out = 30)

vax_prcnt <- seq(0, 1, 0.01)

data.frame(vax_prcnt, betas_0 = rep(my_betas_0, each = length(vax_prcnt)), betas_1 = my_beta_1)  %>% 
  mutate(surv_prob = (1 / (1 + exp(-(betas_0 + betas_1 * vax_prcnt))))) %>% 
  ggplot(aes(x = vax_prcnt, y = surv_prob, color = betas_0, group = betas_0)) +
  geom_line(linewidth = 1) +
  scale_color_viridis_c(option = "A") +
  scale_fill_viridis_c(option = "A") +
  theme_bw() +
  scale_x_continuous(labels = scales::label_percent(), breaks = seq(0, 1, 0.1)) +
  labs(x = "Vaccination rate", y = "Survival Probability", color = expression(beta[0]),
       title = expression("Survival probability curves with fixed"~beta[1]~"=20")) -> C

my__legend <- get_legend(C)


cowplot::plot_grid(A + theme(legend.position = "none"), 
                   B + theme(legend.position = "none"),
                   C + theme(legend.position = "none"),
                   my__legend, rel_widths = c(1,1,1,0.2), nrow=1)






# From the immunity paper:
vax_prcnt = 0

(1 / (1 + exp(-(betas_0 + betas_1 * vax_prcnt))))

calc_beta_0 <- function(zeta) {
  -log((1/zeta)-1)
}


calc_beta_0(.52)

# Disease dispersion effort -

gamma_pdf <- function(x = a, shape = tau, rate = theta) {
  # shape parameter here is tau
  # rate parameter is theta
  # scale is 1/theta
  theta = rate
  tau = shape
  a = x
  ((theta^tau) / gamma(tau)) * a^(tau - 1) * exp(-theta*a)
  
  
}

ggplot() +
  stat_function(fun = gamma_pdf, args = list(shape = 10, rate = 100), aes(color = "low dispersion"), linewidth = 1) +
  stat_function(fun = gamma_pdf, args = list(shape = 10, rate = 10), aes(color = "medium dispersion"), linewidth = 1) +
  stat_function(fun = gamma_pdf, args = list(shape = 10, rate = 1), aes(color = "high dispersion"), linewidth = 1) +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  lims(x = c(0,15))


# Seasonality


b.season	<-	function(b0,b1,period,t){
  
  exp(b0*(1+b1*cos((2*pi*t)/period)))
  
}

nweeks <- 10*52
my_b_fixed <- 0.001
my_b_season <- b.season(rep(-30, nweeks), rep(0.85, nweeks), 3*52, 1:nweeks)

data.frame(week = 1:nweeks, b = my_b_season) %>% 
  ggplot(aes(x = week, y = b)) +
  geom_path() +
  geom_hline(yintercept = my_b_fixed, color = "red") +
  theme_bw() +
  labs(y = "b - probability of transmission")
  






# Testing the stochastic simulation

# Figures

library(tidyverse)


pdf("Stochsims.pdf")
for(i in 1:100){
  do.call(cbind.data.frame, stochastic.sims[[i]]) %>% 
    rownames_to_column(
      var = "time"
    ) %>% 
    pivot_longer(-time, names_to = "compartment") -> my_df
  
  print(  my_df %>% 
            ggplot(aes(x = time, y = value, group = compartment)) +
            facet_wrap(~compartment, scale = "free") +
            geom_line() +
            labs(title = paste("Sim", i)))
}
dev.off()


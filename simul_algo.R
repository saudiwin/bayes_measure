# Simulation and evaluation of algorithm pricing strategy
# Robert Kubinec
# December 6, 2016

require(Rsolnp)
require(ggplot2)
require(dplyr)
require(tidyr)

set.seed(12062016)

n <- 1000
g <- n/10
algos <- seq(0,5,length.out=1000)


# Objective function
exp_revenue <- function(budget=NULL,out_data=NULL,sum_data=NULL) {
  out_data <- mutate(out_data,e_rev = outcome * (2 + budget*soft_latent*5 - 2.5*(budget*soft_latent)^2))
  e_profit <- sum(out_data$e_rev) - budget*n - sum(out_data$e_rev)*.25 - 1500
  return(e_profit)
}


for(i in 1:1000) {
# Create an algorithm with considerable predictive power
latent <- rnorm(n*g,algos[i],sd=0.5) 
outcome <- runif(n,0,1) < plogis(latent + rnorm(n,0,2))

# Sort so we can segment the market
out_data <- data_frame(latent = plogis(latent),as.integer(outcome)) %>% arrange(desc(latent))
out_data$groups <- factor(rep(1:10, each=n/10))
out_sum <- group_by(out_data,groups) %>% summarize(mean_latent=mean(latent)) %>% ungroup %>% 
  mutate(soft_latent=mean_latent)
out_data <- left_join(out_data,out_sum,by='groups')

max_value <- optimize(exp_revenue,interval=c(0,10000),out_data=out_data,sum_data=sum_data,
                      maximum=TRUE)
output_budget[i] <- max_value$maximum
output_profit[i] <- max_value$objective
}

data_frame(algos,Profit=output_profit,Advertising=output_budget*1000) %>%  gather(outcome_var,amount,-algos) %>%  ggplot(aes(y=amount,x=algos)) +
  geom_point() + theme_minimal() + facet_wrap(~outcome_var,scales='free_y')
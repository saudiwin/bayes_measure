# Run a bayesian measurement model on probability estimates of a binomial outcome (customer purchase)
# V0.1 Robert Kubinec
# 10/18/2016

require(dplyr)
require(tidyr)
require(magrittr)
require(rstan)
require(shinystan)
require(ggplot2)
require(plotly)
require(archivist)
require(loo)

# Data directory for archivist R model collection

repoDir <- 'data/'

# Whether to use samples from the data or the full data

sample_it <- FALSE

# Whether to run a hierarchical model

hier <- TRUE

# Name of stan code to run

stan_file <- 'measure_v1.stan'

algo_data <- data.table::fread('data/all_combine_200k.csv') %>% as_data_frame

if(sample_it==TRUE) {
  algo_data <- sample_n(algo_data,10000)
}

outcome <- algo_data$realpurch

# Standardize algorithmic output to the same scale

algos <- select(algo_data,algo1,algo2,algo3,algo4,algo5) %>% mutate_all(scale) 

# Look at raw distributions of algorithms on a standardized scale

algos %>% gather(algo_type,prob_score) %>% ggplot(aes(x=prob_score)) + geom_histogram(bins=60) +
  facet_wrap(~algo_type,scales='free')

stan_file <- stan_model(stan_file,model_name='Binomial Bayes Measurement')

# Do 800 iterations as it is hard to get this much data to compute
# Non-centered parameterization: 9500 seconds with 2 divergent transitions in 800 iterations
# Centered parameterization: 5000 seconds with no divergent transitions, although low ESS on tau
# Centered wins with more data

time1 <- Sys.time()

bayes_model <- sampling(stan_file,data=list(y=outcome,
                                            x=as.matrix(algos),
                                            N=length(outcome),
                                            J=length(algos),
                                            threshold_num=100,
                                            thresh_real=100,
                                            hier=hier),iter=800,chains=2,cores=2)
                        
time2 <- Sys.time()

total_time <- difftime(time2,time1)
print(total_time)

# Check loo for hierarchical model
# Loop over all bayes_model functions that match certain tags in archivist

to_loop <- 

loo1 <- get_log_lik(stan_sample=bayes_model,outcome=outcome,algo_data=algos,nwarmup=400,niters=800)

# it seems that some of the observations are very very unlikely -- i.e., the model is almost certainly wrong


improbable <- which(-0.5>loo1[,8])
loo_check1 <- loo(loo1[-improbable,])

loo_data <- loo1 %>% as_data_frame
names(loo_data) <- paste0('Iter',1:length(loo_data))
loo_data <- mutate(loo_data,observation=paste0('Obs_',row_number()))
loo_data <- sample_n(loo_data,100)
loo_data <- gather(loo_data,iteration,log_density,-observation)

loo_data %>% mutate(log_density=exp(log_density)) %>% 
  ggplot(aes(y=log_density,x=observation)) + stat_smooth() + 
  theme_minimal()

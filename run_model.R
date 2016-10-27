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


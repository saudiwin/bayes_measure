# R Script which runs a hiearchical logit and evaluates each algorithm in turn
# Models are evaluated along 3 criteria: binary log loss, ROC curves and PSIS-LOO

require(dplyr)
require(tidyr)
require(magrittr)
require(rstan)
require(shinystan)
require(ggplot2)
require(plotly)
require(archivist)
require(loo)
require(algoeval)

# Set random number seed

seed <- 
set.seed(seed)

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

# Null Model

null_model <- sampling(stan_file,data=list(y=outcome,
                                           x=as.matrix(algos),
                                           N=length(outcome),
                                           J=length(algos),
                                           threshold_num=100,
                                           thresh_real=100,
                                           hier=hier),iter=800,chains=2,cores=2,seed=seed)

# Full model

bayes_model <- sampling(stan_file,data=list(y=outcome,
                                            x=as.matrix(algos),
                                            N=length(outcome),
                                            J=length(algos),
                                            threshold_num=100,
                                            thresh_real=100,
                                            hier=hier),iter=800,chains=2,cores=2,seed=seed)

# Pick columns to do analysis on 

outcols <- sample(400,100)

# Check whether loo works or not

loo1 <- get_log_lik(stan_sample=bayes_model,outcome=outcome,algo_data=algos,nwarmup=400,niters=800)
mean_loo <- colMeans(loo1)
# There is a problem with arithmetic underflow where very certain outcomes of the logit model come out as probability 1,
# Which screws up the loo
check_zeroes <- mean_loo==0
loo1[,check_zeroes] <- loo1[,check_zeroes] + runif(n=400,max=-2*.Machine$double.neg.eps,min=-10*.Machine$double.neg.eps)
#improbable <- -2*sd(mean_loo)>mean_loo
this_loo <- loo(loo1)
bad_pareto <- this_loo$pareto_k>0.7
if(sum(bad_pareto)==0) {
  print('No problematic outliers, proceed as usual')
} else if(sum(bad_pareto)>1) {
  print('There is a problem with the model, need to diagnose before proceeding further. There are some bad Pareto K estimates, but not all.')
}
  # } else if(mean(bad_pareto)>0.5) {
#   print('Extreme outliers detected. Using a dummy variable to control for outliers.')
#   algos <- mutate(algos,outliers=as.numeric(improbable))
#   bayes_model <- sampling(stan_file,data=list(y=outcome,
#                                               x=as.matrix(algos),
#                                               N=length(outcome),
#                                               J=length(algos),
#                                               threshold_num=100,
#                                               thresh_real=100,
#                                               hier=hier),iter=800,chains=2,cores=2)
#   loo1 <- get_log_lik(stan_sample=bayes_model,outcome=outcome,algo_data=algos,nwarmup=400,niters=800)
#   this_loo <- loo(loo1)
#   bad_pareto <- this_loo$pareto_k>0.7
#   if(sum(bad_pareto)==0) {
#     print('Outlier detection worked. PSIS-LOO estimate is stable')
#   } else {
#     print('Outlier detection failed. Model is still mis-specified.')
#   }
# }
  

# Loop over all algos

all_models <- lapply(names(algos), function(x) {
  
  print(paste0('Now running ',x))
  algos <- select(algos,-one_of(x))
  bayes_model <- sampling(stan_file,data=list(y=outcome,
                                            x=as.matrix(algos),
                                            N=length(outcome),
                                            J=length(algos),
                                            threshold_num=100,
                                            thresh_real=100,
                                            hier=hier),iter=800,chains=2,cores=2)
  
  return(bayes_model)

})
all_loos <- lapply(1:length(names(algos)),function(x) {
  this_data <- select(algos,-one_of(names(algos)[x]))
  out_log <- get_log_lik(stan_sample=all_models[[x]],outcome=outcome,algo_data=this_data,nwarmup=400,niters=800)
  mean_loo <- colMeans(out_log)
  # There is a problem with arithmetic underflow where very certain outcomes of the logit model come out as probability 1,
  # Which screws up the loo
  check_zeroes <- (mean_loo==0 | mean_loo>-200*.Machine$double.neg.eps)
  out_log[,check_zeroes] <- out_log[,check_zeroes] + runif(n=400,max=-100*.Machine$double.neg.eps,min=-200*.Machine$double.neg.eps)
  out_loo <- loo(out_log)
  print(out_loo)
  return(out_loo)
})
compare_loos <- lapply(all_loos,function(x) as.numeric(compare(this_loo,x))) %>% unlist %>% 
  matrix(nrow = length(all_loos),ncol=2,byrow=TRUE) %>% as_data_frame
names(compare_loos) <- c('ELPD','SE')
compare_loos %<>% mutate(high_ci=ELPD + 1.96*SE,low_ci=ELPD - 1.96*SE,algos=names(algos))

compare_loos %>% ggplot(aes(y=ELPD,x=reorder(algos,-ELPD))) + geom_errorbar(aes(ymin=low_ci,ymax=high_ci),width=0.2) + 
  geom_point(colour='red',size=2) + theme_minimal() + ylab("Gain in ELPD") + xlab("")

#OK, let's look at binary log-loss
#Certainly easier to calculate than the former

all_rocs <- lapply(1:length(names(algos)),make_roc,outcome=outcome,all_models=all_models,algos=algos,outcols=outcols)

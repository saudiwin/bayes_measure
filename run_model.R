# Run a bayesian measurement model on probability estimates of a binomial outcome (customer purchase)
# V0.1 Robert Kubinec
# 10/18/2016

require(dplyr)
require(tidyr)
require(magrittr)
require(rstan)
require(shinystan)
require(ggplot2)

sample_it <- TRUE

algo_data <- data.table::fread('data/all_combine_200k.csv') %>% as_data_frame

if(sample_it==TRUE) {
  algo_data <- sample_n(algo_data,1000)
}

#Individual binomial trials with probability-related predictors

outcome <- algo_data$realpurch

algos <- select(algo_data,algo1,algo2,algo3,algo4,algo5)

stan_file <- stan_model('measure_v1.stan',model_name='Binomial Bayes Measurement')

bayes_model <- sampling(stan_file,data=list(y=outcome,
                                            x=as.matrix(algos),
                                            N=length(outcome),
                                            J=length(algos),
                                            threshold_num=100,
                                            thresh_real=100),iter=2000,chains=2,cores=2)
                        #control=list(stepsize=0.01, adapt_delta=0.99))

# Do some posterior prediction

roc_curves <- extract(bayes_model,pars='roc_graph')[[1]]
roc_curves <- apply(roc_curves,1,cbind) %>% t %>% as_data_frame
names(roc_curves) <- c(paste0('TruePositive_',1:101),paste0('FalsePositive_',1:101))
roc_curves <- roc_curves %>% mutate(iter=1:n())
roc_curves <- gather(roc_curves,key = "typeofrate",value='Estimate',-iter) %>% separate(col = typeofrate,
                                                                                        into=c('rate','num'),
                                                                                        sep='\\_') %>% 
  spread(key=rate,value=Estimate)

roc_curves %>% ggplot(aes(y=TruePositive,x=FalsePositive,group=iter)) + geom_line()

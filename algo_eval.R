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

seed <- 1827813124
set.seed(seed)

# Pick columns to do analysis on 

outcols <- sample(400,100)

# Data directory for archivist R model collection

repoDir <- 'data/'

# Whether to use samples from the data or the full data

sample_it <- FALSE

# Whether to run a hierarchical model

hier <- TRUE

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



# Do 800 iterations as it is hard to get this much data to compute
# Non-centered parameterization: 9500 seconds with 2 divergent transitions in 800 iterations
# Centered parameterization: 5000 seconds with no divergent transitions, although low ESS on tau
# Centered wins with more data

# Null Model

# stan_file <- stan_model('measure_null.stan',model_name='Binomial Bayes Measurement')
# 
# 
# null_model <- sampling(stan_file,data=list(y=outcome,
#                                            x=as.matrix(algos),
#                                            N=length(outcome),
#                                            J=length(algos),
#                                            threshold_num=100,
#                                            thresh_real=100,
#                                            hier=hier),iter=800,chains=2,cores=2,seed=seed)

# Full model
bayes_model <- readRDS('bayes_model.rds')
# stan_file <- stan_model('measure_v1.stan',model_name='Binomial Bayes Measurement')
# 
# bayes_model <- sampling(stan_file,data=list(y=outcome,
#                                             x=as.matrix(algos),
#                                             N=length(outcome),
#                                             J=length(algos),
#                                             threshold_num=100,
#                                             thresh_real=100,
#                                             hier=hier),iter=800,chains=2,cores=2,seed=seed)

# Check whether loo works or not

loo1 <- get_log_lik(stan_sample=bayes_model,outcome=outcome,algo_data=algos,nwarmup=400,niters=800)
mean_loo <- colMeans(loo1)
# There is a problem with arithmetic underflow where very certain outcomes of the logit model come out as probability 1,
# Which screws up the loo
check_zeroes <- mean_loo==0
loo1[,check_zeroes] <- loo1[,check_zeroes] + runif(n=400,max=-100*.Machine$double.neg.eps,min=-200*.Machine$double.neg.eps)
#improbable <- -2*sd(mean_loo)>mean_loo
this_loo <- loo(loo1)
bad_pareto <- this_loo$pareto_k>0.7
if(sum(bad_pareto,na.rm=TRUE)==0 && sum(is.na(bad_pareto))==0) {
  print('No problematic outliers, proceed as usual')
} else {
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
# Use a non-hierarchical logit (much faster!)
stan_file <- stan_model('measure_singleton.stan',model_name='Binomial Bayes Measurement')
all_models <- lapply(names(algos), function(x) {
  
  print(paste0('Now running ',x))
  algos <- select(algos,one_of(x))
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
  this_data <- select(algos,one_of(names(algos)[x]))
  out_log <- get_log_lik(stan_sample=all_models[[x]],outcome=outcome,algo_data=this_data,nwarmup=400,niters=800)
  mean_loo <- colMeans(out_log)
  # There is a problem with arithmetic underflow where very certain outcomes of the logit model come out as probability 1,
  # Which screws up the loo
  check_zeroes <- (mean_loo==0 | mean_loo>-200*.Machine$double.neg.eps)
  out_log[,check_zeroes] <- out_log[,check_zeroes] + runif(n=400,max=-100*.Machine$double.neg.eps,min=-200*.Machine$double.neg.eps)
  check_neg_inf <- mean_loo < log(100*.Machine$double.neg.eps)
  # It seems like smallest number you can log is 1e-300 (an infitesimal probability)
  #see if there is a number besides infinity
  # get the lowest non-infinite number (mean and sd)
  min_mean <- min(mean_loo[!is.infinite(mean_loo)])
  min_sd <- sd(out_log[,which(mean_loo==min_mean)])
  out_log[,check_neg_inf] <- apply(out_log[,check_neg_inf,drop=FALSE],2,function(x) {
      x <- log(runif(n=400,min=100*.Machine$double.neg.eps,max=200*.Machine$double.neg.eps))
    return(x)
  })
  # out_log[,check_neg_inf] <- ifelse(is.infinite(out_log[,check_neg_inf]),0,out_log[,check_neg_inf])
  # out_log[,check_neg_inf] <- out_log[,check_neg_inf] + runif(n=400,min=-(40+100*.Machine$double.neg.eps),
  #                                                            max= -(40-100*.Machine$double.neg.eps) )
  out_loo <- loo(out_log)
  print(out_loo)
  return(out_loo)
})
compare_loos <- lapply(all_loos,function(x) as.numeric(compare(this_loo,x))) %>% unlist %>% 
  matrix(nrow = length(all_loos),ncol=2,byrow=TRUE) %>% as_data_frame
names(compare_loos) <- c('ELPD','SE')
compare_loos %<>% mutate(high_ci=ELPD + 1.96*SE,low_ci=ELPD - 1.96*SE,algos=names(algos))
  
compare_loos %>% ggplot(aes(y=ELPD,x=reorder(algos,-ELPD))) + geom_errorbar(aes(ymin=low_ci,ymax=high_ci),width=0.2) + 
  geom_point(colour='red',size=2) + theme_minimal() + ylab("Performance Decrease Compared to Full Model") + xlab("")

#OK, let's look at binary log-loss
#Certainly easier to calculate than the former

full_roc <- make_roc(algo_num = 'all',
                     outcome=outcome,models=bayes_model,
                     algos=algos,outcols=outcols,threshold=1000) %>% select(-Algorithm)

full_roc %<>% group_by(Threshold) %>% summarize(mean_true=mean(True_Positives_Rate),
                        mean_false=mean(False_Positive_Rate),
                        high_ci=quantile(True_Positives_Rate,probs=.95),
                        low_ci=quantile(True_Positives_Rate,probs=.05)) %>% ungroup

all_rocs <- lapply(1:length(names(algos)),make_roc,outcome=outcome,models=all_models,algos=algos,outcols=outcols,
                   threshold=1000)

all_rocs <- all_rocs %>% bind_rows

all_rocs <- all_rocs %>% group_by(Algorithm,Threshold) %>% 
  summarize(mean_true=mean(True_Positives_Rate),
            mean_false=mean(False_Positive_Rate),
            high_ci=quantile(True_Positives_Rate,probs=.95),
            low_ci=quantile(True_Positives_Rate,probs=.05)) %>% ungroup

all_rocs %>% ggplot(aes(y=mean_true,x=mean_false)) + facet_wrap(~Algorithm) +
   geom_ribbon(data=full_roc,aes(ymin=low_ci,ymax=high_ci),fill='red',colour='red',alpha=0.5) + geom_ribbon(aes(ymin=low_ci,ymax=high_ci),fill='blue',colour='blue',alpha=0.5) +
  theme_minimal()


full_logloss <- binary_log_loss(algo_num='all',outcome=outcome,models=bayes_model,algos=algos)


names(full_logloss) <- paste0("Iter_",401:800)
full_logloss$obs_num <- paste0("Obs_",1:nrow(full_logloss))

#Let's sample observations, otherwise it gets unwieldy

sample_obs <- sample(1:nrow(full_logloss),1000)

full_logloss <- full_logloss %>% slice(sample_obs) %>% gather(iter_number,logloss,-obs_num)

logloss_obs_full <- full_logloss %>% group_by(obs_num) %>% 
  summarize(mean_loss = mean(logloss),low_ci=quantile(logloss,probs = 0.05),
            high_ci=quantile(logloss,probs=0.95)) %>% ungroup 


all_logloss <- lapply(1:length(names(algos)),binary_log_loss,outcome=outcome,models=all_models,algos=algos)

all_logloss <- lapply(all_logloss,as_data_frame)
all_logloss <- lapply(all_logloss,slice,sample_obs)
all_logloss <- bind_rows(all_logloss)

names(all_logloss) <- paste0("Iter_",401:800)
all_logloss$obs_num <- rep(paste0("Obs_",sample_obs),length(algos))
all_logloss$type <- paste0("Algo_",rep(1:length(algos),each=length(sample_obs)))

#Let's sample observations, otherwise it gets unwieldy

all_logloss <- all_logloss %>% gather(iter_number,logloss,-obs_num,-type)

logloss_obs_all <- all_logloss %>% group_by(type,obs_num) %>% 
  summarize(mean_loss = mean(logloss),low_ci=quantile(logloss,probs = 0.05),
            high_ci=quantile(logloss,probs=0.95)) %>% ungroup 

logloss_obs_all %>% ggplot(aes(x=reorder(obs_num,mean_loss))) + geom_ribbon(data=logloss_obs_full,aes(ymin=low_ci,ymax=high_ci),colour="red") +
  geom_ribbon(aes(ymin=low_ci,ymax=high_ci),colour='blue') + facet_wrap(~type) + xlab("") +
  ylab("Total Prediction Gain") + theme_minimal() + theme(axis.text.x=element_blank(),plot.background = element_blank(),
                                               panel.background=element_blank()) +ylim(0,20)


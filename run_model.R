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

sample_it <- TRUE

algo_data <- data.table::fread('data/all_combine_200k.csv') %>% as_data_frame

if(sample_it==TRUE) {
  algo_data <- sample_n(algo_data,1000)
}

#Individual binomial trials with probability-related predictors

outcome <- algo_data$realpurch

algos <- select(algo_data,algo1,algo2,algo3,algo4,algo5) %>% mutate_all(scale) 

algos %>% gather(algo_type,prob_score) %>% ggplot(aes(x=prob_score)) + geom_histogram(bins=60) +
  facet_wrap(~algo_type,scales='free')

stan_file <- stan_model('measure_v1.stan',model_name='Binomial Bayes Measurement')

bayes_model <- sampling(stan_file,data=list(y=outcome,
                                            x=as.matrix(algos),
                                            N=length(outcome),
                                            J=length(algos),
                                            threshold_num=100,
                                            thresh_real=100),iter=2000,chains=2,cores=2)
                        #control=list(stepsize=0.01, adapt_delta=0.99))

# Does this model really have an 80% positive prediction rate with no false positives? If so, that's pretty darn good

replicates <- extract(bayes_model,pars='pred_success')[[1]] %>% t %>% as_data_frame %>% sample(100)
names(replicates) <- paste0('iter_',1:length(names(replicates)))
replicates$obs_num <- paste0('Obs_',1:1000)
replicates$original <- outcome
replicates %<>% gather(iterates,iter_num,-obs_num,-original) %>% group_by(obs_num) %>% summarize(num_right=sum(iter_num==original),
                                                                                                 num_wrong=sum(iter_num!=original)) 

Sys.setenv("plotly_username" = "bobkubinec")
Sys.setenv("plotly_api_key" = "8q00qm53km")

output <- replicates %>% plot_ly(x=~obs_num) %>% add_markers(y=~num_right,color=I('red')) %>% layout(yaxis=list(title='Number Correct Predicted'),
                                                                                           xaxis=list(showgrid=F,title="Observed Outcomes (1 or 0)",
                                                                                                      showticklabels=F))
plotly_POST(output)
roc_curves <- extract(bayes_model,pars='roc_graph')[[1]]
roc_curves <- apply(roc_curves,1,cbind) %>% t %>% as_data_frame
names(roc_curves) <- c(paste0('TruePositive_',1:101),paste0('FalsePositive_',1:101))
roc_curves <- roc_curves %>% mutate(iter=1:n())
roc_curves <- gather(roc_curves,key = "typeofrate",value='Estimate',-iter) %>% separate(col = typeofrate,
                                                                                        into=c('rate','num'),
                                                                                        sep='\\_') %>% 
  spread(key=rate,value=Estimate)

# # Add in ROC curve from existing data (algorithms)
# 
make_roc <- function(algo_num=NULL,outcome=NULL,threshold=101) {
  algo <- algos[[algo_num]]
  steps <- seq(from=0,to=1,length.out=threshold)
  roc_data <- matrix(nrow=threshold,ncol=2)
  num_pos <- sum(outcome)
  num_neg <- length(outcome) - sum(outcome)
  for(i in 1:threshold) {
    tp <- (plogis(algo[outcome==1])>steps[i]) %>% sum
    fp <- (plogis(algo[outcome==0])>steps[i]) %>% sum
    roc_data[i,] <- c(tp/num_pos,fp/num_neg)
  }
  roc_data <- roc_data %>%  as_data_frame
  roc_data$algo_type <- paste0("algo_",algo_num)
  return(roc_data)
}


roc_out <-  lapply(1:length(algos),make_roc,outcome) %>% bind_rows

roc_curves %>% ggplot(aes(y=TruePositive,x=FalsePositive,group=iter)) + geom_point(alpha=0.5,colour='grey50') + geom_abline(slope=1,
                                                                                                           intercept=0) +
  geom_line(data=roc_out,aes(y=V1,x=V2,group=algo_type,colour=algo_type),size=1) + theme_minimal() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) + annotate('rect',xmax=0.18,xmin=0,ymin=0.91,ymax=0.97,alpha=0.1,colour='pink') +
  geom_label(data=roc_out,aes(y=V1,x=V2,group=algo_type,label=algo_type))





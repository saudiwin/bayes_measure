
require(dplyr)
require(tidyr)
require(magrittr)
require(rstan)
require(shinystan)
require(ggplot2)
require(plotly)
require(archivist)
require(loo)


algo_data <- data.table::fread('data/all_combine_200k.csv') %>% as_data_frame

outcome <- algo_data$realpurch

# Standardize algorithmic output to the same scale

algos <- select(algo_data,algo1,algo2,algo3,algo4,algo5) %>% mutate_all(scale) 

# Need to work on computing binary log loss so we can figure out which observations are problematic and why


base_model <- loadFromLocalRepo('41e3e006ed220f3601f71585b19fdab6',repoDir='data/',value=TRUE)

all_logloss <- binary_log_loss(base_model,outcome=outcome,algo_data=algos,nwarmup=400,niters=800) %>% as_data_frame

names(all_logloss) <- paste0("Iter_",401:800)
all_logloss$obs_num <- paste0("Obs_",1:nrow(all_logloss))
all_logloss <- all_logloss %>% gather(iter_number,logloss,-obs_num)

# Collapse dataset down to mean and sd of logloss by each observation for plotting/summarization

logloss_obs <- all_logloss %>% group_by(obs_num) %>% 
  summarize(mean_loss = mean(logloss),low_ci=quantile(logloss,probs = 0.05),
            high_ci=quantile(logloss,probs=0.95)) %>% ungroup 

logloss_obs %>% ggplot(aes(y=mean_loss,x=reorder(obs_num,-mean_loss),ymax=high_ci,ymin=low_ci)) + geom_pointrange(colour='grey50') +
  xlab("Observations") + ylab("Log Loss Over True Outcome") + theme(axis.text.x=element_blank(),
                   axis.ticks.x=element_blank())

# plot tau against mean logloss per iter

tau <- extract(base_model,'tau')[[1]][401:800] %>% as_data_frame %>% mutate(iter_number=paste0('Iter_',401:800))
logloss_iter <- all_logloss %>% group_by(iter_number) %>% summarize(mean_loss=mean(logloss)) %>% ungroup

logloss_iter %>% left_join(tau) %>% ggplot(aes(y=mean_loss,x=value)) + geom_point() + ylab("Mean Logloss") +
  xlab("Values of Tau (Between Effect)") + stat_smooth()+ theme_minimal()

# plot tau against mean logloss for bottom 1 percent 

bottom1percent <- filter(logloss_obs,mean_loss<percent_rank(1))
filter(all_logloss,obs_num %in% bottom1percent$obs_num) %>% group_by(iter_number) %>% 
  summarize(mean_loss=mean(logloss)) %>% left_join(tau) %>% ggplot(aes(y=mean_loss,x=value)) + geom_point() + ylab("Mean Logloss Bottom 1 Percent") +
  xlab("Values of Tau (Between Effect)") + stat_smooth()+ theme_minimal()
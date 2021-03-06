# We can compute ranks of the different algorithms by drawing from the hyperparameters
# In particular, we can examine how often algorithm 5, the highest-performing algorithm, has a higher prob. value
# Than the combined model [any theta_j]


# Does this model really have an 80% positive prediction rate with no false positives? If so, that's pretty darn good

replicates <- extract(bayes_model,pars='pred_success')[[1]] %>% t %>% as_data_frame %>% sample(200)
names(replicates) <- paste0('iter_',1:length(names(replicates)))
replicates$obs_num <- paste0('Obs_',1:1000)
replicates$original <- outcome
replicates %<>% gather(iterates,iter_num,-obs_num,-original) %>% group_by(obs_num) %>% summarize(num_right=sum(iter_num==original),
                                                                                                 num_wrong=sum(iter_num!=original)) 

Sys.setenv("plotly_username" = "bobkubinec")
Sys.setenv("plotly_api_key" = "8q00qm53km")

# Output raw predictive performance

output <- replicates %>% plot_ly(x=~obs_num) %>% add_markers(y=~num_right,color=I('red')) %>% layout(yaxis=list(title='Number Correct Predicted'),
                                                                                                     xaxis=list(showgrid=F,title="Observed Outcomes (1 or 0)",
                                                                                                                showticklabels=F))
plotly_POST(output)

# Look at ROC curve differences between models and algorithm 


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
make_roc <- function(algo_num=NULL,outcome=NULL,threshold=10000) {
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

#Put labels where the lines are closest to 50% false positives

roc_out <- roc_out %>% group_by(algo_type) %>% mutate(from_mean = (0.5 - V1)^2 + (0.5-V2)^2,
                                                      line_labels = ifelse(from_mean==min(from_mean),algo_type,NA)) %>% 
  ungroup %>% arrange(V1,V2)

roc_curves %>% ggplot(aes(y=TruePositive,x=FalsePositive,group=iter)) + geom_point(alpha=0.2,colour='grey50') + geom_abline(slope=1,
                                                                                                                            intercept=0) +
  geom_path(data=roc_out,aes(y=V1,x=V2,group=algo_type,colour=algo_type),size=0.8) + theme_minimal() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) + annotate('rect',xmax=0.10,xmin=0,ymin=0.91,ymax=0.97,alpha=0.1,colour='pink') +
  geom_text(data=roc_out,aes(y=V1,x=V2,group=algo_type,label=line_labels),hjust='inward',vjust='inward',check_overlap=TRUE) +
  scale_color_brewer(palette='Set1')

# plotly_out <- roc_curves  %>% plot_ly(x=~FalsePositive,y=~TruePositive,type='scatter') %>% add_markers()
# 
# plotly_POST(plotly_out)

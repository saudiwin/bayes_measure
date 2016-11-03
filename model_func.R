# Helper functions for model evaluation

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

get_log_lik <- function(stan_sample=NULL,outcome=NULL,algo_data=NULL,nwarmup=NULL,
                        niters=NULL) {
  
  predictors <- rstan::extract(stan_sample,pars='theta_raw')[[1]][(nwarmup+1):niters,]
  intercept <- rstan::extract(stan_sample,pars='alpha')[[1]][(nwarmup+1):niters]
  algo_data <- as.matrix(algo_data)
  
  raw_predict <- algo_data %*% t(predictors)
  rm(predictors)
  raw_predict <- apply(raw_predict,1,function(x) x + intercept) %>% t
  # apply loop is too memory hungry
  # switch to data.table to modify one column at a time in-place
  raw_predict <- data.table::as.data.table(raw_predict)
  
  #log_lik <- apply(raw_predict,2,function(x) dbinom(x=outcome,size=1,prob=plogis(x)))
  #increases memory by about double, but it still able to return the object
  #runs much faster than apply loop
  raw_predict <- raw_predict[,lapply(.SD,
                                 function(x) dbinom(x=outcome,size=1,prob=plogis(x),log=TRUE))]
  
  return(as.matrix(raw_predict))
}

binary_log_loss <- function(stan_sample=NULL,outcome=NULL,algo_data=NULL,nwarmup=NULL,
                            niters=NULL) {
  
  predictors <- rstan::extract(stan_sample,pars='theta_raw')[[1]][(nwarmup+1):niters,]
  intercept <- rstan::extract(stan_sample,pars='alpha')[[1]][(nwarmup+1):niters]
  algo_data <- as.matrix(algo_data)
  
  raw_predict <- algo_data %*% t(predictors)
  rm(predictors)
  raw_predict <- apply(raw_predict,1,function(x) x + intercept) %>% t
  # apply loop is too memory hungry
  # switch to data.table to modify one column at a time in-place
  raw_predict <- data.table::as.data.table(raw_predict)
  
  #calculate predicted probabilities per observation
  raw_predict <- raw_predict[,lapply(.SD,plogis)]
  
  log_loss <- raw_predict[,lapply(.SD,function(x) {
    output <- ifelse(outcome==0,-log(x),-log(1-x))
    # Overflow zero prediction errors to 50
    output[is.infinite(output)] <- 50
    return(output)
    }
    )]
  
  return(log_loss)
  
}
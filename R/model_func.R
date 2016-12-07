# Helper functions for model evaluation

# # Add in ROC curve from existing data (algorithms)
#' @export
#' @import dplyr
#' @import tidyr
make_roc <- function(algo_num='all',outcome=NULL,models=NULL,algos=NULL,threshold=1000,
                     nwarmup=400,niters=800,outcols=NULL) {
  
  if(algo_num!='all') {
    print(paste0('Now processing algorithm ',names(algos)[algo_num]))
    algo_data <-  select(algos,one_of(names(algos)[algo_num])) %>% as.matrix
    stan_sample <- models[[algo_num]]
    
    predictors <- rstan::extract(stan_sample,pars='theta_raw')[[1]][(nwarmup+1):niters,]
    intercept <- rstan::extract(stan_sample,pars='alpha')[[1]][(nwarmup+1):niters]
    predictors <- predictors[outcols]
    intercept <- intercept[outcols]
    
    # Convert predictors to matrix to enable matrix multiplication
    # then it is conformable if transposed (NX1) x (1XS) where N=data and S=iterations
     
    raw_predict <- algo_data %*% t(as.matrix(predictors))
    rm(predictors)
    raw_predict <- apply(raw_predict,1,function(x) x + intercept) %>% t
    # apply loop is too memory hungry
    # switch to data.table to modify one column at a time in-place
    raw_predict <- data.table::as.data.table(raw_predict)
    
  } else {
    stan_sample <- models
    algo_data <- as.matrix(algos)
    predictors <- rstan::extract(stan_sample,pars='theta_raw')[[1]][(nwarmup+1):niters,]
    intercept <- rstan::extract(stan_sample,pars='alpha')[[1]][(nwarmup+1):niters]
    predictors <- predictors[outcols,]
    intercept <- intercept[outcols]
    
    raw_predict <- algo_data %*% t(predictors)
    rm(predictors)
    raw_predict <- apply(raw_predict,1,function(x) x + intercept) %>% t
    # apply loop is too memory hungry
    # switch to data.table to modify one column at a time in-place
    raw_predict <- data.table::as.data.table(raw_predict)
  }
  
  
  #calculate predicted probabilities per observation
  raw_predict <- raw_predict[,lapply(.SD,plogis)]
  
  steps <- seq(from=0,to=1,length.out=threshold)
  roc_data <- matrix(nrow=threshold,ncol=2)
  num_pos <- sum(outcome)
  num_neg <- length(outcome) - sum(outcome)
  col <- 1
  rocs <- lapply(raw_predict,function(x) {
    roc_data <- matrix(nrow=threshold,ncol=2)
    for(i in 1:threshold) {
      tp <- (x[outcome==1]>steps[i]) %>% sum
      fp <- (x[outcome==0]>steps[i]) %>% sum
      roc_data[i,] <- c(tp/num_pos,fp/num_neg)
    }
    roc_data <- roc_data %>%  as_data_frame
    roc_data <- mutate(roc_data,algo_type= paste0("algo_",algo_num),step=steps,iteration=paste0('iteration_',col))
    print(paste0('Finished iteration ',col))
    col <<- col + 1
    return(roc_data)
  }) %>% bind_rows
  names(rocs) <- c('True_Positives_Rate','False_Positive_Rate','Algorithm','Threshold','Iteration')

  return(rocs)
}

#' @import data.table
#' @export
get_log_lik <- function(stan_sample=NULL,outcome=NULL,algo_data=NULL,nwarmup=NULL,
                        niters=NULL) {
  
  predictors <- rstan::extract(stan_sample,pars='theta_raw')[[1]][(nwarmup+1):niters,]
  intercept <- rstan::extract(stan_sample,pars='alpha')[[1]][(nwarmup+1):niters]
  algo_data <- select(algo_data,matches('algo')) %>% as.matrix
  raw_predict <- algo_data %*% t(predictors)
  rm(predictors)
  raw_predict <- apply(raw_predict,1,function(x) x + intercept) %>% t
  # apply loop is too memory hungry
  # switch to data.table to modify one column at a time in-place
  raw_predict <- as.data.table(raw_predict)
  
  #log_lik <- apply(raw_predict,2,function(x) dbinom(x=outcome,size=1,prob=plogis(x)))
  #increases memory by about double, but it still able to return the object
  #runs much faster than apply loop
  raw_predict <- raw_predict[,lapply(.SD,
                                 function(x) dbinom(x=outcome,size=1,prob=plogis(x),log=TRUE))]
  
  return(t(as.matrix(raw_predict)))
}

#' @export
binary_log_loss <- function(algo_num='all',outcome=NULL,models=NULL,algos=NULL,
                            nwarmup=400,niters=800) {
  if(algo_num!='all') {
    print(paste0('Now processing algorithm ',names(algos)[algo_num]))
    algo_data <-  select(algos,one_of(names(algos)[algo_num])) %>% as.matrix
    stan_sample <- models[[algo_num]]
    
    predictors <- rstan::extract(stan_sample,pars='theta_raw')[[1]][(nwarmup+1):niters]
    intercept <- rstan::extract(stan_sample,pars='alpha')[[1]][(nwarmup+1):niters]
    algo_data <- as.matrix(algo_data)
    
    raw_predict <- algo_data %*% t(as.matrix(predictors))
    rm(predictors)
    raw_predict <- apply(raw_predict,1,function(x) x + intercept) %>% t  
  } else {
    stan_sample <- models
    predictors <- rstan::extract(stan_sample,pars='theta_raw')[[1]][(nwarmup+1):niters,]
    intercept <- rstan::extract(stan_sample,pars='alpha')[[1]][(nwarmup+1):niters]
    algo_data <- as.matrix(algos)
    
    raw_predict <- algo_data %*% t(predictors)
    rm(predictors)
    raw_predict <- apply(raw_predict,1,function(x) x + intercept) %>% t
  }
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


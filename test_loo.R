# Test psis-loo calc

beta1 <- 1.1
X <- rnorm(100)
y <- rbinom(100,1,prob=plogis(beta1 * X))


# instead of beta1, we have random estimates of beta1

beta_est <- rnorm(100,beta1,0.2)

log_lik <- sapply(beta_est,function(beta) {
  
  log_est <- dbinom(y,size=1,prob=plogis(beta*X),log=TRUE)
})

check_loo1 <- loo(log_lik)

beta_est <- rnorm(100,beta1,1)

log_lik <- sapply(beta_est,function(beta) {
  
  log_est <- dbinom(y,size=1,prob=plogis(beta*X),log=TRUE)
})

check_loo2 <- loo(log_lik)

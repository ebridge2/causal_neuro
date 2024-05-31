require(copula)
require(causalBatch)

expit <- function(x) {
  1/(1 + exp(-x))
}

collider.sim_simple <- function(n=200, predictor_prob=0.25, eff_sz=1, a=2,
                                corruption=.2, confounding=0) {
  risk_outcome <- rbinom(n, 1, prob=predictor_prob)
  # confounder
  xs <- 2*rbeta(n, 5, 5) - 1
  
  # true measurement edge-weight
  #t_true <- (1 - confounding)*(2*sapply(1:n, function(i) {
  #  ifelse(risk_outcome[i]==1, rbeta(1, 10, 2), rbeta(1, 15, 15))
  #}) - 1) + confounding*xs
  # true measurement edge-weight
  tx <- generate_correlated_betas_pergroup(risk_outcome, alphas1=c(4, 5), betas1=c(4, 5),
                                          alphas2=c(4, 5), betas2=c(2, 5), corr.mtx=confounding)
  t_true <- tx[,1]; xs <- tx[,2]
  
  # whether or not someone has outcome
  p_y <- expit(eff_sz*(t_true) + a*xs - 2.5)
  ys <- rbinom(n, 1, p_y)
  
  # usability degree
  p_usable <- 1 - expit(t_true + 2*xs + 2*ys - 2.5)
  usable <- as.logical(rbinom(n, 1, p_usable))
  
  # bias measurements due to unusability
  t_obs <- (1 - corruption)*t_true + corruption*(sapply(1:n, function(i) {
    ifelse(usable[i], rnorm(1, mean=0, sd=0.5), rnorm(1, mean=0.3, sd=0.5))
  }))
  
  return(list(Ys=ys, Ts=t_obs, Xs=xs, Risk.group=risk_outcome, usable=usable, p_usable=p_usable, Ttrue=t_true,
              Effect=eff_sz))
}

collider.sim_linear <- function(n=200, predictor_prob=0.25, eff_sz=1, covar_eff=2, corruption=.2, confounding=0) {
  risk_outcome <- rbinom(n, 1, prob=predictor_prob)
  # confounder
  xs <- 2*rbeta(n, 5, 5) - 1
  
  t_true <- 2*sapply(1:n, function(i) {
    ifelse(risk_outcome[i]==1, rbeta(1, + 4*(1 + confounding*xs[i]), 2),
           rbeta(1, 4*(1 + confounding*xs[i]), 4*(1 + confounding*xs[i])))
  }) - 1
  
  # whether or not someone has outcome
  p_y <- expit(eff_sz*(t_true) + covar_eff*xs - 1.5)
  ys <- rbinom(n, 1, p_y)
  
  # usability degree
  p_usable <- 1 - expit(t_true + 2*xs + 2*ys - 1.5)
  usable <- as.logical(rbinom(n, 1, p_usable))
  
  # bias measurements due to unusability
  t_obs <- (1 - corruption)*t_true + corruption*(2*sapply(1:n, function(i) {
    ifelse(usable[i], rbeta(1, 10, 10), rbeta(1, 10, 2))
  })-1)
  
  return(list(Ys=ys, Ts=t_obs, Xs=xs, Risk.group=risk_outcome, usable=usable, p_usable=p_usable, Ttrue=t_true,
              Effect=eff_sz))
}

collider.sim_nonlinear <- function(n=200, predictor_prob=0.25, eff_sz=1, covar_eff=2, b=2, corruption=.2, confounding=0) {
  risk_outcome <- rbinom(n, 1, prob=predictor_prob)
  # confounder
  xs <- 2*rbeta(n, 5, 5) - 1
  
  t_true <- 2*sapply(1:n, function(i) {
    ifelse(risk_outcome[i]==1, rbeta(1, + 4*(1 + confounding*xs[i]), 2),
           rbeta(1, 4*(1 + confounding*xs[i]), 4*(1 + confounding*xs[i])))
  }) - 1
  
  # whether or not someone has outcome
  p_y <- expit(eff_sz*t_true + covar_eff*causalBatch:::sigmoid(b*xs) - 1.5)
  ys <- rbinom(n, 1, p_y)
  
  # usability degree
  p_usable <- 1 - expit(t_true + 2*xs + 2*ys - 1.5)
  usable <- as.logical(rbinom(n, 1, p_usable))
  
  # bias measurements due to unusability
  t_obs <- (1 - corruption)*t_true + corruption*(2*sapply(1:n, function(i) {
    ifelse(usable[i], rbeta(1, 10, 10), rbeta(1, 10, 2))
  })-1)
  
  return(list(Ys=ys, Ts=t_obs, Xs=xs, Risk.group=risk_outcome, usable=usable, p_usable=p_usable, Ttrue=t_true,
              Effect=eff_sz))
}

collider.sim_nonmonotone <- function(n=200, predictor_prob=0.25, eff_sz=1, covar_eff=2, corruption=.2, confounding=0) {
  risk_outcome <- rbinom(n, 1, prob=predictor_prob)
  # confounder
  xs <- 2*rbeta(n, 5, 5) - 1
  
  t_true <- 2*sapply(1:n, function(i) {
    ifelse(risk_outcome[i]==1, rbeta(1, + 4*(1 + confounding*xs[i]), 2),
           rbeta(1, 4*(1 + confounding*xs[i]), 4*(1 + confounding*xs[i])))
  }) - 1
  
  # whether or not someone has outcome
  p_y <- expit(eff_sz*t_true + covar_eff*dnorm(xs, mean=.5, sd=1/2) - 1.5)
  ys <- rbinom(n, 1, p_y)
  # usability degree
  p_usable <- 1 - expit(t_true + 2*xs + 2*ys - 1.5)
  usable <- as.logical(rbinom(n, 1, p_usable))
  
  # bias measurements due to unusability
  t_obs <- (1 - corruption)*t_true + corruption*(2*sapply(1:n, function(i) {
    ifelse(usable[i], rbeta(1, 10, 10), rbeta(1, 10, 2))
  })-1)
  
  return(list(Ys=ys, Ts=t_obs, Xs=xs, Risk.group=risk_outcome, usable=usable, p_usable=p_usable, Ttrue=t_true,
              Effect=eff_sz))
}
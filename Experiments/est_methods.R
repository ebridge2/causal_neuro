require(survival)
require(twang)
require(survey)
require(MASS)
require(causalBatch)

apply_poisson <- function(Ys, Ts, Xs, models) {
  retain.ids <- which(Xs >= min(models$xlims) & Xs <= max(models$xlims))
  Y.prime <- Ys[retain.ids,,drop=FALSE]; T.prime <- Ts[retain.ids]; X.prime <- Xs[retain.ids]
  Y.adj <- matrix(NA, nrow=nrow(Y.prime), ncol=ncol(Y.prime))
  for (k in 1:length(models$fit)) {
    adj.coef <- sapply(T.prime, function(batch) {
      ifelse(batch == 1, exp(models$fit[[k]]$coefficients["Batch1"]), 1)
    })
    Y.adj[,k] <- Y.prime/adj.coef
  }
  return(list(Ys=Y.adj, Ts=T.prime, Xs=X.prime))
}

ass.pois <- function(Ys, Ts, Xs, ...) {
  mod <- list()
  for (k in 1:dim(Ys)[2]) {
    dat <- data.frame(Y=Ys[,k], Batch=Ts, Covar=Xs)
    mod[[k]] <- glm(Y ~ Batch, family=poisson, data=dat)
  }
  return(list(fit=mod, apply=apply_poisson, xlims=c(min(Xs), max(Xs))))
}

cond.pois <- function(Ys, Ts, Xs, ...) {
  mod <- list()
  for (k in 1:dim(Ys)[2]) {
    dat <- data.frame(Y=Ys[,k], Batch=Ts, Covar=Xs)
    mod[[k]] <- glm(Y ~ Batch + Covar, family=poisson, data=dat)
  }
  return(list(fit=mod, apply=apply_poisson, xlims=c(min(Xs), max(Xs))))
}

sigmoid.xfm <- function(Xs, a=4, b=8) {
  sapply(Xs, function(x) {
    a*causalBatch:::sigmoid(b*x)
  })
}

linear.xfm <- function(Xs, a=2, b=1) {
  return(a*Xs + b)
}

impulse.xfm <- function(Xs, a=1/2, b=1/2, c=4) {
  return(dnorm(Xs, mean=a, sd=b)*c)
}

oracle.pois <- function(Ys, Ts, Xs, xfm=c(), xfm.params=list(), ...) {
  Xs.oracle <- do.call(xfm, c(list(apply(Xs, 1, mean)), xfm.params))
  
  mod <- list()
  for (k in 1:dim(Ys)[2]) {
    dat <- data.frame(Y=Ys[,k], Batch=Ts, Covar=Xs.oracle)
    mod[[k]] <- glm(Y ~ Batch + Covar, family=poisson, data=dat)
  }
  return(list(fit=mod, apply=apply_poisson, xlims=c(min(Xs, max(Xs)))))
}

strat.pois <- function(Ys, Ts, Xs, nstrata=5, ...) {
  strata.cuts <- seq(min(Xs), max(Xs), length.out=nstrata)
  Xs.strata <- cut(Xs, breaks=strata.cuts)
  
  mod <- list()
  for (k in 1:dim(Ys)[2]) {
    dat <- data.frame(Y=Ys[,k], Batch=Ts, Covar=Xs.strata)
    mod[[k]] <- glm(Y ~ Batch + strata(Covar), family=poisson, data=dat)
  }
  return(list(fit=mod, apply=apply_poisson, xlims=c(min(Xs), max(Xs)), is.stratified=TRUE,
              strata.cuts=strata.cuts))
}

ipw.pois <- function(Ys, Ts, Xs, ...) {
  retain.ids <- causalBatch::cb.align.vm_trim(Ts, Xs)
  Ys.prime <- Ys[retain.ids,,drop=FALSE]; Ts.prime <- Ts[retain.ids]
  Xs.prime <- Xs[retain.ids]
  
  dat <- data.frame(Batch=Ts.prime, Covar=Xs.prime)
  # get propensity scores and weights
  ps_mod <- glm("Batch ~ Covar", family=binomial, data=dat)
  probs <- predict(ps_mod, dat, type="response")
  dat$wts <- sapply(1:length(Ts.prime), function(i) {
    if (Ts.prime[i] == 1) {
      1/probs[i]
    } else {
      1/(1 - probs[i])
    }
  })
  
  mod <- list()
  for (k in 1:dim(Ys)[2]) {
    dat$Y <- Ys.prime[,k]
    mod[[k]] <- svyglm(Y ~ Batch, design=svydesign(id=~1, weights=~wts, data=dat), 
                  data=dat, family=poisson)
  }
  return(list(fit=mod, apply=apply_poisson, xlims=c(min(Xs.prime), max(Xs.prime))))
}

match.pois <- function(Ys, Ts, Xs, ...) {
  retain.ids.pt <- causalBatch::cb.align.vm_trim(Ts, Xs)
  Ys.prime <- Ys[retain.ids.pt,,drop=FALSE]; Ts.prime <- Ts[retain.ids.pt]
  Xs.prime <- Xs[retain.ids.pt]
  
  retain.ids.match <- causalBatch::cb.align.kway_match(Ts.prime, data.frame(Xs.prime),
                                          match.form = "Xs.prime",
                                          match.args = list(method="nearest", distance="mahalanobis", 
                                                            caliper=c(Xs.prime=.2)))
  Ys.prime <- Ys.prime[retain.ids.match,,drop=FALSE]; Ts.prime <- Ts.prime[retain.ids.match]
  Xs.prime <- Xs.prime[retain.ids.match]
  
  for (k in 1:dim(Ys)[2]) {
    dat <- data.frame(Y=Ys.prime[,k], Batch=Ts.prime, Covar=Xs.prime)
    mod[[k]] <- glm(Y ~ Batch, family=poisson, data=dat)
  }
  return(list(fit=mod, apply=apply_poisson, xlims=c(min(Xs.prime), max(Xs.prime))))
}


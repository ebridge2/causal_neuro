require(causalBatch)
source("./utilities.R")

sim_sigmoid <- function(n=100, pi=.5, eff_sz=1, alpha=2, d=2, unbalancedness=1, 
                        a=4, b=8, err=1/2, unobs.cor=NULL, ncovar=1, unobs.alpha=1,
                        unobs.beta=3, nbreaks=200) {
  batches <- rbinom(n=n, size=1, prob=pi)
  
  beta <- alpha*unbalancedness
  if (!is.null(unobs.cor)) {
    corr.mtx <- make_corr_mtx(ncovar,unobs.cor)
    alphas <- c(alpha, rep(2, ncovar - 1), unobs.alpha); betas <- c(beta, rep(2, ncovar - 1), unobs.beta)
  } else {
    corr.mtx = diag(ncovar)
    alphas <- alpha; betas <- beta
  }
  Xs <- generate_correlated_betas_pergroup(batches, alphas, betas, betas, alphas, corr.mtx=corr.mtx)
  scale_off <- seq(1/d, 1/d^2, length.out=d)
  if (dim(Xs)[2] == 1) {
    covar.offset <- 1
  } else {
    covar.offset <- c(1, rep(1, dim(Xs)[2] - 1))
  }
  ys <- sapply(1:d, function(k) {
    sapply(1:dim(Xs)[1], function(i) {
      rpois(n=1, lambda=exp(scale_off[k]*(a*mean(covar.offset*causalBatch:::sigmoid(b*Xs[i,])) + eff_sz*batches[i])))
    })
  })
  batches <- factor(batches, levels=c("0", "1"))
  
  xtrue_tmp <- seq(from=-1, to=1, length.out=nbreaks)
  xtrue <- c(xtrue_tmp, xtrue_tmp)
  batch_true <- c(rep(0, nbreaks), rep(1, nbreaks))
  ytrue <- sapply(1:d, function(k) {
    sapply(xtrue, function(x) {
      exp(scale_off[k]*(a*causalBatch:::sigmoid(b*x)))
    })*exp(eff_sz*batch_true)
  })
  if (d == 1) {
    ys <- matrix(ys, ncol=1)
    ytrue <- matrix(ytrue, ncol=1)
  }
  
  return(list(Ys=ys, Ts=factor(batches), Xs=Xs, Effect=scale_off*eff_sz,
              x.bounds=c(-1, 1), Ytrue=ytrue, Xtrue=xtrue, 
              Ttrue=factor(batch_true,levels=levels(batches)),
              Overlap=causalBatch:::cb.sims.get_beta_overlap(alpha, beta, beta, alpha)))
}

sim_linear <- function(n=100, pi=.5, eff_sz=1, alpha=2, d=2, unbalancedness=1, 
                       a=2, b=1, err=1/2, unobs.cor = NULL, ncovar=1, unobs.alpha=1,
                       unobs.beta=3, nbreaks=200) {
  batches <- rbinom(n=n, size=1, prob=pi)
  
  beta <- alpha*unbalancedness
  if (!is.null(unobs.cor)) {
    corr.mtx <- make_corr_mtx(ncovar,unobs.cor)
    alphas <- c(rep(alpha, ncovar), unobs.alpha); betas <- c(rep(beta, ncovar), unobs.beta)
  } else {
    corr.mtx = diag(ncovar)
    alphas <- alpha; betas <- beta
  }
  Xs <- generate_correlated_betas_pergroup(batches, alphas, betas, betas, alphas, corr.mtx=corr.mtx)
  scale_off <- seq(1/d, 1/d^2, length.out=d)
  if (dim(Xs)[2] == 1) {
    covar.offset <- 1
  } else {
    covar.offset <- c(1, rep(1, dim(Xs)[2] - 1))
  }
  ys <- sapply(1:d, function(k) {
    sapply(1:dim(Xs)[1], function(i) {
      rpois(n=1, lambda=exp(scale_off[k]*(a*mean(covar.offset*(Xs[i,] + b)) + eff_sz*batches[i])))
    })
  })
  batches <- factor(batches, levels=c("0", "1"))
  
  xtrue_tmp <- seq(from=-1, to=1, length.out=nbreaks)
  xtrue <- c(xtrue_tmp, xtrue_tmp)
  batch_true <- c(rep(0, nbreaks), rep(1, nbreaks))
  ytrue <- sapply(1:d, function(k) {
    sapply(xtrue, function(x) {
      exp(scale_off[k]*(a*(x + b)))
    })*exp(eff_sz*batch_true)
  })
  if (d == 1) {
    ys <- matrix(ys, ncol=1)
    ytrue <- matrix(ytrue, ncol=1)
  }
  
  return(list(Ys=ys, Ts=factor(batches), Xs=Xs, Effect=scale_off*eff_sz,
              x.bounds=c(-1, 1), Ytrue=ytrue, Xtrue=xtrue, 
              Ttrue=factor(batch_true,levels=levels(batches)),
              Overlap=causalBatch:::cb.sims.get_beta_overlap(alpha, beta, beta, alpha)))
}

sim_impulse <- function(n=100, pi=.5, eff_sz=1, alpha=2, d=2, unbalancedness=1, 
                       a=.5, b=1/2, c=4, unobs.cor=NULL, ncovar=1, unobs.alpha=1,
                       unobs.beta=3, nbreaks=200) {
  batches <- rbinom(n=n, size=1, prob=pi)
  
  beta <- alpha*unbalancedness
  if (!is.null(unobs.cor)) {
    corr.mtx <- make_corr_mtx(ncovar,unobs.cor)
    alphas <- c(rep(alpha, ncovar), unobs.alpha); betas <- c(rep(beta, ncovar), unobs.beta)
  } else {
    corr.mtx = diag(ncovar)
    alphas <- alpha; betas <- beta
  }
  Xs <- generate_correlated_betas_pergroup(batches, alphas, betas, betas, alphas, corr.mtx=corr.mtx)
  scale_off <- seq(1/d, 1/d^2, length.out=d)
  if (dim(Xs)[2] == 1) {
    covar.offset <- 1
  } else {
    covar.offset <- c(1, rep(1, dim(Xs)[2] - 1))
  }
  ys <- sapply(1:d, function(k) {
    sapply(1:dim(Xs)[1], function(i) {
      rpois(n=1, lambda=exp(scale_off[k]*(mean(dnorm(covar.offset*Xs[i,], mean=a, sd=b)*c)) + eff_sz*batches[i]))
    })
  })
  batches <- factor(batches, levels=c("0", "1"))
  
  xtrue_tmp <- seq(from=-1, to=1, length.out=nbreaks)
  xtrue <- c(xtrue_tmp, xtrue_tmp)
  batch_true <- c(rep(0, nbreaks), rep(1, nbreaks))
  ytrue <- sapply(1:d, function(k) {
    sapply(xtrue, function(x) {
      exp(scale_off[k]*(dnorm(x, mean=a, sd=b)*c))
    })*exp(eff_sz*batch_true)
  })
  if (d == 1) {
    ys <- matrix(ys, ncol=1)
    ytrue <- matrix(ytrue, ncol=1)
  }
  
  return(list(Ys=ys, Ts=factor(batches), Xs=Xs, Effect=scale_off*eff_sz,
              x.bounds=c(-1, 1), Ytrue=ytrue, Xtrue=xtrue, 
              Ttrue=factor(batch_true,levels=levels(batches)),
              Overlap=causalBatch:::cb.sims.get_beta_overlap(alpha, beta, beta, alpha)))
}
sim_sigmoid <- function(n=100, pi=.5, eff_sz=1, alpha=2, d=2, unbalancedness=1, 
                        a=4, b=8, err=1/2, nbreaks=200) {
  batches <- rbinom(n=n, size=1, prob=pi)
  
  beta <- alpha*unbalancedness
  xs <- causalBatch:::cb.sims.covar_generator(batches, alpha, beta, beta, alpha)
  scale_off <- seq(1/d, 1/d^2, length.out=d)
  ys <- sapply(1:d, function(k) {
    sapply(1:length(xs), function(i) {
      rpois(n=1, lambda=exp(scale_off[k]*(a*causalBatch:::sigmoid(b*xs[i]) + eff_sz*batches[i])))
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
  
  return(list(Ys=ys, Ts=factor(batches), Xs=xs,
              x.bounds=c(-1, 1), Ytrue=ytrue, Xtrue=xtrue, 
              Ttrue=factor(batch_true,levels=levels(batches)),
              Overlap=causalBatch:::cb.sims.get_beta_overlap(alpha, beta, beta, alpha)))
}

sim_linear <- function(n=100, pi=.5, eff_sz=1, alpha=2, d=2, unbalancedness=1, 
                       a=2, b=1, err=1/2, nbreaks=200) {
  batches <- rbinom(n=n, size=1, prob=pi)
  
  beta <- alpha*unbalancedness
  xs <- causalBatch:::cb.sims.covar_generator(batches, alpha, beta, beta, alpha)
  scale_off <- seq(1/d, 1/d^2, length.out=d)
  ys <- sapply(1:d, function(k) {
    sapply(1:length(xs), function(i) {
      rpois(n=1, lambda=exp(scale_off[k]*(a*(xs[i] + b) + eff_sz*batches[i])))
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
  
  return(list(Ys=ys, Ts=factor(batches), Xs=xs,
              x.bounds=c(-1, 1), Ytrue=ytrue, Xtrue=xtrue, 
              Ttrue=factor(batch_true,levels=levels(batches)),
              Overlap=causalBatch:::cb.sims.get_beta_overlap(alpha, beta, beta, alpha)))
}


sim_impulse <- function(n=100, pi=.5, eff_sz=1, alpha=2, d=2, unbalancedness=1, 
                       a=.5, b=1/2, c=4, nbreaks=200) {
  batches <- rbinom(n=n, size=1, prob=pi)
  
  beta <- alpha*unbalancedness
  xs <- causalBatch:::cb.sims.covar_generator(batches, alpha, beta, beta, alpha)
  scale_off <- seq(1/d, 1/d^2, length.out=d)
  ys <- sapply(1:d, function(k) {
    sapply(1:length(xs), function(i) {
      rpois(n=1, lambda=exp(scale_off[k]*(dnorm(xs[i], mean=a, sd=b)*c + eff_sz*batches[i])))
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
  
  return(list(Ys=ys, Ts=factor(batches), Xs=xs,
              x.bounds=c(-1, 1), Ytrue=ytrue, Xtrue=xtrue, 
              Ttrue=factor(batch_true,levels=levels(batches)),
              Overlap=causalBatch:::cb.sims.get_beta_overlap(alpha, beta, beta, alpha)))
}
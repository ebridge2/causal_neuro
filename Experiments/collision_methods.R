require(MatchIt)
require(WeightIt)
require(SuperLearner)
require(drtmle)
require(causalBatch)

conditional.collision.lo <- function(Ys, Ts, Xs, ...) {
  sim.dat <- data.frame(Ys=Ys, Ts=Ts, Xs=Xs)
  mod <- glm("Ys ~ Ts + Xs", family="binomial", data=sim.dat)
  return(list(est=mod$coefficients["Ts"], pvalue=summary(mod)$coefficients["Ts", "Pr(>|z|)"]))
}

filtering.collision.lo <- function(Ys, Ts, Xs, usable, ...) {
  sim.dat <- data.frame(Ys=Ys, Ts=Ts, Xs=Xs)
  mod <- glm("Ys ~ Ts + Xs", family="binomial", data=sim.dat[which(usable),])
  return(list(est=mod$coefficients["Ts"], pvalue=summary(mod)$coefficients["Ts", "Pr(>|z|)"]))
}

sigmoid.xfm <- function(Xs, a=1, b=1) {
  a*causalBatch:::sigmoid(b*Xs)
}

linear.xfm <- function(Xs, a=2, b=1) {
  return(a*Xs + b)
}

impulse.xfm <- function(Xs, a=1/2, b=1/2, c=2) {
  return(dnorm(Xs, mean=a, sd=b)*c)
}

oracle.collision.lo <- function(Ys, Ts, Xs, usable, Ttrue=NULL, xfm=linear.xfm, ...) {
  sim.dat <- data.frame(Ys=Ys, Ts=Ttrue, Xs=do.call(xfm, list(Xs)))
  mod <- glm("Ys ~ Ts + Xs", family="binomial", data=sim.dat)
  return(list(est=mod$coefficients["Ts"], pvalue=summary(mod)$coefficients["Ts", "Pr(>|z|)"]))
}

causal.stratified.lo <- function(Ys, Ts, Xs, usable, ...) {
  aug.data <- data.frame(Ts=Ts, Xs=Xs, Ys=Ys, Usable=usable)
  W <- weightit(Ts ~ Xs + Usable, data = aug.data,
                moments = 2, int = TRUE,
                method = "ebal")
  
  fit <- glm_weightit(Ys ~ Ts + Xs + Usable, family=binomial(), 
                      data=aug.data, weightit=W)
}

causal.motion.regress.matching <- function(Ys, Ts, Xs, usable) {
  # match unusable points (cases) to usable points (controls)
  usability.df <- data.frame(Unusable=as.numeric(!usable), Xs=Xs, Ys=Ys, Ts=Ts)
  # sample with replacement
  matches <- matchit(Unusable ~ Ts + Xs + Ys, exact=~Ys, data=usability.df, method="optimal", distance="mahalanobis",
                     replace=FALSE, estimand = "ATT", ratio=1)
  
  idx.retain <- c(as.numeric(rownames(matches$match.matrix)), as.numeric(matches$match.matrix))
  
  matched.df <- usability.df[idx.retain,]
  match.fit.mod <- glm("Ts ~ Xs + Ys + Unusable", family=gaussian(), data=matched.df)
  
  Ts.tau <- Ts
  Ts.tau[which(!usable)] <- Ts[which(!usable)] - match.fit.mod$coefficients["Unusable"]
  return(Ts.tau)
}

causal.motion.regress.ipw <- function(Ys, Ts, Xs, usable) {
  # match unusable points (cases) to usable points (controls)
  usability.df <- data.frame(Unusable=as.numeric(!usable), Xs=Xs, Ys=Ys, Ts=Ts)
  # weight to estimate ATT
  W <- weightit(Unusable ~ Xs + Ys, data=usability.df, estimand = "ATT")
  
  ipw.fit.mod <- glm_weightit("Ts ~ Xs + Ys + Unusable", family=gaussian(), data=usability.df,
                              weightit=W)
  Ts.tau <- Ts
  Ts.tau[which(!usable)] <- Ts[which(!usable)] - ipw.fit.mod$coefficients["Unusable"]
  return(Ts.tau)
}

causal.ipw.lo <- function(Ys, Ts, Xs, usable, ...) {
  Ts.tau <- causal.motion.regress.ipw(Ys, Ts, Xs, usable)
  
  # entropy balancing weights from Vegetabile (2021)
  aug.data <- data.frame(Ts.tau=Ts.tau, Xs=Xs, Ys=Ys)
  W <- weightit(Ts.tau ~ Xs, data = aug.data,
                moments = 2, int = TRUE,
                method = "ebal")
  
  fit <- glm_weightit(Ys ~ Ts.tau + Xs, family=binomial(), data=aug.data, weightit=W)
  
  return(list(est=fit$coefficients["Ts.tau"], pvalue=summary(fit)$coefficients["Ts.tau", "Pr(>|z|)"]))
}

causal.matching.lo <- function(Ys, Ts, Xs, usable, ...) {
  Ts.tau <- causal.motion.regress.matching(Ys, Ts, Xs, usable)
  
  # entropy balancing weights from Vegetabile (2021)
  aug.data <- data.frame(Ts.tau=Ts.tau, Xs=Xs, Ys=Ys)
  W <- weightit(Ts.tau ~ Xs, data = aug.data,
                moments = 2, int = TRUE,
                method = "ebal")
  
  fit <- glm_weightit(Ys ~ Ts.tau + Xs, family=binomial(), data=aug.data, weightit=W)
  
  return(list(est=fit$coefficients["Ts.tau"], pvalue=summary(fit)$coefficients["Ts.tau", "Pr(>|z|)"]))
}

conditional.collision.md <- function(Ys, Ts, Xs, ...) {
  sim.dat <- data.frame(Ys=Ys, Ts=Ts, Xs=Xs)
  mod <- glm("Ts ~ Ys + Xs", data=sim.dat)
  return(list(est=mod$coefficients["Ys"], pvalue=summary(mod)$coefficients["Ys", "Pr(>|t|)"]))
}

filtering.collision.md <- function(Ys, Ts, Xs, usable, ...) {
  sim.dat <- data.frame(Ys=Ys, Ts=Ts, Xs=Xs)
  mod <- glm("Ts ~ Ys + Xs", data=sim.dat[which(usable),])
  return(list(est=mod$coefficients["Ys"], pvalue=summary(mod)$coefficients["Ys", "Pr(>|t|)"]))
}

oracle.collision.md <- function(Ys, Ts, Xs, usable, Ttrue=NULL, xfm=c(), ...) {
  sim.dat <- data.frame(Ys=Ys, Ts=Ttrue, Xs=do.call(xfm, list(Xs)))
  mod <- glm("Ts ~ Ys + Xs", data=sim.dat)
  return(list(est=mod$coefficients["Ys"], pvalue=summary(mod)$coefficients["Ys", "Pr(>|t|)"]))
}

oracle.naive.collision.md <- function(Ys, Ts, Xs, usable, Ttrue=NULL, xfm=c(), ...) {
  sim.dat <- data.frame(Ys=Ys, Ts=Ttrue, Xs=do.call(xfm, list(Xs)))
  mod <- glm("Ts ~ Ys", data=sim.dat)
  return(list(est=mod$coefficients["Ys"], pvalue=summary(mod)$coefficients["Ys", "Pr(>|t|)"]))
}

nebel.collision.md <- function(Ys, Ts, Xs, usable, ...) {
  Xmat <- data.frame(Ys, Xs)
  my.SL.libs.gn = c("SL.earth","SL.glmnet","SL.gam","SL.glm","SL.ranger","SL.step","SL.step.interaction","SL.xgboost","SL.mean")
  (gn.est = mcSuperLearner(Y = as.numeric(usable), X = Xmat, family=binomial(link='logit'),SL.library = my.SL.libs.gn, cvControl = list(V = 10), method='method.CC_nloglik'))
  
  my.SL.libs.Qbar = c("SL.earth","SL.glmnet","SL.gam","SL.glm","SL.ranger","SL.ridge","SL.step","SL.step.interaction","SL.svm","SL.xgboost","SL.mean")
  (Qbar.est =   mcSuperLearner(Y = Ts[usable], X=Xmat[usable, ], family=gaussian(), 
                               SL.library = my.SL.libs.Qbar, cvControl = list(V = 10), method = drtmle:::tmp_method.CC_LS))
  
  Qbar.predict = predict(Qbar.est, newdata = Xmat)[[1]]
  
  mean_fconn_asd.SL <- drtmle(Y = Ts[Ys==1],
                              A = as.numeric(usable[Ys==1]), # apologies for the notation conflict -- use Delta here
                              W = NULL, # does not do anything with user-input propensities and outcomes
                              a_0 = 1, # set this to one to correspond to counterfactual that all Delta=1
                              Qn = list(Qbar.predict[Ys==1]), # pass in fitted outcome values
                              gn = list(gn.est$SL.predict[Ys==1]), # pass in fitted propensities
                              SL_Qr = "SL.npreg", # uses non-parametric regression in the drtmle step
                              SL_gr = "SL.npreg",
                              maxIter = 1
  )
  
  mean_fconn_td.SL <- drtmle(Y = Ts[Ys==0],
                             A = as.numeric(usable[Ys==0]),
                             W = NULL,
                             a_0 = 1,
                             Qn = list(Qbar.predict[Ys==0]),
                             gn = list(gn.est$SL.predict[Ys==0]),
                             SL_Qr = "SL.npreg",
                             SL_gr = "SL.npreg",
                             maxIter = 1
  )
  
  difference = mean_fconn_asd.SL[[1]]$est - mean_fconn_td.SL[[1]]$est
  
  
  se.diff = sqrt(mean_fconn_asd.SL[[1]]$cov + mean_fconn_td.SL[[1]]$cov)
  z.diff = as.numeric(difference/se.diff)
  return(list(est=difference, pvalue=2*(1-pnorm(abs(z.diff)))))
}
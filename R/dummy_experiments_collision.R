sim <- collider.sim_linear(n=500)

Ys <- sim$Ys; Ts <- sim$Ts; Xs <- sim$Xs; usable <- sim$usable; Ttrue=sim$Ttrue
risk <- sim$Risk.group

Ts.tau <- causal.motion.regress(Ys, Ts, Xs, usable)

sum(abs(Ts.tau - Ttrue))
sum(abs(Ts - Ttrue))

df <- data.frame(Variable=rep(c("Ttrue", "Ttrue", "Ts", "Ts", "Ttau", "Ttau"), 2), Mean=NA, usability=rep(c(FALSE, TRUE, FALSE, TRUE, FALSE, TRUE),2), Ys=c(c(0,0,0,0,0,0), c(1,1,1,1,1,1)))
df[1,2] = mean(Ttrue[usable == 0 & Ys == 0])
df[2,2] = mean(Ttrue[usable == 1 & Ys == 0])

df[3,2] = mean(Ts[usable == 0 & Ys == 0])
df[4,2] = mean(Ts[usable == 1 & Ys == 0])

df[5,2] = mean(Ts.tau[usable == 0 & Ys == 0])
df[6,2] = mean(Ts.tau[usable == 1 & Ys == 0])

df[7,2] = mean(Ttrue[usable == 0 & Ys == 1])
df[8,2] = mean(Ttrue[usable == 1 & Ys == 1])

df[9,2] = mean(Ts[usable == 0 & Ys == 1])
df[10,2] = mean(Ts[usable == 1 & Ys == 1])

df[11,2] = mean(Ts.tau[usable == 0 & Ys == 1])
df[12,2] = mean(Ts.tau[usable == 1 & Ys == 1])
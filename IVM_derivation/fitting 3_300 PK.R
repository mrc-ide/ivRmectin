rm(list = ls())
library(odin)

# fit PK model

#### basic PK model ####

pk <- odin::odin({

  ka <- 0.474
  ke <- (10.9*(WT/60))/(161.7*(WT/60)^0.75)
  k23 <- 0.123
  k32 <- 0.033

  init_X1 <- user()
  init_X2 <- user()
  init_X3 <- user()
  WT <- user()

  initial(X1) <- init_X1
  initial(X2) <- init_X2
  initial(X3) <- init_X3

  deriv(X1) <- -ka*X1
  deriv(X2) <- ka*X1 - (ke + k23)*X2 + k32*X3
  deriv(X3) <- k23*X2 - k32*X3


})


#### function for running pk model with three doses ####


dose_fit = function(dose, WT=60){

  pars = list(init_X1 = dose, init_X2 = 0, init_X3 = 0, WT=WT)

  mod <- pk(user = pars)
  t <- seq(0, 24*30, by=0.1)
  y <- data.frame(mod$run(t))

  inits_d1 =as.numeric(y[round(y$t,1)==23.9,])

  pars2 = list(init_X1 = inits_d1[2]+dose, init_X2 = inits_d1[3], init_X3 = inits_d1[4], WT=WT)

  mod2 <- pk(user = pars2)
  y2 <- data.frame(mod2$run(t))


  inits_d2 =as.numeric(y2[round(y2$t,1)==23.9,])

  pars3 = list(init_X1 = inits_d2[2]+dose, init_X2 = inits_d2[3], init_X3 = inits_d2[4], WT=WT)

  mod3 <- pk(user = pars3)
  y3 <- data.frame(mod3$run(t))

  yy = c(y$X3[which(y$t==0):which(round(y$t,1)==23.9)],
         y2$X3[which(y2$t==0):which(round(y2$t,1)==23.9)],
         y3$X3)

  yy2 = c(y$X2[which(y$t==0):which(round(y$t,1)==23.9)],
          y2$X2[which(y2$t==0):which(round(y2$t,1)==23.9)],
          y3$X2)

  tind = min(which(y$t>=48))
  tmax = y$t[which.max(yy2[tind:length(y$t)])]
  cmax = max(yy2[1:length(y$t)])
  return(list(c(tmax=tmax, cmax=cmax), yy2))
}


# the cmax from IVERMAL with 3x300 dose is 69.4
# table 1 in https://ascpt.onlinelibrary.wiley.com/doi/epdf/10.1002/cpt.1219


mm = matrix(NA, ncol=2, nrow=100)

dose_i = seq(95, 100, l=100)

for(i in 1:100){
  mm[i,] = dose_fit(dose_i[i], WT=60)[[1]]
}

# find the dose that gives the correct Cmax value
dose_new = dose_i[which.min(abs(mm[,2] - 69.4))]

tt <- seq(0, 24*30, by=0.1)
yy2_new = dose_fit(dose_new)[[2]]

plot(tt/24, yy2_new[1:length(tt)], type="l", las=1, ylab="IVM concentration (ng/ml)", xlab="Time (days)")


write.table(data.frame(time = tt, pk3_300 = yy2_new[1:length(tt)]), file = "ivm_derivation/pk_300_hs.txt", quote = F)

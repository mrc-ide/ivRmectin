# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # #                                                                         # # #
# # # # # #       DERIVING HAZARD RATIOS FROM SLATER ET AL. 2020 LID IVERMECTIN     # # #
# # # # # #                                                                         # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rm(list = ls())

#### 1) raw hazard ratio data from Smit et al. ####
# the hazard ratios are takes from Smit et al. Supplementary table S3
# https://www.sciencedirect.com/science/article/pii/S1473309918301634
# https://ars.els-cdn.com/content/image/1-s2.0-S1473309918301634-mmc1.pdf

#times (days)
vt = c(2+4/24, 7, 10, 14, 21, 28)
# hazard ratios - column IVM-3x300 vs Placebo
m300 = c(9.76, 4.23, 2.56, 2.02,1.54 ,1.38)

plot(vt, m300, xlab = "Day", ylab = "Hazard ratio", pch=21, bg = "firebrick2", cex=1.5,
     ylim = c(0,12), xlim = c(0, 28), las = 1, main = "Hazard ratios from Smit et al. ")

# a line was fitting to these based on a coxph model detailed in Table S19 of
# same supplementary material linked above
haz_ivermal = read.table("ivm_derivation/hazards_IVERMAL_FINAL.txt", header=T)

lines(haz_ivermal$Day, haz_ivermal$m300)


#### 2) PK model fit for 3 x 300 ####

# full details on PK derivation are in 'ivm_derivation/fitting_3_300_PK.R'
pks = read.table("ivm_derivation/pk_300_hs.txt", header=T)

plot(pks$time/24, pks$pk3_300, type="l", xlab = "Time (hours)", lwd=3, col="firebrick3",
     las = 1, ylab = "Ivermectin concentration (ng/ml)")


#### 3) mean concentration for each day post-ingestion ####

# take the average concentration across each 24 hour period

time = seq(0,30*24, l=200)

mean_conc = function(time, pk){
  mc = rep(NA, length(1:max(time))/24)
  for(i in 1:length(mc)){
    ind = which(time/24 >= (i-1) & time/24 < i)
    mc[i] = mean(pk[ind])
  }
  return(mc)
}

mc3_300 = mean_conc(time = pks$time, pk = pks$pk3_300)

#### 4)

# fit a model between the ivermal hazards and the pk concentrations
# fit using data from day 3 onwards
yfit = loess(haz_ivermal$m300[2:27] ~ log10(mc3_300[3:28]))

# predict over the daily concentration values
# a benefit of this approach is that you get meaning hazard ratios for day 1 and 2
# it also means you can apply to other doses by deriving new PK profiles for other doses
y3_300 = predict(yfit, log10(mc3_300))

write.table(data.frame(day=1:30, d300 = y3_300), "IVM_derivation/hazards_300.txt", row.names=FALSE)



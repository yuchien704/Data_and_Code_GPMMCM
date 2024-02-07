load(paste(SPATH, "/Data/CM_1.RData", sep = ""))
load(paste(SPATH, "/Data/CM_2.RData", sep = ""))
load(paste(SPATH, "/Data/CM_3.RData", sep = ""))
load(paste(SPATH, "/Data/CM_4.RData", sep = ""))
load(paste(SPATH, "/Data/CM_5.RData", sep = ""))

total_sum_1 = total_sum_2 = total_sum_3 = total_sum_4 = total_sum_5 = matrix(0, 5,14)  # 初始总和为0
for (i in 1:20) {
  total_sum_1 <- total_sum_1 + (SIM_VVE_20_0609$BIC)[[i]]
  total_sum_2 <- total_sum_2 + (SIM_VVE_20_2$BIC)[[i]]
  total_sum_3 <- total_sum_3 + (SIM_VVE_20_3$BIC)[[i]]
  total_sum_4 <- total_sum_4 + (SIM_VVE_20_4$BIC)[[i]]
  total_sum_5 <- total_sum_5 + (SIM_VVE_CM_1$BIC)[[i]]
}


BIC.mean = (total_sum_1+total_sum_2+total_sum_3+total_sum_4+total_sum_5)/100

bestModel_BIC_all = rbind(SIM_VVE_20_0609$bestModel,SIM_VVE_20_2$bestModel,SIM_VVE_20_3$bestModel,SIM_VVE_20_4$bestModel,SIM_VVE_CM_1$bestModel)
Tab_bestModel_BIC = table(bestModel_BIC_all)
freq_best = matrix(0, 5,14, dimnames = list(c(1:k),c("EII","VII","EEI","VEI","EVI","VVI","EEE","VEE","EVE","VVE","EEV","VEV","EVV","VVV")))
freq_best[3,12]=freq_best[4, c(10,11)]=1; freq_best[3,8]=20; freq_best[3,10]=77

correct.rate_BIC = sum(SIM_VVE_20_0609$correct.rate,SIM_VVE_20_2$correct.rate,SIM_VVE_20_3$correct.rate,
                        SIM_VVE_20_4$correct.rate,SIM_VVE_CM_1$correct.rate)/5

Table7 = rbind(BIC.mean[1,], freq_best[1,], BIC.mean[2,], freq_best[2,], BIC.mean[3,], freq_best[3,], 
               BIC.mean[4,], freq_best[4,], BIC.mean[5,], freq_best[5,])
rownames(Table7) = c("BIC_g=1", "freq_=g=1", "BIC_g=2", "freq_=g=2", "BIC_g=3", "freq_=g=3", "BIC_g=4", "freq_=g=4","BIC_g=5", "freq_=g=5")

write.csv(Table7, paste(SPATH,'/Results/Table7.csv',sep=""), row.names = TRUE)

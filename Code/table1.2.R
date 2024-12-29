library(mvtnorm)
library(tmvtnorm)
library(mclust)
library(MomTrunc)
load(paste(SPATH, "/Data/hawks_new.RData", sep = ""))

Y = Y.na
g = list()
num.na = mean = sd = matrix(NA, 5,3)
for(i in 1:5){
  for(j in 1:3){
    na.posi = which(is.na(Y[true.clus==j,i]))
    num.na[i,j] = length(which(is.na(Y[true.clus==j,i])))
    g[[j]] = which(true.clus==j)
    mean[,j] = colMeans(Y[g[[j]],],na.rm = T)
    sd[,j] = apply(na.omit(Y[g[[j]],]), 2, sd)
    Y[g[[j]][na.posi],i] = mean[i,j] 
  }
}
Table1 = cbind(round(mean,3)[,1], round(sd,3)[,1], num.na[,1], 
              round(mean,3)[,2], round(sd,3)[,2], num.na[,2],
              round(mean,3)[,3], round(sd,3)[,3], num.na[,3])
colnames(Table1) = c("CH_mean", "CH_sd", "CH_num.na","RT_mean", "RT_sd", "RT_num.na", "SS_mean", "SS_sd", "SS_num.na")
write.csv(Table1, paste(SPATH,'/Results/Table1.csv',sep=""), row.names = TRUE)


Table2=rbind(fit_VVV$model.inf[c(8,9,1)], fit_EEE$model.inf[c(8,9,1)], fit_EEV$model.inf[c(8,9,1)], fit_EVE$model.inf[c(8,9,1)],
             fit_EVV$model.inf[c(8,9,1)], fit_VEE$model.inf[c(8,9,1)], fit_VEV$model.inf[c(8,9,1)], fit_VVE$model.inf[c(8,9,1)],
             fit_VEI$model.inf[c(8,9,1)], fit_EEI$model.inf[c(8,9,1)], fit_VVI$model.inf[c(8,9,1)], fit_EVI$model.inf[c(8,9,1)],
             fit_VII$model.inf[c(8,9,1)], fit_EII$model.inf[c(8,9,1)],c(fit.mnc.ARI, fit.mnc.CCR, fit.mnc$res$time[[1]]),
             c(fit.msnc.ARI, fit.msnc.CCR, fit.msnc$res$time[[1]]), c(ARI1, CCR1,time_Amelia), c(ARI2, CCR2,time_mice),
             c(ARI3, CCR3,time_mi), c(ARI4, CCR4,time_kpod))
rownames(Table2)=c('VVV','EEE','EEV','EVE','EVV','VEE','VEV','VVE','VEI','EEI','VVI','EVI','VII','EII',
                   'FM-MNC','FM-MSNC','Amelia','mice','mi','kpod')
# t(round(Table2,4))
write.csv(t(round(Table2,4)), paste(SPATH,'/Results/Table2.csv',sep=""), row.names = TRUE)

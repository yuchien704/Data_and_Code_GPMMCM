library(mvtnorm)
library(tmvtnorm)
library(mclust)
library(MomTrunc)
load(paste(SPATH, "/Data/hawks.RData", sep = ""))

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
Table2 = cbind(round(mean,3)[,1], round(sd,3)[,1], num.na[,1], 
              round(mean,3)[,2], round(sd,3)[,2], num.na[,2],
              round(mean,3)[,3], round(sd,3)[,3], num.na[,3])
colnames(Table2) = c("CH_mean", "CH_sd", "CH_num.na","RT_mean", "RT_sd", "RT_num.na", "SS_mean", "SS_sd", "SS_num.na")
write.csv(Table2, paste(SPATH,'/Results/Table2.csv',sep=""), row.names = TRUE)

a3 = round(a3)
n = nrow(Ycm)
pclus = matrix(NA, n, 14)
pclus[,1] = fit_VVV3$post.clus
pclus[,2] = fit_EEE3$post.clus
pclus[,3] = fit_EEV3$post.clus
pclus[,4] = fit_EVE3$post.clus
pclus[,5] = fit_EVV3$post.clus
pclus[,6] = fit_VEE3$post.clus
pclus[,7] = fit_VEV3$post.clus
pclus[,8] = fit_VVE3$post.clus
pclus[,9] = fit_VEI3$post.clus
pclus[,10] = fit_EEI3$post.clus
pclus[,11] = fit_VVI3$post.clus
pclus[,12] = fit_EVI3$post.clus
pclus[,13] = fit_VII3$post.clus
pclus[,14] = fit_EII3$post.clus

ari_hawks = CCR_hawks = c()
for(i in 1:14){
  ari_hawks[i] = adjustedRandIndex(true.clus, pclus[,i])
  CCR_hawks[i] = 1-classError(true.clus, pclus[,i])$errorRate
}
names(ari_hawks) = names(CCR_hawks) = c("VVV", "EEE", "EEV","EVE","EVV","VEE","VEV","VVE","VEI","EEI","VVI","EVI","VII","EII")

Table3 = rbind(round(ari_hawks,4),round(CCR_hawks,4))
rownames(Table3) = c("ARI", "CCR")
write.csv(Table3, paste(SPATH,'/Results/Table3.csv',sep=""), row.names = TRUE)

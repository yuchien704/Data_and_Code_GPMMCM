library(mvtnorm)
library(tmvtnorm)
library(mclust)
library(MomTrunc)

source(paste(SPATH, '/Function/gen_cen_na.r', sep=''))
source(paste(SPATH, '/Function/F-G.R', sep=''))
source(paste(SPATH, '/Function/Both_gmixcm/EEE_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/Both_gmixcm/EEI_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/Both_gmixcm/EEV_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/Both_gmixcm/EII_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/Both_gmixcm/EVI_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/Both_gmixcm/EVV_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/Both_gmixcm/VEE_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/Both_gmixcm/VEI_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/Both_gmixcm/VEV_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/Both_gmixcm/VII_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/Both_gmixcm/VVI_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/Both_gmixcm/VVV_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/Both_gmixcm/EVE_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/Both_gmixcm/VVE_GMIXCMB.R', sep=''))

Data = read.csv(paste(SPATH, '/Data/source/Hawks.csv', sep=''))
Y.na = Data[,2:6]
true.clus = Data[,1]
set.seed(0307)
NC=gener.cen.na(Data = Y.na,cen.type = c(0,0,0,1,1),cen.rate=c(0,0,0,0.1,0.1),na.rate = 0)
Ycm = NC$Data
cen = NC$cen
cen[is.na(Y.na)]=NA

cen.value_4 = Ycm[which(cen[,4]==1), 4][1]
cen.value_5 = Ycm[which(cen[,5]==1), 5][1]

cen.num = colSums(cen, na.rm = T)
na.num = colSums(is.na(cen))

g=3;clus3 = true.clus
tol = 1e-5;neq.S = T;max.iter=1000;tol.FG=1e-3;max.iter.FG=10; per=1
fit_VVV3 = VVV.GMIX.MVNCM.B.ECM(Ycm, cen, g, clus3, tol, neq.S = T, max.iter, per)
fit_EEE3 = EEE.GMIX.MVNCM.B.ECM(Ycm, cen, g, clus3, tol, neq.S = F, max.iter, per)
fit_EEV3 = EEV.GMIX.MVNCM.B.ECM(Ycm, cen, g, clus3, tol, neq.S = T, max.iter, per)
fit_EVE3 = EVE.GMIX.MVNCM.B.ECM(Ycm, cen, g, clus3, tol, neq.S = T, max.iter, tol.FG, max.iter.FG, per)
fit_EVV3 = EVV.GMIX.MVNCM.B.ECM(Ycm, cen, g, clus3, tol, neq.S = T, max.iter, per)
fit_VEE3 = VEE.GMIX.MVNCM.B.ECM(Ycm, cen, g, clus3, tol, neq.S = T, max.iter, per)
fit_VEV3 = VEV.GMIX.MVNCM.B.ECM(Ycm, cen, g, clus3, tol, neq.S = T, max.iter, per)
fit_VVE3 = VVE.GMIX.MVNCM.B.ECM(Ycm, cen, g, clus3, tol, neq.S = T, max.iter, tol.FG, max.iter.FG, per)

fit_VEI3 = VEI.GMIX.MVNCM.B.ECM(Ycm, cen, g, clus3, tol, neq.S = T, max.iter, per)
fit_EEI3 = EEI.GMIX.MVNCM.B.ECM(Ycm, cen, g, clus3, tol, neq.S = F, max.iter, per)
fit_VVI3 = VVI.GMIX.MVNCM.B.ECM(Ycm, cen, g, clus3, tol, neq.S = T, max.iter, per)
fit_EVI3 = EVI.GMIX.MVNCM.B.ECM(Ycm, cen, g, clus3, tol, neq.S = T, max.iter, per)

fit_VII3 = VII.GMIX.MVNCM.B.ECM(Ycm, cen, g, clus3, tol, neq.S = T, max.iter, per)
fit_EII3 = EII.GMIX.MVNCM.B.ECM(Ycm, cen, g, clus3, tol, neq.S = F, max.iter, per)
a3=round(c(fit_VVV3$model.inf[[5]],fit_EEE3$model.inf[[5]],fit_EEV3$model.inf[[5]],fit_EVE3$model.inf[[5]],
      fit_EVV3$model.inf[[5]],fit_VEE3$model.inf[[5]],fit_VEV3$model.inf[[5]],fit_VVE3$model.inf[[5]],
      fit_VEI3$model.inf[[5]],fit_EEI3$model.inf[[5]],fit_VVI3$model.inf[[5]],fit_EVI3$model.inf[[5]],
      fit_VII3$model.inf[[5]],fit_EII3$model.inf[[5]]))
names(a3) = paste(c("VVV", "EEE", "EEV","EVE","EVV","VEE","VEV","VVE","VEI","EEI","VVI","EVI","VII","EII"),3,sep = ",")
Best_hawks = sort(a3)[1:3]

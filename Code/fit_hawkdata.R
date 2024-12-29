library(mvtnorm)
library(tmvtnorm)
library(mclust)
library(MomTrunc)
library(Stat2Data)
library(VIM)
library(CensMFM)
library(kpodclustr)
library(Amelia)
library(mice)
library(mi)
source(paste(SPATH, '/Function/gen_cen_na_new.r', sep=''))
source(paste(SPATH, '/Function/F-G.R', sep=''))
source(paste(SPATH, '/Function/fn/EEE_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/fn/EEI_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/fn/EEV_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/fn/EII_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/fn/EVI_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/fn/EVV_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/fn/VEE_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/fn/VEI_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/fn/VEV_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/fn/VII_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/fn/VVI_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/fn/VVV_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/fn/EVE_GMIXCMB.R', sep=''))
source(paste(SPATH, '/Function/fn/VVE_GMIXCMB.R', sep=''))

data("Hawks")
Y.na = as.matrix(Hawks[,10:14])
true.clus = as.numeric(Hawks[,7])
p=ncol(Y.na);n=nrow(Y.na)
NC=gener.cen.na(Data=Y.na, cen.type=c(0,0,0,1,1), cen.rate=c(0,0,0,0.1,0.1), na.rate=0)
Ycm = NC$Data
cen = NC$cen
cen[is.na(Y.na)]=NA


# fit the GPMMCM model
set.seed(0307)
g=3; Y.knn = kNN(Ycm, k = 5)[,1:p]
fit_VVV = VVV.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-5, neq.S = T, max.iter=500, per=20)
fit_EEE = EEE.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-5, neq.S = T, max.iter=500, per=20)
fit_EEV = EEV.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-5, neq.S = T, max.iter=500, per=20)
fit_EVE = EVE.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-5, neq.S = T, tol.FG=1e-3, max.iter.FG=10, max.iter=500, per=20)
fit_EVV = EVV.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-5, neq.S = T, max.iter=500, per=20)
fit_VEE = VEE.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-5, neq.S = T, max.iter=500, per=20)
fit_VEV = VEV.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-5, neq.S = T, max.iter=500, per=20)
fit_VVE = VVE.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-5, neq.S = T, tol.FG=1e-3, max.iter.FG=10, max.iter=500, per=20)
fit_VEI = VEI.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-5, neq.S = T, max.iter=500, per=20)
fit_EEI = EEI.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-5, neq.S = T, max.iter=500, per=20)
fit_VVI = VVI.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-5, neq.S = T, max.iter=500, per=20)
fit_EVI = EVI.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-5, neq.S = T, max.iter=500, per=20)
fit_VII = VII.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-5, neq.S = T, max.iter=500, per=20)
fit_EII = EII.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-5, neq.S = T, max.iter=500, per=20)


# === other method ===#
#Amelia
start_time_1 <- Sys.time()
a1.out <- amelia(x = Ycm)
X1 = matrix(0,n,p)
A1 = a1.out$imputations
m1 = a1.out$m
for(j in 1:m1){
  X1 = X1 + A1[[j]]/m1
}
kmeansclusters1 <- kmeans(as.matrix(X1),g)$cluster
ARI1 = adjustedRandIndex(kmeansclusters1,true.clus)
CCR1 = 1-classError(kmeansclusters1,true.clus)$errorRate
end_time_1 <- Sys.time()
time_Amelia = end_time_1 - start_time_1

#mice
start_time_2 <- Sys.time()
b.out = mice(Ycm) 
X2 = b.out$data
B1.posi = B1.num = list()
for(j in 1:p) {
  B1.posi[[j]] = as.numeric(lapply(b.out$imp, row.names)[[j]])
  B1.num[[j]] = lapply(b.out$imp, rowMeans)[[j]]
  X2[B1.posi[[j]],j] = B1.num[[j]]
}
kmeansclusters2 <- kmeans(as.matrix(X2),g)$cluster
ARI2 = adjustedRandIndex(kmeansclusters2,true.clus)
CCR2 = 1-classError(kmeansclusters2,true.clus)$errorRate
end_time_2 <- Sys.time()
time_mice = end_time_2 - start_time_2

#mi
start_time_3 <- Sys.time()
c.out <- missing_data.frame(Ycm) # warnings about missingness patterns
c.out <- change(c.out, y = c("Wing","Weight","Culmen","Hallux","Tail"), what = "transformation", to = "identity")
imputations <- mi(c.out)
mi1 = (imputations@data$`chain:1`@X[,2:6])
mi2 = (imputations@data$`chain:2`@X[,2:6])
mi3 = (imputations@data$`chain:3`@X[,2:6])
mi4 = (imputations@data$`chain:4`@X[,2:6])
X3 = (mi1+mi2+mi3+mi4)/4
kmeansclusters3 <- kmeans(as.matrix(X3),g)$cluster
ARI3 = adjustedRandIndex(kmeansclusters3,true.clus)
CCR3 = 1-classError(kmeansclusters3,true.clus)$errorRate
end_time_3 <- Sys.time()
time_mi = end_time_3 - start_time_3

#kpod
start_time_4 <- Sys.time()
kpod_result <- kpod(as.matrix(Ycm),g)
kpodclusters <- kpod_result$cluster
ARI4 = adjustedRandIndex(kpodclusters,true.clus)
CCR4 = 1-classError(kpodclusters,true.clus)$errorRate
end_time_4 <- Sys.time()
time_kpod = end_time_4 - start_time_4


# fit the FM-MNC and FM-MSNC model
if(g==1){clus=rep(1,n)}  else{clus=kmeans(Y.knn, g)$cluster}
mu <- Sigma <- shape <- list()
pii=length(g)
n=nrow(Ycm)
for(i in 1:g)
  {
    pii[i] = nrow(Y.knn[clus == i,]) / n
    mu[[i]] = colMeans(Y.knn[clus == i,], na.rm = T)
    Sigma[[i]] = var(Y.knn, na.rm = T)
    shape[[i]] <- rep(0.1, g)
  }
#nu <- rep(0, g)


LI=cen
LI[is.na(cen)!=0]= -Inf
LI[cen==1]= -Inf
LI[cen==2 & is.na(cen)==0]= Ycm[cen==2 & is.na(cen)==0]
LS=cen
LS[is.na(cen)!=0]= Inf
LS[cen==2]= Inf
LS[cen==1 & is.na(cen)==0]=  Ycm[cen==1 & is.na(cen)==0]
cc=cen
cc[cen==2]=1
cc[is.na(cen)!=0]=1

start_time_mnc <- Sys.time()
fit.mnc = fit.FMMSNC(cc, LI, LS, Ycm, mu=mu,
                          Sigma = Sigma, shape=shape, pii = pii, g = g, get.init = F,
                          criteria = TRUE, family = "Normal", error = 0.00001,
                          iter.max = 350, uni.Gama = FALSE, cal.im = F)
fit.mnc.CCR=1-classError(true.clus, fit.mnc$res$group)$errorRate
fit.mnc.ARI=adjustedRandIndex(true.clus, fit.mnc$res$group)
end_time_mnc <- Sys.time()
# cat('FM-MNC running time=', end_time_mnc - start_time_mnc, 'CCR=', fit.mnc.CCR, 'ARI=', fit.mnc.ARI, '\n')
start_time_msnc <- Sys.time()
fit.msnc = fit.FMMSNC(cc, LI, LS, Ycm, mu=mu,
                          Sigma = Sigma, shape=shape, pii = pii, g = g, get.init = T,
                          criteria = TRUE, family = "SN", error = 0.00001,
                          iter.max = 350, uni.Gama = FALSE, cal.im = F)
fit.msnc.CCR=1-classError(true.clus, fit.msnc$res$group)$errorRate
fit.msnc.ARI=adjustedRandIndex(true.clus, fit.msnc$res$group)
end_time_msnc <- Sys.time()
# cat('FM-MSNC running time=', end_time_msnc - start_time_msnc, 'CCR=', fit.msnc.CCR, 'ARI=', fit.msnc.ARI, '\n')

Table2=rbind(fit_VVV$model.inf[c(8,9,1)], fit_EEE$model.inf[c(8,9,1)], fit_EEV$model.inf[c(8,9,1)], fit_EVE$model.inf[c(8,9,1)],
          fit_EVV$model.inf[c(8,9,1)], fit_VEE$model.inf[c(8,9,1)], fit_VEV$model.inf[c(8,9,1)], fit_VVE$model.inf[c(8,9,1)],
          fit_VEI$model.inf[c(8,9,1)], fit_EEI$model.inf[c(8,9,1)], fit_VVI$model.inf[c(8,9,1)], fit_EVI$model.inf[c(8,9,1)],
          fit_VII$model.inf[c(8,9,1)], fit_EII$model.inf[c(8,9,1)],c(fit.mnc.ARI, fit.mnc.CCR, fit.mnc$res$time[[1]]),
          c(fit.msnc.ARI, fit.msnc.CCR, fit.msnc$res$time[[1]]), c(ARI1, CCR1,time_Amelia), c(ARI2, CCR2,time_mice),
          c(ARI3, CCR3,time_mi), c(ARI4, CCR4,time_kpod))
rownames(Table2)=c('VVV','EEE','EEV','EVE','EVV','VEE','VEV','VVE','VEI','EEI','VVI','EVI','VII','EII',
                   'FM-MNC','FM-MSNC','Amelia','mice','mi','kpod')
t(round(Table2,4))


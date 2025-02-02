library(mvtnorm)
library(tmvtnorm)
library(mclust)
library(MomTrunc)
library(VIM)
library(CensMFM)
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

data=read.csv(paste(SPATH, "/Data/source/DMSdata.csv", sep = ""))
true.clus=as.numeric(factor(data[,1]))
Ycm=as.matrix(data[,-1])
p=ncol(Ycm); n=nrow(Ycm)
cen=matrix(0, n, p)
cen[is.na(Ycm[,1]),1]=NA
cen[,2]=(Ycm[,2]==2)*1
cen[,3]=(Ycm[,3]==4)*1
cen[,4]=(Ycm[,4]==11.3)*1
cen[,5]=(Ycm[,5]==0.22)*1

# fit the GPMMCM model
set.seed(7)
Y.knn = kNN(Ycm, k =5)[,1:p]
g=2
fit_VVV = VVV.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-6, neq.S = T, max.iter=500, per=20)
fit_EEE = EEE.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-6, neq.S = T, max.iter=500, per=20)
fit_VEE = VEE.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-6, neq.S = T, max.iter=500, per=20)
fit_EEV = EEV.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-6, neq.S = T, max.iter=500, per=20)
fit_VEV = VEV.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-6, neq.S = T, max.iter=500, per=20)
fit_EVV = EVV.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-6, neq.S = T, max.iter=500, per=20)
fit_EVE = EVE.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-6, neq.S = T, tol.FG=1e-3, max.iter.FG=10, max.iter=500, per=20)
fit_VVE = VVE.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-6, neq.S = T, tol.FG=1e-3, max.iter.FG=10, max.iter=500, per=20)
fit_EEI = EEI.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-6, neq.S = T, max.iter=500, per=20)
fit_VVI = VVI.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-6, neq.S = T, max.iter=500, per=20)
fit_VEI = VEI.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-6, neq.S = T, max.iter=500, per=20)
fit_EVI = EVI.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-6, neq.S = T, max.iter=500, per=20)
fit_VII = VII.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-6, neq.S = T, max.iter=500, per=20)
fit_EII = EII.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-6, neq.S = T, max.iter=500, per=20)
fit_EVE = EVE.GMIXCMB.ECM.new(Ycm, Y.knn, cen, g=g, true.clus=true.clus, tol = 1e-6, neq.S = T, tol.FG=1e-3, max.iter.FG=10, max.iter=500, per=20)


if(g==1){clus=rep(1,n)}  else{clus=kmeans(Y.knn, g)$cluster}
mu <- Sigma <- shape <- list()
pii=length(g)
n=nrow(Ycm)
for(i in 1:g)
{
  pii[i] = nrow(Y.knn[clus == i,]) / n
  mu[[i]] = colMeans(Y.knn[clus == i,], na.rm = T)
  Sigma[[i]] = var(Y.knn, na.rm = T)
  shape[[i]] <- rnorm(p)
}

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

# fit the FM-MNC and FM-MSNC model
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
                          Sigma = Sigma, shape=shape, pii = pii, g = g, get.init =F,
                          criteria = TRUE, family = "SN", error = 0.00001,
                          iter.max = 350, uni.Gama = FALSE, cal.im = F)
fit.msnc.CCR=1-classError(true.clus, fit.msnc$res$group)$errorRate
fit.msnc.ARI=adjustedRandIndex(true.clus, fit.msnc$res$group)
end_time_msnc <- Sys.time()
# cat('FM-MSNC running time=', end_time_msnc - start_time_msnc, 'CCR=', fit.msnc.CCR, 'ARI=', fit.msnc.ARI, '\n')


idx=c(10,3,5,8,9,1)
Table4=round(rbind(fit_VVV$model.inf[idx], fit_EEE$model.inf[idx], fit_EEV$model.inf[idx], fit_EVE$model.inf[idx],
                   fit_EVV$model.inf[idx], fit_VEE$model.inf[idx], fit_VEV$model.inf[idx], fit_VVE$model.inf[idx],
                   fit_VEI$model.inf[idx], fit_EEI$model.inf[idx], fit_VVI$model.inf[idx], fit_EVI$model.inf[idx],
                   fit_VII$model.inf[idx], fit_EII$model.inf[idx], 
                   c(fit_VVV$model.inf[10], fit.mnc$res$loglik, fit.mnc$res$bic, fit.mnc.ARI, fit.mnc.CCR, fit.mnc$res$time[[1]]),
                   c(fit_VVV$model.inf[10]+10, fit.msnc$res$loglik, fit.msnc$res$bic, fit.msnc.ARI, fit.msnc.CCR, fit.msnc$res$time[[1]])),3)
rownames(Table4) = c("VVV","EEE","EEV","EVE","EVV","VEE","VEV","VVE","VEI","EEI","VVI","EVI","VII","EII",'FM-MNC','FM-MSNC')
Table4


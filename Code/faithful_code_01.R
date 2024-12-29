library(mvtnorm)
library(tmvtnorm)
library(mclust)
library(MomTrunc)
library(mixtools)

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

#choose best model function
B_model = function(Y,Ycm,cen,k,tol = 1e-5,neq.S = T,max.iter=1000,tol.FG=1e-3,max.iter.FG=10, per=1)
{
  table = matrix(NA, k,14)
  model = c()
  for(g in 1:k){
    clus=kmeans(Y,g)$cluster
    
    fit_VVV = VVV.GMIXCMB.ECM.new(Ycm, Y.knn=Y, cen, g, true.clus=clus, tol = 1e-5, neq.S = T, max.iter, per)
    fit_EEE = EEE.GMIXCMB.ECM.new(Ycm, Y.knn=Y, cen, g, true.clus=clus, tol = 1e-5, neq.S = T, max.iter, per)
    fit_EEV = EEV.GMIXCMB.ECM.new(Ycm, Y.knn=Y, cen, g, true.clus=clus, tol = 1e-5, neq.S = T, max.iter, per)
    fit_EVE = EVE.GMIXCMB.ECM.new(Ycm, Y.knn=Y, cen, g, true.clus=clus, tol = 1e-5, neq.S = T, max.iter, tol.FG, max.iter.FG, per)
    fit_EVV = EVV.GMIXCMB.ECM.new(Ycm, Y.knn=Y, cen, g, true.clus=clus, tol = 1e-5, neq.S = T, max.iter, per)
    fit_VEE = VEE.GMIXCMB.ECM.new(Ycm, Y.knn=Y, cen, g, true.clus=clus, tol = 1e-5, neq.S = T, max.iter, per)
    fit_VEV = VEV.GMIXCMB.ECM.new(Ycm, Y.knn=Y, cen, g, true.clus=clus, tol = 1e-5, neq.S = T, max.iter, per)
    fit_VVE = VVE.GMIXCMB.ECM.new(Ycm, Y.knn=Y, cen, g, true.clus=clus, tol = 1e-5, neq.S = T, max.iter, tol.FG, max.iter.FG, per)
    
    fit_VEI = VEI.GMIXCMB.ECM.new(Ycm, Y.knn=Y, cen, g, true.clus=clus, tol = 1e-5, neq.S = T, max.iter, per)
    fit_EEI = EEI.GMIXCMB.ECM.new(Ycm, Y.knn=Y, cen, g, true.clus=clus, tol = 1e-5, neq.S = T, max.iter, per)
    fit_VVI = VVI.GMIXCMB.ECM.new(Ycm, Y.knn=Y, cen, g, true.clus=clus, tol = 1e-5, neq.S = T, max.iter, per)
    fit_EVI = EVI.GMIXCMB.ECM.new(Ycm, Y.knn=Y, cen, g, true.clus=clus, tol = 1e-5, neq.S = T, max.iter, per)
    
    fit_VII = VII.GMIXCMB.ECM.new(Ycm, Y.knn=Y, cen, g, true.clus=clus, tol = 1e-5, neq.S = T, max.iter, per)
    fit_EII = EII.GMIXCMB.ECM.new(Ycm, Y.knn=Y, cen, g, true.clus=clus, tol = 1e-5, neq.S = T, max.iter, per)
    a=(c(fit_VVV$model.inf[[5]],fit_EEE$model.inf[[5]],fit_EEV$model.inf[[5]],fit_EVE$model.inf[[5]],
         fit_EVV$model.inf[[5]],fit_VEE$model.inf[[5]],fit_VEV$model.inf[[5]],fit_VVE$model.inf[[5]],
         fit_VEI$model.inf[[5]],fit_EEI$model.inf[[5]],fit_VVI$model.inf[[5]],fit_EVI$model.inf[[5]],
         fit_VII$model.inf[[5]],fit_EII$model.inf[[5]]))
    names(a) = paste(c("VVV", "EEE", "EEV","EVE","EVV","VEE","VEV","VVE","VEI","EEI","VVI","EVI","VII","EII"),g,sep = ",")
    model = c(model,a)
    table[g,] = a
  }
  result = sort(model,decreasing = F)[1:3]
  colnames(table) = c("VVV", "EEE", "EEV","EVE","EVV","VEE","VEV","VVE","VEI","EEI","VVI","EVI","VII","EII")
  rownames(table) = c(1:k)
  return(list(Best_model = result, BIC = table))
}

#=============== Best model for g = 1 to 5 without scale ==========================#
set.seed(0307)
Y = (faithful)
NC=gener.cen.na(Data = Y,cen.type = c(0,0),cen.rate=c(0,0),na.rate = 0)
Ycm = NC$Data
cen = NC$cen


Best5_B = B_model(Y,Ycm,cen,k=5,tol = 1e-5,neq.S = T,max.iter=1000,tol.FG=1e-3,max.iter.FG=10, per=1)
Bestmodel = Best5_B$Best_model

g=3;clus0=kmeans(Y,g)$cluster
fit0 = EEE.GMIXCMB.ECM.new(Y, cen, g, clus0, tol = 1e-5, neq.S = T, max.iter=1000, per=1)
n=nrow(Y)
p.clus0 = fit0$post.clus

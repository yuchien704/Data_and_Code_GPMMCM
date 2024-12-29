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
B_model = function(Y,Ycm,cen,k,tol = 1e-5,neq.S = T,max.iter=500,tol.FG=1e-3,max.iter.FG=10, per=1)
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

#=============== Best model for g = 1 to g = 5 without scale ==========================#
set.seed(0307)
Y = (faithful)
NC=gener.cen.na(Data = Y,cen.type = c(2,1),cen.rate=c(0.1,0.1),na.rate = 0.15)
Ycm = NC$Data
cen = NC$cen


Best5_B = B_model(Y,Ycm,cen,k=5,tol = 1e-5,neq.S = T,max.iter=1000,tol.FG=1e-3,max.iter.FG=10, per=100)
BIC = Best5_B$BIC
#================================ plot faithful=====================================#

x=faithful[,1]
y=faithful[,2]
x_cm=Ycm[,1]
y_cm=Ycm[,2]
n=nrow(Ycm)
p=ncol(Ycm)
cen_11 = which(cen[,1]==2)
cen_21 = which(cen[,2]==1)
cen_11_value = Ycm[which(cen[,1]==2)[1],1] 
cen_21_value = Ycm[which(cen[,2]==1)[1],2] 
na.posi = which(rowSums(is.na(cen))==1) 


#contour
set.seed(0307)
g=3;clus=kmeans(Y,g)$cluster
fit1 = EVI.GMIXCMB.ECM.new(Ycm, cen, g, clus, tol = 1e-5, neq.S = T, max.iter=1000, per=1)
m=30
tmp1=range(na.omit(Ycm[,1]))
tmp2=range(na.omit(Ycm[,2]))
xx=expand.grid(x1<-seq(tmp1[1],tmp1[2],length=m),x2<-seq(tmp2[1],tmp2[2],length=m))
den=0
for(i in 1:g) den=den+mvtnorm::dmvnorm(xx, mean=fit1$para$mu[,i], sigma=fit1$para$Sigma[,,i])
den=matrix(den,m,m)
p.clus = fit1$post.clus
y.pred = fit1$post.pred

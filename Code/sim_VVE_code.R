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

source(paste(SPATH, '/Function/Both_gmixcm/VVV_GMIXCMB_1.R', sep=''))
source(paste(SPATH, '/Function/Both_gmixcm/EVV_GMIXCMB_1.R', sep=''))
source(paste(SPATH, '/Function/Both_gmixcm/EVE_GMIXCMB_1.R', sep=''))
source(paste(SPATH, '/Function/Both_gmixcm/VVE_GMIXCMB_1.R', sep=''))
source(paste(SPATH, '/Function/Both_gmixcm/VEV_GMIXCMB_1.R', sep=''))
source(paste(SPATH, '/Function/Both_gmixcm/EEV_GMIXCMB_1.R', sep=''))


gen.nmix = function(ni, mu, Sigma)                # generate Y directly
{
  g = ncol(mu); 
  Y = list(g)
  for(i in 1: g) Y[[i]] = rmvnorm(ni[i], mu[,i], Sigma[,,i]) 
  y = NULL
  for(i in 1: g) y = rbind(y, Y[[i]])
  Y.gen = cbind(y)
  return(Y.gen) 
}

sim.GPMMCM = function(n, p, G, k, max.sim, cen.type, cen.rate, na.rate, tol ,max.iter, tol.FG, max.iter.FG, per)
{
  begin = proc.time()[1]
  # ----- initial values ----------- #
  ni=c(n/3, n/3, n/3)
  w=c(1/3,1/3,1/3)
  mu = matrix(c(-6, -5, 6, -1, -4, 4), p, G)
  Sigma = array(NA, dim = c(p, p, G))
  lambda1 = 16
  lambda2 = 9
  lambda3 = 4
  theta1 = pi/4
  D = matrix(c(cos(theta1),-sin(theta1),sin(theta1),cos(theta1)),p,p)
  
  A1 = diag(c(2,0.5))
  A2 = diag(c(3,1/3))
  A3 = diag(c(4,1/4))
  Sigma[, , 1] = lambda1*D%*%A1%*%t(D)
  Sigma[, , 2] = lambda2*D%*%A2%*%t(D)
  Sigma[, , 3] = lambda3*D%*%A3%*%t(D)
  
  g=1:k; nn=length(g)
  AIC = BIC = list()
  select_model = select_group = bestModel = c()
  select_model_AIC = select_group_AIC = bestModel_AIC = c()
  
  for(j in 1:max.sim){
    ni=c(n/3, n/3, n/3)
    w=c(1/3,1/3,1/3)
    mu = matrix(c(-6, -5, 6, -1, -4, 4), p, G)
    Sigma = array(NA, dim = c(p, p, G))
    lambda1 = 16
    lambda2 = 9
    lambda3 = 4
    theta1 = pi/4
    D = matrix(c(cos(theta1),-sin(theta1),sin(theta1),cos(theta1)),p,p)
    A1 = diag(c(2,0.5))
    A2 = diag(c(3,1/3))
    A3 = diag(c(4,1/4))
    Sigma[, , 1] = lambda1*D%*%A1%*%t(D)
    Sigma[, , 2] = lambda2*D%*%A2%*%t(D)
    Sigma[, , 3] = lambda3*D%*%A3%*%t(D)
    Y = matrix(NA, n, p)
    Ycm = cen = matrix(NA, n, p) 
    Data = list()
    Y = gen.nmix(ni , mu, Sigma) 
    Data = gener.cen.na(Y, cen.type = c(1, 2),cen.rate=c(0.1,0.1),na.rate = 0.1) 
    Ycm = as.matrix(Data$Data)
    cen = Data$cen
    VEE.BIC=c();EVV.BIC=c();VVE.BIC=c();EVE.BIC=c();EEE.BIC=c();EEI.BIC=c();EEV.BIC=c()
    EII.BIC=c();EVI.BIC=c();VEI.BIC=c();VEV.BIC=c();VII.BIC=c();VVI.BIC=c();VVV.BIC=c()
    VEE.AIC=c();EVV.AIC=c();VVE.AIC=c();EVE.AIC=c();EEE.AIC=c();EEI.AIC=c();EEV.AIC=c()
    EII.AIC=c();EVI.AIC=c();VEI.AIC=c();VEV.AIC=c();VII.AIC=c();VVI.AIC=c();VVV.AIC=c()
    
    cat('This is the ', j, ' time. \n')
    for(i in 1:nn)
    {
      if(i==1){clus=rep(1,n)}else{clus=kmeans(Y, g[i])$cluster}
      repeat{
        fit.GMIXCM.VVV= try(VVV.GMIX.MVNCM.B.ECM(Ycm, cen, g[i], clus, tol=1e-5, neq.S=T, max.iter=1000, per=1), silent=T)
        if(class(fit.GMIXCM.VVV) != "try-error") break
        else{
          fit.GMIXCM.VVV = try(VVV.GMIX.MVNCM.B.ECM.renew(Y, Ycm, cen, g[i], tol=1e-5, neq.S=T, max.iter=1000, per=1), silent=T)
          if(class(fit.GMIXCM.VVV) != "try-error") break
        }
      }
      
      fit.GMIXCM.EEE=EEE.GMIX.MVNCM.B.ECM(Ycm, cen, g[i], clus, tol=1e-5, neq.S=F, max.iter=1000, per=1)
      fit.GMIXCM.VEE=VEE.GMIX.MVNCM.B.ECM(Ycm, cen, g[i], clus, tol=1e-5, neq.S=T, max.iter=1000, per=1)
      repeat{
        fit.GMIXCM.EEV= try(EEV.GMIX.MVNCM.B.ECM(Ycm, cen, g[i], clus, tol=1e-5, neq.S=T, max.iter=1000, per=1), silent=T)
        if(class(fit.GMIXCM.EEV) != "try-error") break
        else{
          fit.GMIXCM.EEV = try(EEV.GMIX.MVNCM.B.ECM.renew(Y, Ycm, cen, g[i], tol=1e-5, neq.S=T, max.iter=1000, per=1), silent=T)
          if(class(fit.GMIXCM.EEV) != "try-error") break
        }
      }
      repeat{
        fit.GMIXCM.VEV= try(VEV.GMIX.MVNCM.B.ECM(Ycm, cen, g[i], clus, tol=1e-5, neq.S=T, max.iter=1000, per=1), silent=T)
        if(class(fit.GMIXCM.VEV) != "try-error") break
        else{
          fit.GMIXCM.VEV = try(VEV.GMIX.MVNCM.B.ECM.renew(Y, Ycm, cen, g[i], tol=1e-5, neq.S=T, max.iter=1000, per=1), silent=T)
          if(class(fit.GMIXCM.VEV) != "try-error") break
        }
      }
      repeat{
        fit.GMIXCM.EVV= try(EVV.GMIX.MVNCM.B.ECM(Ycm, cen, g[i], clus, tol=1e-5, neq.S=T, max.iter=1000, per=1), silent=T)
        if(class(fit.GMIXCM.EVV) != "try-error") break
        else{
          fit.GMIXCM.EVV = try(EVV.GMIX.MVNCM.B.ECM.renew(Y, Ycm, cen, g[i], tol=1e-5, neq.S=T, max.iter=1000, per=1), silent=T)
          if(class(fit.GMIXCM.EVV) != "try-error") break
        }
      }
      ##FG
      repeat{
        fit.GMIXCM.EVE = try(EVE.GMIX.MVNCM.B.ECM(Ycm, cen, g[i], clus, tol = 1e-5, neq.S = T, max.iter, tol.FG=1e-5, max.iter.FG = 100, per=1), silent=T)
        if(class(fit.GMIXCM.EVE) != "try-error") break
        else{
          fit.GMIXCM.EVE = try(EVE.GMIX.MVNCM.B.ECM.renew(Y, Ycm, cen, g[i], tol = 1e-5, neq.S = T, max.iter=1000, tol.FG=1e-5, max.iter.FG = 100, per=1), silent=T)
          if(class(fit.GMIXCM.EVE) != "try-error") break
        }
      }
      repeat{
        fit.GMIXCM.VVE = try(VVE.GMIX.MVNCM.B.ECM(Ycm, cen, g[i], clus, tol = 1e-5, neq.S = T, max.iter, tol.FG=1e-5, max.iter.FG = 100, per=1), silent=T)
        if(class(fit.GMIXCM.VVE) != "try-error") break
        else{
          fit.GMIXCM.VVE = try(VVE.GMIX.MVNCM.B.ECM.renew(Y, Ycm, cen, g[i], tol = 1e-5, neq.S = T, max.iter=1000, tol.FG=1e-5, max.iter.FG = 100, per=1), silent=T)
          if(class(fit.GMIXCM.VVE) != "try-error") break
        }
      }

      fit.GMIXCM.EEI=EEI.GMIX.MVNCM.B.ECM(Ycm, cen, g[i], clus, tol=1e-5, neq.S=F, max.iter=1000, per=1)
      fit.GMIXCM.VEI=VEI.GMIX.MVNCM.B.ECM(Ycm, cen, g[i], clus, tol=1e-5, neq.S=T, max.iter=1000, per=1)
      fit.GMIXCM.EVI=EVI.GMIX.MVNCM.B.ECM(Ycm, cen, g[i], clus, tol=1e-5, neq.S=T, max.iter=1000, per=1)
      fit.GMIXCM.VVI=VVI.GMIX.MVNCM.B.ECM(Ycm, cen, g[i], clus, tol=1e-5, neq.S=T, max.iter=1000, per=1)
      
      fit.GMIXCM.EII=EII.GMIX.MVNCM.B.ECM(Ycm, cen, g[i], clus, tol=1e-5, neq.S=F, max.iter=1000, per=1)
      fit.GMIXCM.VII=VII.GMIX.MVNCM.B.ECM(Ycm, cen, g[i], clus, tol=1e-5, neq.S=T, max.iter=1000, per=1)
      
      ###########BIC#########
      EEE.BIC[i]=fit.GMIXCM.EEE$model.inf[[5]]
      EEI.BIC[i]=fit.GMIXCM.EEI$model.inf[[5]]
      EEV.BIC[i]=fit.GMIXCM.EEV$model.inf[[5]]
      EII.BIC[i]=fit.GMIXCM.EII$model.inf[[5]]
      EVI.BIC[i]=fit.GMIXCM.EVI$model.inf[[5]]
      VEI.BIC[i]=fit.GMIXCM.VEI$model.inf[[5]]
      VEV.BIC[i]=fit.GMIXCM.VEV$model.inf[[5]]
      VII.BIC[i]=fit.GMIXCM.VII$model.inf[[5]]
      VVI.BIC[i]=fit.GMIXCM.VVI$model.inf[[5]]
      VVV.BIC[i]=fit.GMIXCM.VVV$model.inf[[5]]
      VEE.BIC[i]=fit.GMIXCM.VEE$model.inf[[5]]
      EVV.BIC[i]=fit.GMIXCM.EVV$model.inf[[5]]
      #FG#
      VVE.BIC[i]=fit.GMIXCM.VVE$model.inf[[5]]  
      EVE.BIC[i]=fit.GMIXCM.EVE$model.inf[[5]]
    }
    BIC[[j]] = matrix(c(EII.BIC,VII.BIC,EEI.BIC,VEI.BIC,EVI.BIC,VVI.BIC,EEE.BIC,VEE.BIC,EVE.BIC,VVE.BIC,EEV.BIC,VEV.BIC,EVV.BIC,VVV.BIC),
                      ncol=14,dimnames = list(c(1:k),c("EII","VII","EEI","VEI","EVI","VVI","EEE","VEE","EVE","VVE","EEV","VEV","EVV","VVV"))) 
    AIC[[j]] = matrix(c(EII.AIC,VII.AIC,EEI.AIC,VEI.AIC,EVI.AIC,VVI.AIC,EEE.AIC,VEE.AIC,EVE.AIC,VVE.AIC,EEV.AIC,VEV.AIC,EVV.AIC,VVV.AIC),
                      ncol=14,dimnames = list(c(1:k),c("EII","VII","EEI","VEI","EVI","VVI","EEE","VEE","EVE","VVE","EEV","VEV","EVV","VVV"))) 
    # BIC select
    Modelname = c("EII","VII","EEI","VEI","EVI","VVI","EEE","VEE","EVE","VVE","EEV","VEV","EVV","VVV")
    select_group[j] = which.min(BIC[[j]])%%k
    if(select_group[j]==0){
      select_group[j]=k
      select_model[j] = Modelname[which.min(BIC[[j]])%/%k]
    }else{
      select_model[j] = Modelname[which.min(BIC[[j]])%/%k+1]
    }
    bestModel[j] = paste(select_model[j], select_group[j],sep = ",")
    # AIC select
    select_group_AIC[j] = which.min(AIC[[j]])%%k
    if(select_group_AIC[j]==0){
      select_group_AIC[j]=k
      select_model_AIC[j] = Modelname[which.min(AIC[[j]])%/%k]
    }else{
      select_model_AIC[j] = Modelname[which.min(AIC[[j]])%/%k+1]
    }
    bestModel_AIC[j] = paste(select_model_AIC[j], select_group_AIC[j],sep = ",")
  }
  end = proc.time()[1]
  run.sec = end - begin
  correct.rate = mean(bestModel=="VVE,3")
  bestModel = matrix(bestModel, ncol = 10, byrow= T) 
  correct.rate_AIC = mean(bestModel_AIC=="VVE,3")
  bestModel_AIC = matrix(bestModel_AIC, ncol = 10, byrow= T)
  return(list(run.sec = run.sec, correct.rate = correct.rate, BIC = BIC, bestModel = bestModel, correct.rate_AIC = correct.rate_AIC, AIC = AIC, bestModel_AIC = bestModel_AIC))
}

n=450; p=2; G=3;k=5; max.sim=20; cen.type = c(1, 2);cen.rate=c(0.1,0.1);na.rate = 0.1 
tol=1e-5; max.iter=1000; tol.FG=1e-3; max.iter.FG = 10; per=1
SIM_VVE_CM_1 = sim.GPMMCM(n, p, G, k, max.sim, cen.type, cen.rate, na.rate, tol, max.iter, tol.FG, max.iter.FG , per)
SIM_VVE_CM_2 = sim.GPMMCM(n, p, G, k, max.sim, cen.type, cen.rate, na.rate, tol, max.iter, tol.FG, max.iter.FG , per)
SIM_VVE_CM_3 = sim.GPMMCM(n, p, G, k, max.sim, cen.type, cen.rate, na.rate, tol, max.iter, tol.FG, max.iter.FG , per)
SIM_VVE_CM_4 = sim.GPMMCM(n, p, G, k, max.sim, cen.type, cen.rate, na.rate, tol, max.iter, tol.FG, max.iter.FG , per)
SIM_VVE_CM_5 = sim.GPMMCM(n, p, G, k, max.sim, cen.type, cen.rate, na.rate, tol, max.iter, tol.FG, max.iter.FG , per)

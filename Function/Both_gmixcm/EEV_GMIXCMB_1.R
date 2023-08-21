library(mvtnorm)
library(tmvtnorm)
library(mclust)
library(MomTrunc)

vech.posi=function(dim) cbind(rep(1:dim, 1:dim), unlist(mapply(':', 1, 1:dim)))

log.like_B = function(Ycm, cen, cen1, w, mu, Sigma, yc, Oj, Cj, Fj, cen.subj, cen.posi)
{
  n = nrow(Ycm) 
  p = ncol(Ycm) 
  g = ncol(mu)
  Ip = diag(p)
  na = is.na(cen1) 
  po = p - rowSums(na)
  
  w.den = matrix(NA, n, g)
  det.f = rep(NA, n)
  delta.f = matrix(NA, n, g)
  cdf = matrix(1, n, g)
  Ycm.imp = Ycm
  Ycm.imp[is.na(Ycm)] = 999
  
  nc = length(cen.subj) #設限值總列數
  cj = numeric(nc) 
  ccj = rep(0, n) 
  if(nc != 0){ #若有設限
    for(j in 1: nc) cj[j] = length(cen.posi[[cen.subj[j]]]) #某列設限值個數
    ccj[cen.subj] = cj #每列設限值的個數
  }#else{ cj = 0 ; ccj = rep(0, n)  } #沒有設限
  cumsum.c = cumsum(cj) #累加cj (計算yc.hat)
  fj = po - ccj #每列完全看得到的個數
  for(i in 1:g){
    cent = t(Ycm.imp) - mu[,i]
    if(nc!=0){
      for(j in 1:nc){
        subj = cen.subj[j]
        if(j == 1){ idx = 1: cumsum.c[j]
        } else idx = (cumsum.c[(j-1)]+1): cumsum.c[j]
        if(cj[j] == po[subj]){
          det.f[subj] = 1
          delta.f[subj,i] = 0
          mu.cf = Cj[[subj]]%*%Oj[[subj]]%*%mu[,i]
          Sig.cc.f = Cj[[subj]] %*%Oj[[subj]] %*% Sigma[,,i] %*% t(Oj[[subj]])%*% t(Cj[[subj]])
          
          if(cj[[j]] == 1)
          {
            h = t(cen)[,subj] # 第j種設限情況
            if(h[which(h!=0)]==1){ #左設限
              cdf[subj,i] = pnorm(yc[idx], mean=mu.cf, sd=sqrt(Sig.cc.f))
            } else cdf[subj,i] = 1-pnorm(yc[idx], mean=mu.cf, sd=sqrt(Sig.cc.f)) #右設限
          } else{
            d0 = t(cen)[,subj]
            d0 = d0[which(d0!=0)]
            bound = yc[idx] #censored value
            A=rep(-Inf,length(d0))
            B=rep(Inf,length(d0))
            for(t in 1:length(d0)){
              if(d0[t]==1){
                B[t]=bound[t]
              }else{
                A[t]=bound[t]
              }
            } #censored bound
            cdf[subj,i] = pmvnorm(lower=A, upper=B, mean=c(mu.cf), sigma=Sig.cc.f)[1]
          }
        } else{
          cent.subj = cent[,subj] 
          FOSOF = Fj[[subj]] %*% Oj[[subj]] %*% Sigma[,,i] %*% t(Oj[[subj]]) %*% t(Fj[[subj]])
          det.f[subj] = det(FOSOF)
          foFOSOFof = t(Fj[[subj]] %*% Oj[[subj]]) %*% solve(FOSOF) %*% Fj[[subj]] %*% Oj[[subj]]
          delta.f[subj,i] = colSums(cent.subj * (foFOSOFof %*% cent.subj))
          mu.cf = Cj[[subj]] %*% Oj[[subj]] %*% (mu[,i] + Sigma[,,i] %*% foFOSOFof %*% cent.subj)
          Sig.cc.f = Cj[[subj]] %*% Oj[[subj]] %*% (Ip - Sigma[,,i] %*% foFOSOFof) %*% Sigma[,,i] %*% t(Oj[[subj]]) %*% t(Cj[[subj]])
          if(cj[[j]] == 1)
          {
            h = t(cen)[,subj] # 第j種設限情況
            if(h[which(h!=0)]==1){ #左設限
              cdf[subj,i] = pnorm(yc[idx], mean=mu.cf, sd=sqrt(Sig.cc.f))
            } else cdf[subj,i] = 1-pnorm(yc[idx], mean=mu.cf, sd=sqrt(Sig.cc.f)) #右設限
          } else{
            d0 = t(cen)[,subj]
            d0 = d0[which(d0!=0)]
            bound = yc[idx] #censored value
            A=rep(-Inf,length(d0))
            B=rep(Inf,length(d0))
            for(t in 1:length(d0)){
              if(d0[t]==1){
                B[t]=bound[t]
              }else{
                A[t]=bound[t]
              }
            } #censored bound
            cdf[subj,i] = pmvnorm(lower=A, upper=B, mean=c(mu.cf), sigma=Sig.cc.f)[1]
          }
        }
      }
      for(j in (1:n)[-cen.subj]){
        cent.subj = cent[,j] 
        OSO = Oj[[j]] %*% Sigma[,,i] %*% t(Oj[[j]])
        det.f[j] = det(OSO)
        Soo = t(Oj[[j]]) %*% solve(OSO) %*% Oj[[j]]
        delta.f[j,i] = colSums(cent.subj * (Soo %*% cent.subj))
      }}else{
        for(j in 1:n){
          cent.subj = cent[,j] 
          OSO = Oj[[j]] %*% Sigma[,,i] %*% t(Oj[[j]])
          det.f[j] = det(OSO)
          Soo = t(Oj[[j]]) %*% solve(OSO) %*% Oj[[j]]
          delta.f[j,i] = colSums(cent.subj * (Soo %*% cent.subj))
        }
      }
    w.den[,i] =  log(w[i]) - 0.5*fj*log(2*pi) - 0.5*log(det.f) - 0.5*delta.f[,i] + log(cdf[,i])
  }
  return(w.den)
}
EEV.GMIX.MVNCM.B.ECM.renew = function(Y, Ycm, cen, g, tol = 1e-5, neq.S = T, max.iter, per)
{
  begin = proc.time()[1]
  iter =0
  # Basic setting
  n = nrow(Ycm) #樣本數
  p = ncol(Ycm) #參數個數
  
  # Specify initial values
  w = rep(NA, g)
  mu = matrix(NA, p, g)
  Sigma = array(NA, dim = c(p,p,g))
  #W=array(NA, dim = c(p,p,g))
  lambda= rep(1, g) 
  
  omega=array(NA, dim = c(p,p,g))
  L=array(NA, dim = c(p,p,g))
  if(g==1){clus=rep(1,n)}else{clus=kmeans(Y, g)$cluster}
  
  for(i in 1:g)
  {
    w[i] = nrow(Ycm[clus == i,]) / n
    mu[,i] = colMeans(Ycm[clus == i,],na.rm = T)
    Sigma[,,i] =  cov(na.omit(Ycm[clus == i,]))
    #cent = t(Ycm) - mu[,i]
  }
  
  np = n*p      #總樣本數
  Ip = diag(p)  #pxp的單位矩陣
  cen1 =ifelse(cen,1,0)
  na = is.na(cen1) #判斷cen是否遺失
  po = p - rowSums(na) #每一列看得見的維度
  na.posi = cen.posi = Oj = Fj = Cj = as.list(n)
  for(j in 1: n){
    Oj[[j]] = Ip
    Ipo = diag(po[j])
    Fj[[j]] = Ipo
    Cj[[j]] = NULL
    if(sum(is.na(Ycm[j,])!=0)){
      na.posi[[j]] = which(is.na(Ycm[j,]))
      Oj[[j]] = matrix(Ip[-na.posi[[j]],], ncol = p)}
    if(sum(cen[j, ], na.rm=T) != 0){
      cen.posi[[j]] = which(na.omit(cen[j, ]) != 0)
      Fj[[j]] = matrix(Ipo[-cen.posi[[j]],], ncol = po[j])
      Cj[[j]] = matrix(Ipo[cen.posi[[j]], ], ncol = po[j])}
  }
  
  ycm = as.vector(t(Ycm)) #每個觀測者觀測資料的向量
  na.ind = which(is.na(as.vector(t(cen1)))) #每個觀測者觀測資料向量的遺失位置
  
  if(length(na.ind)!=0){ #若有遺失
    yo = ycm[-na.ind]
    cen.ind = which(as.vector(t(cen1))[-na.ind]==1) #每個觀測者觀測資料向量的設限位置
  }else { #沒有遺失
    yo = ycm ; cen.ind = which(as.vector(t(cen1))==1)}
  
  if(length(cen.ind)!=0){ #若有設限
    yf = yo[-cen.ind]
    yc = yo[cen.ind]
  }else { #沒有設限
    yf = yo ; yc = 0}
  
  cumsum.p = cumsum(rep(p, n)) #累加p (後面M-step 計算sigma)
  cumsum.o = cumsum(po) #累加po 
  no = sum(po) #每列看得到的個數加總
  na.subj = which(rowSums(na) != 0) #哪幾列遺失
  cen.subj = which(apply(cen1, 1, sum, na.rm=T) != 0) #哪幾列設限
  #==== 計算cj,fj ====#
  nc = length(cen.subj) #設限值總列數
  cj = numeric(nc) 
  ccj = rep(0, n) 
  if(nc != 0){ #若有設限
    for(j in 1: nc) cj[j] = length(cen.posi[[cen.subj[j]]]) #某列設限值個數
    ccj[cen.subj] = cj #每列設限值的個數
  }
  cumsum.c = cumsum(cj) #累加cj (計算yc.hat)
  fj = po - ccj #每列完全看得到的個數
  cumsum.f = cumsum(fj) #累加fj (計算 log-likelihood)
  nf = sum(fj) #每列完全看得到的個數加總(計算observed log-likelihood)
  #===================#
  if(length(na.ind)!=0){ #若有遺失
    TO = diag(np)[-na.ind, ] #輔助矩陣Oj
    TM = diag(np)[na.ind, ] #輔助矩陣Mj
  }else{ #沒有遺失
    TO = diag(np)}
  if(length(cen.ind)!=0){ #若有設限
    TF = diag(no)[-cen.ind, ] #輔助矩陣Fj
    TC = diag(no)[cen.ind, ] #輔助矩陣Cj
  }else{ #沒有設限
    TF = diag(no)}
  num.na = length(na.ind) #遺失值個數
  num.cen = length(cen.ind) #設限值總個數
  vechS = vech.posi(p) #pxp矩陣下三角位置
  
  # observed data log-likelihood
  
  w.den = log.like_B(Ycm, cen, cen1, w, mu, Sigma, yc, Oj, Cj, Fj, cen.subj, cen.posi)
  
  max.wden = apply(w.den, 1, max)
  w.den = exp(w.den - max.wden)
  indv.den = rowSums(w.den)
  log.indv.den = log(indv.den) + max.wden
  old.loglik = sum(log.indv.den)  
  
  cat(paste(rep('-', 20), sep = '', collapse = ''), 'EEV.GMIX.MVNCM: ', 'g =', g, paste(rep('-', 20), sep = '', collapse = ''), '\n')
  #cat("iter =", iter, "\t loglik =", old.loglik, "\n") 
  
  repeat
  {
    iter = iter+1
    
    # E-step
    Z = w.den / indv.den
    ni = colSums(Z) 
    w = ni / n
    
    y.hat = matrix(NA, nrow=np, ncol= g)
    y2.hat = array(0, dim=c(np, np, g))
    yc.hat = matrix(NA, nrow=num.cen, ncol=1)
    yc2.hat = matrix(0, nrow=num.cen, ncol=num.cen)
    sum.Ome = array(NA, dim=c(p, p, g))
    
    for(i in 1:g)
    {
      Tmu = rep(mu[,i], n) #mu向量
      TSig = kronecker(diag(n), Sigma[,,i]) #Sigma矩陣
      if(num.na!=0){#若有遺失
        TSig.oo = TSig[-na.ind, -na.ind] # Oj%*%Sigma%*%T(Oj) 以下敘述以SOO表示
        TSig.mo = TSig[na.ind, -na.ind]  # Mj%*%Sigma%*%T(Oj) 
        TSig.mm = TSig[na.ind, na.ind]   # Mj%*%Sigma%*%T(Mj) 
        if(num.na == 1) TSig.mo = t(TSig.mo) 
      }else{ #沒有遺失
        TSig.oo = TSig }
      if(num.cen!=0){ #若有設限
        TSig.ff = TSig.oo[-cen.ind, -cen.ind] # Fj%*%SOO%*%t(Fj)
        TSig.cf = TSig.oo[cen.ind, -cen.ind]  # Cj%*%SOO%*%t(Fj)
        TSig.cc = TSig.oo[cen.ind, cen.ind]   # Fj%*%SOO%*%t(Cj)
        if(num.cen == 1){
          TSig.cf = t(TSig.cf); TSig.cc = as.matrix(TSig.cc)}
      }else{ #沒有設限
        TSig.ff = TSig.oo}
      TSig.ff.inv = matrix(0, nf, nf)
      for(j in 1:n){
        if(j == 1){
          idx.f = 1: cumsum.f[1]
        }else{
          if((cumsum.f[j-1]+1)>cumsum.f[j]){
            j = j + 1
          }else{
            idx.f = (cumsum.f[j-1]+1): cumsum.f[j]
          }}
        TSig.ff.inv[idx.f,idx.f] = solve(TSig.ff[idx.f,idx.f])
      }
      #TSig.ff.inv = solve(TSig.ff) #the inverse of Fj%*%SOO%*%t(Fj)
      
      if(num.na!=0)
      { #遺失情況
        if(num.cen!=0){ #有遺失、有設限(M+C,M+C+F)
          yf.cent = yf - Tmu[-na.ind][-cen.ind] 
          mu.cf = Tmu[-na.ind][cen.ind] + TSig.cf %*% TSig.ff.inv %*% yf.cent
          Sig.cc.f = TSig.cc - TSig.cf %*% TSig.ff.inv %*% t(TSig.cf)
          #evaluate right or left censored
          d0 = t(cen)[,cen.subj]
          d0 = d0[which(d0!=0)]
          bound = yc  #censored value
          A=rep(-Inf,length(d0))
          B=rep(Inf,length(d0))
          for(t in 1:length(d0)){
            if(d0[t]==1){
              B[t]=bound[t]
            }else{
              A[t]=bound[t]
            }
          }
          ycj.hat = meanvarTMD(lower = A[1:cumsum.c[1]],upper = B[1:cumsum.c[1]] ,mu = mu.cf[1:cumsum.c[1]] ,Sigma = Sig.cc.f[1:cumsum.c[1],1:cumsum.c[1]],dist = 'normal')
          yc.hat[1:cumsum.c[1]] = ycj.hat$mean
          yc2.hat[1:cumsum.c[1], 1:cumsum.c[1]] = ycj.hat$EYY
          if(nc != 1){
            for(j in 2: nc){
              ycj.hat = meanvarTMD(lower=A[(cumsum.c[j-1]+1):cumsum.c[j]], upper=B[(cumsum.c[j-1]+1):cumsum.c[j]] ,mu=c(mu.cf[(cumsum.c[j-1]+1):cumsum.c[j]]) ,Sigma=Sig.cc.f[(cumsum.c[j-1]+1):cumsum.c[j],(cumsum.c[j-1]+1):cumsum.c[j]],dist = 'normal')
              yc.hat[(cumsum.c[j-1]+1):cumsum.c[j]] = ycj.hat$mean
              yc2.hat[(cumsum.c[j-1]+1):cumsum.c[j], (cumsum.c[j-1]+1):cumsum.c[j]] = ycj.hat$EYY
            }}
          TSig.oo.inv = matrix(0, no, no)
          for(j in 1:n){
            if(j == 1){
              idx1 = 1: cumsum.o[1]
            }else{
              idx1 = (cumsum.o[j-1]+1): cumsum.o[j]
            }
            TSig.oo.inv[idx1,idx1] = solve(TSig.oo[idx1,idx1])
          }
          #TSig.oo.inv = solve(TSig.oo)
          yhat.cent = t(TF)%*%yf + t(TC)%*%yc.hat - Tmu[-na.ind]
          yo.cent = t(TF)%*%yf - Tmu[-na.ind]
          ym.hat = Tmu[na.ind] + TSig.mo %*% TSig.oo.inv %*% yhat.cent
          y.hat[,i] =  t(TO)%*%t(TF)%*%yf + t(TO)%*%t(TC)%*%yc.hat + t(TM)%*%ym.hat
          Sig.mm.o = TSig.mm - TSig.mo %*% TSig.oo.inv %*% t(TSig.mo)
          ym2.hat = Tmu[na.ind]%*%t(Tmu[na.ind]) + Tmu[na.ind]%*%t(yo.cent)%*%TSig.oo.inv%*%t(TSig.mo) + Tmu[na.ind]%*%t(yc.hat)%*%TC%*%TSig.oo.inv%*%t(TSig.mo)+ TSig.mo %*% TSig.oo.inv %*% yo.cent%*%t(Tmu[na.ind]) +TSig.mo %*% TSig.oo.inv %*% yo.cent%*%t(yo.cent)%*%TSig.oo.inv%*%t(TSig.mo) + TSig.mo %*% TSig.oo.inv %*% yo.cent%*%t(yc.hat)%*%TC%*%TSig.oo.inv%*%t(TSig.mo)+ TSig.mo %*% TSig.oo.inv %*% t(TC)%*%yc.hat%*%t(Tmu[na.ind]) + TSig.mo %*% TSig.oo.inv %*% t(TC)%*%yc.hat%*%t(yo.cent)%*%TSig.oo.inv%*%t(TSig.mo)+TSig.mo %*% TSig.oo.inv %*% t(TC)%*%yc2.hat%*%TC%*%TSig.oo.inv%*%t(TSig.mo)+ Sig.mm.o
          #ym2.hat = Tmu[na.ind]%*%t(Tmu[na.ind]) + Tmu[na.ind]%*%t(yhat.cent)%*%TSig.oo.inv%*%t(TSig.mo) +TSig.mo %*% TSig.oo.inv %*% yhat.cent%*%t(Tmu[na.ind]) + TSig.mo %*% TSig.oo.inv %*% yhat.cent%*%t(yhat.cent)%*%TSig.oo.inv%*%t(TSig.mo) + Sig.mm.o
          ycm.hat = yc.hat %*% t(Tmu[na.ind]) + (yc.hat%*%t(yf)%*%TF + yc2.hat %*% TC - yc.hat%*%t(Tmu[-na.ind])) %*% TSig.oo.inv %*% t(TSig.mo)
          y2.hat[,,i] = t(TO)%*%t(TF)%*%yf %*% (t(yf)%*%TF%*%TO + t(yc.hat)%*%TC%*%TO + t(ym.hat)%*%TM) + t(TO)%*%t(TC) %*% (yc.hat%*%t(yf)%*%TF%*%TO + yc2.hat%*%TC%*%TO + ycm.hat %*% TM) + t(TM)%*%(ym.hat%*%t(yf)%*%TF%*%TO + t(ycm.hat)%*%TC%*%TO + ym2.hat%*%TM)
          #y2.hat[,,i] = t(TM)%*%ym2.hat%*%TM + t(TM)%*%ym.hat%*%t(yf)%*%TF%*%TO + t(TM)%*%t(ycm.hat)%*%TC%*%TO + t(TO)%*%t(TF)%*%yf%*%t(ym.hat)%*%TM + t(TO)%*%t(TF)%*%yf%*%t(yf)%*%TF%*%TO + t(TO)%*%t(TF)%*%yf%*%t(yc.hat)%*%TC%*%TO + t(TO)%*%t(TC)%*%ycm.hat%*%TM + t(TO)%*%t(TC)%*%yc.hat%*%t(yf)%*%TF%*%TO + t(TO)%*%t(TC)%*%yc2.hat%*%TC%*%TO
        }else{ #有遺失、沒設限(M+F)
          TSig.oo.inv = matrix(0, no, no)
          for(j in 1:n){
            if(j == 1){
              idx1 = 1: cumsum.o[1]
            }else{
              idx1 = (cumsum.o[j-1]+1): cumsum.o[j]
            }
            TSig.oo.inv[idx1,idx1] = solve(TSig.oo[idx1,idx1])
          }
          #TSig.oo.inv = solve(TSig.oo)
          yhat.cent = t(TF)%*%yf - Tmu[-na.ind]
          ym.hat = Tmu[na.ind] + TSig.mo %*% TSig.oo.inv %*% yhat.cent
          y.hat[,i] =  t(TO)%*%t(TF)%*%yf + t(TM)%*%ym.hat
          Sig.mm.o = TSig.mm - TSig.mo %*% TSig.oo.inv %*% t(TSig.mo)
          ym2.hat = Tmu[na.ind]%*%t(Tmu[na.ind]) + Tmu[na.ind]%*%t(yhat.cent)%*%TSig.oo.inv%*%t(TSig.mo) + Tmu[na.ind]%*%t(yc.hat)%*%TC%*%TSig.oo.inv%*%t(TSig.mo)+ TSig.mo %*% TSig.oo.inv %*% yhat.cent%*%t(Tmu[na.ind]) +TSig.mo %*% TSig.oo.inv %*% yhat.cent%*%t(yhat.cent)%*%TSig.oo.inv%*%t(TSig.mo) + TSig.mo %*% TSig.oo.inv %*% yhat.cent%*%t(yc.hat)%*%TC%*%TSig.oo.inv%*%t(TSig.mo)+ TSig.mo %*% TSig.oo.inv %*% t(TC)%*%yc.hat%*%t(Tmu[na.ind]) + TSig.mo %*% TSig.oo.inv %*% t(TC)%*%yc.hat%*%t(yhat.cent)%*%TSig.oo.inv%*%t(TSig.mo)+TSig.mo %*% TSig.oo.inv %*% t(TC)%*%yc2.hat%*%TC%*%TSig.oo.inv%*%t(TSig.mo)+ Sig.mm.o
          #ym2.hat = Tmu[na.ind]%*%t(Tmu[na.ind]) + Tmu[na.ind]%*%t(yhat.cent)%*%TSig.oo.inv%*%t(TSig.mo) + TSig.mo %*% TSig.oo.inv %*% yhat.cent%*%t(Tmu[na.ind]) + TSig.mo %*% TSig.oo.inv %*% yhat.cent%*%t(yhat.cent)%*%TSig.oo.inv%*%t(TSig.mo) + Sig.mm.o
          y2.hat[,,i] = t(TM)%*%ym2.hat%*%TM + t(TM)%*%ym.hat%*%t(yf)%*%TF%*%TO + t(TO)%*%t(TF)%*%yf%*%t(ym.hat)%*%TM + t(TO)%*%t(TF)%*%yf%*%t(yf)%*%TF%*%TO 
        }
      }else{ #沒遺失情況
        if(num.cen!=0){ #沒遺失、有設限(C,C+F)
          yf.cent = yf - Tmu[-cen.ind] 
          mu.cf = Tmu[cen.ind] + TSig.cf %*% TSig.ff.inv %*% yf.cent
          Sig.cc.f = TSig.cc - TSig.cf %*% TSig.ff.inv %*% t(TSig.cf)
          #evaluate right or left censored
          d0 = t(cen)[,cen.subj]
          d0 = d0[which(d0!=0)]
          bound = yc #censored value
          A=rep(-Inf,length(d0))
          B=rep(Inf,length(d0))
          for(t in 1:length(d0)){
            if(d0[t]==1){
              B[t]=bound[t]
            }else{
              A[t]=bound[t]
            }
          }
          ycj.hat = meanvarTMD(lower = A[1:cumsum.c[1]],upper = B[1:cumsum.c[1]] ,mu = mu.cf[1:cumsum.c[1]] ,Sigma = Sig.cc.f[1:cumsum.c[1],1:cumsum.c[1]],dist = 'normal')
          yc.hat[1:cumsum.c[1]] = ycj.hat$mean
          yc2.hat[1:cumsum.c[1], 1:cumsum.c[1]] = ycj.hat$EYY
          if(nc != 1){
            for(j in 2: nc){
              ycj.hat = meanvarTMD(lower=A[(cumsum.c[j-1]+1):cumsum.c[j]], upper=B[(cumsum.c[j-1]+1):cumsum.c[j]] ,mu=c(mu.cf[(cumsum.c[j-1]+1):cumsum.c[j]]) ,Sigma=Sig.cc.f[(cumsum.c[j-1]+1):cumsum.c[j],(cumsum.c[j-1]+1):cumsum.c[j]],dist = 'normal')
              yc.hat[(cumsum.c[j-1]+1):cumsum.c[j]] = ycj.hat$mean
              yc2.hat[(cumsum.c[j-1]+1):cumsum.c[j], (cumsum.c[j-1]+1):cumsum.c[j]] = ycj.hat$EYY
            }}
          #TSig.oo.inv = solve(TSig.oo)
          yhat.cent = t(TF)%*%yf + t(TC)%*%yc.hat - Tmu
          y.hat[,i] =  t(TO)%*%t(TF)%*%yf + t(TO)%*%t(TC)%*%yc.hat
          y2.hat[,,i] = t(TO)%*%t(TF)%*%yf %*% (t(yf)%*%TF%*%TO + t(yc.hat)%*%TC%*%TO) + t(TO)%*%t(TC) %*% (yc.hat%*%t(yf)%*%TF%*%TO + yc2.hat%*%TC%*%TO)
          #y2.hat = t(TO)%*%t(TF)%*%yf%*%t(yf)%*%TF%*%TO + t(TO)%*%t(TF)%*%yf%*%t(yc.hat)%*%TC%*%TO + t(TO)%*%t(TC)%*%yc.hat%*%t(yf)%*%TF%*%TO + t(TO)%*%t(TC)%*%yc2.hat%*%TC%*%TO
        }else{ #沒遺失、沒設限 (F)
          #TSig.oo.inv = solve(TSig.oo)
          yhat.cent = t(TF)%*%yf - Tmu
          y.hat[,i] = t(TO)%*%t(TF)%*%yf 
          y2.hat[,,i] = t(TO)%*%t(TF)%*%yf%*%t(yf)%*%TF%*%TO 
        }
      }
      
      Omega = rep(Z[,i],each=p) *(y2.hat[,,i] - y.hat[,i] %*% t(Tmu) - Tmu %*% t(y.hat[,i]) + Tmu %*% t(Tmu))
      
      # CM-step
      mu[,i] = colSums(Z[,i] * matrix(y.hat[,i] ,ncol = p, byrow =T))/ni[i] #更新mu值
      sum.Ome[,,i] = Omega[1:p, 1:p]
      for(j in 2: n) sum.Ome[,,i] = sum.Ome[,,i] + Omega[(cumsum.p[j-1]+1):cumsum.p[j], (cumsum.p[j-1]+1):cumsum.p[j]]
      #Sigma[,,i] = sum.Ome/ ni[i] #更新Sigma值
      L[,,i]=eigen(sum.Ome[,,i])$vectors
      omega[,,i]=diag(eigen(sum.Ome[,,i])$values)
    } 
    Aj = apply(omega, c(1,2), sum)/ det(apply(omega, c(1,2), sum))^(1/p) #Ahat 
    lambda=det(apply(omega, c(1,2), sum))^(1/p)/n #Lamdahat    
    
    for (i in 1:g) {
      Sigma[,,i]=lambda * (L[,,i] %*% Aj %*% t(L[,,i]))
    }
    # new log-likelihood
    w.den = log.like_B(Ycm, cen, cen1, w, mu, Sigma, yc, Oj, Cj, Fj, cen.subj, cen.posi)
    
    max.wden = apply(w.den, 1, max)
    w.den = exp(w.den - max.wden)
    indv.den = rowSums(w.den)
    log.indv.den = log(indv.den) + max.wden
    new.loglik = sum(log.indv.den)
    
    diff.lnL = (new.loglik-old.loglik)
    #if(iter%%per == 0) cat('iter = ', iter, ',\t obs.loglik = ', new.loglik, ',\t diff = ', diff.lnL, sep = ' ', '\n')
    if(diff.lnL < tol || iter >= max.iter) break
    old.loglik = new.loglik 
  }
  #cat('iter = ', iter, ',\t obs.loglik = ', new.loglik, sep = '', '\n')
  end = proc.time()[1]
  run.sec = end - begin #總執行時間
  #cat('It took', run.sec, 'seconds.\n')
  #cat(paste(rep('-', 50), sep = '', collapse = ''), '\n')
  #cat('mu =', mu, '\n')
  #cat('Sigma =', '\n')  
  #print(Sigma)
  para = list(mu=mu, Sigma = Sigma)
  est = c(mu, Sigma[vechS])
  #===========================#
  Z = w.den / indv.den 
  post.clus = matrix(apply(Z, 1, order),nrow=g)[g,] 
  ### predicting missing value(?w?????ƿ򥢭?) ###
  Y.pred = array(NA, dim=c(n,p,g))
  for(i in 1:g){
    Y.pred[,,i]=t(matrix(y.hat[,i],p,n))
  }
  post.pred = apply(rep(Z, p) * aperm(Y.pred, perm = c(1,3,2)), 3, rowSums)
  #=========================#
  alpha=g*p+g-1; beta=p*(p+1)/2; m=alpha+g*beta -(g-1)*p    #總參數個數
  aic = 2 * m - 2 * new.loglik
  bic = m * log(n) - 2 * new.loglik
  icl = bic - 2*sum(Z*log(Z+1e-300))
  awe = 2*(3/2+log(n))- 2 * new.loglik
  model.inf = c(run.sec = run.sec, iter = iter,loglik = new.loglik, aic = aic, bic = bic, icl=icl, awe=awe, m=m)
  data.inf = list(cen.subj = cen.subj, na.subj=na.subj, cen.rate=num.cen/np, na.rate=num.na/np)
  print(round(rbind(model.inf),4))
  return(list(run.sec = run.sec, iter = iter, model.inf = model.inf, para = para, est=est, data.inf=data.inf, post.pred = post.pred, post.clus = post.clus))    
}

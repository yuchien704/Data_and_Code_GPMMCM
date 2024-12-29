vech.posi=function(dim) cbind(rep(1:dim, 1:dim), unlist(mapply(':', 1, 1:dim)))
log.like.B = function(Ycm, cen, na, w, mu, Sigma, Oj, Cj, Fj)
{
  p=ncol(Ycm)
  n=nrow(Ycm)
  g=ncol(mu)
  pjo = p - rowSums(na) #每一列看得見的維度
  cen1 =ifelse(cen,1,0)
  pjf = pjo - rowSums(cen1, na.rm = T)
  na.subj = which(rowSums(na)!=0) #只有遺失列
  cen.subj = which(apply(cen, 1, sum, na.rm=T)!=0) #只有設限列
  num.na = length(na.subj) #遺失值個數
  num.cen = length(cen.subj) #設限值總個數
  log.cdf = matrix(0, n, g)
  w.den = delta.f = log.det.Sig.inv = matrix(NA, n, g)
  for(i in 1:g){
    if(num.cen!=0){
      for(j in c(cen.subj)){
        if(pjf[j] == pjo[j]){
          yf.cent = Oj[[j]]%*%Ycm[j,] - Oj[[j]]%*%mu[,i]
          log.det.Sig.inv[j,i] = log(det(solve(Oj[[j]]%*%Sigma[,,i]%*%t(Oj[[j]]))))  # log of det(TSig.ff.inv)
          delta.f[j,i] = t(yf.cent) %*% solve(Fj[[j]]%*%Oj[[j]]%*%Sigma[,,i]%*%t(Oj[[j]])%*%t(Fj[[j]])) %*% yf.cent
        }else{
          if(pjf[j] == 0){
            yf.cent = 0  
            log.det.Sig.inv[j,i] = 0
            delta.f[j,i] = 0
            mu.cf = Oj[[j]]%*%mu[,i] 
            Sig.cc.f = Oj[[j]]%*%Sigma[,,i]%*%t(Oj[[j]])
            if(sum(cen[j,]!=0,na.rm = T) == 1){
              h = cen[j,] # 第j種設限情況
              if(h[which(h!=0)]==1){ #左設限
                log.cdf[j,i] = log(pnorm(Cj[[j]]%*%Oj[[j]]%*%Ycm[j,], mean=mu.cf, sd=sqrt(Sig.cc.f)))
              }else log.cdf[j,i] = log(1-pnorm(Cj[[j]]%*%Oj[[j]]%*%Ycm[j,], mean=mu.cf, sd=sqrt(Sig.cc.f))) #右設限
            }else{
              d0 = cen[j,]
              d0 = d0[which(d0!=0)]
              bound = Cj[[j]]%*%Oj[[j]]%*%Ycm[j,] #censored value
              A=rep(-Inf,length(d0))
              B=rep(Inf,length(d0))
              for(t in 1:length(d0)){
                if(d0[t]==1){
                  B[t]=bound[t]
                }else{
                  A[t]=bound[t]
                }
              } #censored bound
              log.cdf[j,i] = log(pmvnorm(lower=A, upper=B, mean=c(mu.cf), sigma=Sig.cc.f))[1]
            }
          }else{
            yf.cent = Fj[[j]]%*%Oj[[j]]%*%Ycm[j,] - Fj[[j]]%*%Oj[[j]]%*%mu[,i]
            mu.cf = Cj[[j]]%*%Oj[[j]]%*%mu[,i] + Cj[[j]]%*%Oj[[j]]%*%Sigma[,,i]%*%t(Oj[[j]])%*%t(Fj[[j]]) %*% solve(Fj[[j]]%*%Oj[[j]]%*%Sigma[,,i]%*%t(Oj[[j]])%*%t(Fj[[j]]))%*% yf.cent
            Sig.cc.f = Cj[[j]]%*%Oj[[j]]%*%Sigma[,,i]%*%t(Oj[[j]])%*%t(Cj[[j]]) - Cj[[j]]%*%Oj[[j]]%*%Sigma[,,i]%*%t(Oj[[j]])%*%t(Fj[[j]]) %*% solve(Fj[[j]]%*%Oj[[j]]%*%Sigma[,,i]%*%t(Oj[[j]])%*%t(Fj[[j]])) %*% Fj[[j]]%*%Oj[[j]]%*%Sigma[,,i]%*%t(Oj[[j]])%*%t(Cj[[j]])
            if(sum(cen[j,]!=0,na.rm = T) == 1){
              h = cen[j,] # 第j種設限情況
              if(h[which(h!=0)]==1){ #左設限
                log.cdf[j,i] = log(pnorm(Cj[[j]]%*%Oj[[j]]%*%Ycm[j,], mean=mu.cf, sd=sqrt(Sig.cc.f)))
              }else log.cdf[j,i] = log(1-pnorm(Cj[[j]]%*%Oj[[j]]%*%Ycm[j,], mean=mu.cf, sd=sqrt(Sig.cc.f))) #右設限
            }else{
              d0 = cen[j,]
              d0 = d0[which(d0!=0)]
              bound = Cj[[j]]%*%Oj[[j]]%*%Ycm[j,] #censored value
              A=rep(-Inf,length(d0))
              B=rep(Inf,length(d0))
              for(t in 1:length(d0)){
                if(d0[t]==1){
                  B[t]=bound[t]
                }else{
                  A[t]=bound[t]
                }
              } #censored bound
              log.cdf[j,i] = log(pmvnorm(lower=A, upper=B, mean=c(mu.cf), sigma=Sig.cc.f))[1]
            }
            log.det.Sig.inv[j,i] = log(det(solve(Fj[[j]]%*%Oj[[j]]%*%Sigma[,,i]%*%t(Oj[[j]])%*%t(Fj[[j]]))))  # log of det(TSig.ff.inv)
            delta.f[j,i] = t(yf.cent) %*% solve(Fj[[j]]%*%Oj[[j]]%*%Sigma[,,i]%*%t(Oj[[j]])%*%t(Fj[[j]])) %*% yf.cent
          }
        }
      }
      for(j in (1:n)[-cen.subj]){
        yf.cent = Oj[[j]]%*%Ycm[j,]-Oj[[j]]%*%mu[,i] ; log.cdf[j,i] = 0 
        log.det.Sig.inv[j,i] = log(det(solve(Oj[[j]]%*%Sigma[,,i]%*%t(Oj[[j]]))))  # log of det(TSig.ff.inv)
        delta.f[j,i] = t(yf.cent) %*% solve(Oj[[j]]%*%Sigma[,,i]%*%t(Oj[[j]])) %*% yf.cent
      }
    }else {
      for(j in 1:n){
        yf.cent = Oj[[j]]%*%Ycm[j,] - Oj[[j]]%*%mu[,i] ; log.cdf[j,i] = 0 
        log.det.Sig.inv[j,i] = log(det(solve(Oj[[j]]%*%Sigma[,,i]%*%t(Oj[[j]]))))  # log of det(TSig.ff.inv)
        delta.f[j,i] = t(yf.cent) %*% solve(Oj[[j]]%*%Sigma[,,i]%*%t(Oj[[j]])) %*% yf.cent
      }
    }
    w.den[,i] = log(w[i]) + log.cdf[,i] -0.5*log(2*pi)*pjf + 0.5*log.det.Sig.inv[,i] - 0.5*delta.f[,i]
  }
  return(w.den)
}

VVI.GMIXCMB.ECM.new = function(Ycm, Y.knn, cen, g, true.clus, tol = 1e-5, neq.S = T, max.iter, per)
{
  begin = proc.time()[1]
  iter =0
  # Basic setting
  n = nrow(Ycm) #樣本數
  p = ncol(Ycm) #參數個數
  # Specify initial values
  mu = matrix(NA, p, g)
  Sigma = array(NA, dim = c(p,p,g))
  lambda= rep(NA, g) 
  W=array(NA, dim = c(p,p,g))
  Bj = array(NA, dim = c(p,p,g))
if(g==1) class=rep(1,n) else class=kmeans(Y.knn, g)$cluster
   w = table(class)/nrow(Y.knn)
  for(i in 1:g)
  {
    mu[,i] = colMeans(Y.knn[class == i,]) 
    Sigma[,,i] = diag(diag(cov(Y.knn)))
  }
  Ip = diag(p)  #pxp的單位矩陣
  cen1 =ifelse(cen,1,0)
  na = is.na(cen1) #判斷cen是否遺失
  pjo = p - rowSums(na) #每一列看得見的維度
  pjf = pjo - rowSums(cen1, na.rm = T)
  na.subj = which(rowSums(na)!=0) #只有遺失列
  cen.subj = which(apply(cen, 1, sum, na.rm=T)!=0) #只有設限列
  num.na = length(na.subj)#遺失值個數
  num.cen = length(cen.subj)#設限值總個數
  na.posi = cen.posi = Oj = Mj = Fj = Cj = as.list(n)
  for(j in 1:n)
  {
    Oj[[j]] = Ip
    Mj[[j]] = NULL
    Ipjo = diag(pjo[j])
    Fj[[j]] = Ipjo
    Cj[[j]] = NULL
    if(sum(is.na(cen[j,]))!=0){
      na.posi[[j]] = which(is.na(cen[j,]))
      Oj[[j]] = matrix(Ip[-na.posi[[j]],], ncol = p)
      Mj[[j]] = matrix(Ip[na.posi[[j]],], ncol = p)
    }
    if(sum(cen[j, ], na.rm=T) != 0){
      cen.posi[[j]] = which(na.omit(cen[j, ]) != 0)
      Fj[[j]] = matrix(Ipjo[-cen.posi[[j]],], ncol = pjo[j])
      Cj[[j]] = matrix(Ipjo[cen.posi[[j]], ], ncol = pjo[j])}
  }
  Ycm[is.na(Ycm)] = 9999
  Ycm = as.matrix(Ycm) 
  vechS = vech.posi(p) #pxp矩陣下三角位置
  
  # observed data log-likelihood
  w.den = log.like.B(Ycm, cen, na, w, mu, Sigma, Oj, Cj, Fj)
  max.wden = apply(w.den, 1, max)
  w.den = exp(w.den - max.wden)
  indv.den = rowSums(w.den)
  log.indv.den = log(indv.den) + max.wden
  old.loglik = sum(log.indv.den)  
  cat(paste(rep('-', 20), sep = '', collapse = ''), 'VVI.GMIXCMB: ', 'g =', g, paste(rep('-', 20), sep = '', collapse = ''), '\n')
  cat("iter =", iter, "\t loglik =", old.loglik, "\n") 
  repeat
  {
    iter = iter+1
    # E-step
    Z = w.den / indv.den
    ni = colSums(Z) 
    w = ni / n
    p2=p^2
    y.hat = array(NA, dim=c(n, p, g))
    y2.hat = array(0, dim=c(p, p, n, g))
    for(i in 1:g)
    {
      for (j in 1:n){
        if(pjo[j]!=p){# 遺失
          mu.o = Oj[[j]]%*%mu[,i]
          mu.m = Mj[[j]]%*%mu[,i]
          Sig.oo = Oj[[j]]%*%Sigma[,,i]%*%t(Oj[[j]])
          Sig.om = Oj[[j]]%*%Sigma[,,i]%*%t(Mj[[j]])
          Sig.mm = Mj[[j]]%*%Sigma[,,i]%*%t(Mj[[j]])
          if(pjf[j]==pjo[j]){ #F+M
            yf = Oj[[j]]%*%Ycm[j,]
            yhat.cent = yf - mu.o
            ym.hat = mu.m + t(Sig.om) %*% solve(Sig.oo) %*% yhat.cent
            y.hat[j,,i] = t(t(Oj[[j]])%*%yf + t(Mj[[j]])%*%ym.hat)
            Sig.mm.o = Sig.mm - t(Sig.om) %*% solve(Sig.oo) %*% Sig.om
            ym2.hat = mu.m%*%t(mu.m) + mu.m%*%t(yhat.cent)%*%solve(Sig.oo)%*%Sig.om +t(Sig.om)%*% solve(Sig.oo) %*% yhat.cent %*%t(mu.m) +t(Sig.om)%*% solve(Sig.oo) %*% yhat.cent%*%t(yhat.cent)%*%solve(Sig.oo)%*%Sig.om + Sig.mm.o
            y2.hat[,,j,i] = t(Oj[[j]])%*%yf %*% (t(yf)%*%Oj[[j]] + t(ym.hat)%*%Mj[[j]]) + t(Mj[[j]])%*%(ym.hat%*%t(yf)%*%Oj[[j]] + ym2.hat%*%Mj[[j]])
          }else if(pjf[j]==0){ #C+M
            mu.cf = mu.o 
            Sig.cc.f = Sig.oo
            d0 = cen[j,]
            d0 = d0[which(d0!=0)]
            bound = Oj[[j]]%*%Ycm[j,] #censored value
            A=rep(-Inf,length(d0))
            B=rep(Inf,length(d0))
            for(t in 1:length(d0)){
              if(d0[t]==1){
                B[t]=bound[t]
              }else{
                A[t]=bound[t]
              }
            } #censored bound
            EX=meanvarTMD(lower=A,upper=B,mu=mu.cf,Sigma=Sig.cc.f,dist="normal")
            yc.hat = EX$mean
            yc2.hat = EX$EYY
            yhat.cent = yc.hat - mu.o
            ym.hat = mu.m + t(Sig.om) %*% solve(Sig.oo) %*% yhat.cent
            y.hat[j,,i] = t(t(Oj[[j]])%*%t(Cj[[j]])%*%yc.hat + t(Mj[[j]])%*%ym.hat)
            Sig.mm.o = Sig.mm - t(Sig.om) %*% solve(Sig.oo) %*% Sig.om
            ym2.hat = mu.m%*%t(mu.m) + mu.m%*%t(yhat.cent)%*%solve(Sig.oo)%*%Sig.om +t(Sig.om)%*% solve(Sig.oo) %*% yhat.cent %*%t(mu.m) +t(Sig.om)%*% solve(Sig.oo) %*% (yc2.hat-yc.hat%*%t(mu.o)-mu.o%*%t(yc.hat)+mu.o%*%t(mu.o))%*%solve(Sig.oo)%*%Sig.om + Sig.mm.o
            ycm.hat = yc.hat %*% t(mu.m) + (yc2.hat - yc.hat%*%t(mu.o)) %*% solve(Sig.oo) %*% Sig.om
            y2.hat[,,j,i] = t(Oj[[j]])%*% (yc2.hat%*%Oj[[j]] + ycm.hat %*% Mj[[j]]) + t(Mj[[j]])%*%(t(ycm.hat)%*%Oj[[j]] + ym2.hat%*%Mj[[j]])
          }else{ #F+C+M
            Sig.ff = Fj[[j]]%*%Sig.oo%*%t(Fj[[j]])
            Sig.cf = Cj[[j]]%*%Sig.oo%*%t(Fj[[j]])
            Sig.cc = Cj[[j]]%*%Sig.oo%*%t(Cj[[j]])
            yf = Fj[[j]]%*%Oj[[j]]%*%Ycm[j,]
            yf.cent = yf - Fj[[j]]%*%mu.o
            mu.cf = Cj[[j]]%*%mu.o + Sig.cf %*% solve(Sig.ff)%*% yf.cent
            Sig.cc.f = Sig.cc - Sig.cf %*% solve(Sig.ff) %*% t(Sig.cf)
            d0 = cen[j,]
            d0 = d0[which(d0!=0)]
            bound = Cj[[j]]%*%Oj[[j]]%*%Ycm[j,] #censored value
            A=rep(-Inf,length(d0))
            B=rep(Inf,length(d0))
            for(t in 1:length(d0)){
              if(d0[t]==1){
                B[t]=bound[t]
              }else{
                A[t]=bound[t]
              }
            } #censored bound
            EX=meanvarTMD(lower=A,upper=B,mu=mu.cf,Sigma=Sig.cc.f,dist="normal")
            yc.hat = EX$mean
            yc2.hat = EX$EYY
            yhat.cent = t(Fj[[j]])%*%yf+ t(Cj[[j]])%*%yc.hat - mu.o
            ym.hat = mu.m + t(Sig.om) %*% solve(Sig.oo) %*% yhat.cent
            y.hat[j,,i] = t(t(Oj[[j]])%*%t(Fj[[j]])%*%yf + t(Oj[[j]])%*%t(Cj[[j]])%*%yc.hat + t(Mj[[j]])%*%ym.hat)
            Sig.mm.o = Sig.mm - t(Sig.om) %*% solve(Sig.oo) %*% Sig.om
            ym2.hat = mu.m%*%t(ym.hat) + t(Sig.om) %*% solve(Sig.oo)%*%(t(Fj[[j]])%*%yf- mu.o)%*%t(ym.hat) + t(Sig.om)%*%solve(Sig.oo)%*%t(Cj[[j]])%*%yc.hat%*%t(mu.m) + t(Sig.om)%*%solve(Sig.oo)%*%t(Cj[[j]])%*%yc.hat%*%t(t(Fj[[j]])%*%yf - mu.o)%*%solve(Sig.oo)%*%Sig.om + t(Sig.om)%*%solve(Sig.oo)%*%t(Cj[[j]])%*%yc2.hat%*%Cj[[j]]%*%solve(Sig.oo)%*%Sig.om + Sig.mm.o
            ycm.hat = yc.hat %*% t(mu.m) + (yc.hat%*%t(yf)%*%Fj[[j]] + yc2.hat %*% Cj[[j]] - yc.hat%*%t(mu.o)) %*% solve(Sig.oo) %*% Sig.om
            y2.hat[,,j,i] = t(Oj[[j]])%*%t(Fj[[j]])%*%yf %*% (t(yf)%*%Fj[[j]]%*%Oj[[j]] + t(yc.hat)%*%Cj[[j]]%*%Oj[[j]] + t(ym.hat)%*%Mj[[j]]) + t(Oj[[j]])%*%t(Cj[[j]]) %*% (yc.hat%*%t(yf)%*%Fj[[j]]%*%Oj[[j]] + yc2.hat%*%Cj[[j]]%*%Oj[[j]] + ycm.hat %*% Mj[[j]]) + t(Mj[[j]])%*%(ym.hat%*%t(yf)%*%Fj[[j]]%*%Oj[[j]] + t(ycm.hat)%*%Cj[[j]]%*%Oj[[j]] + ym2.hat%*%Mj[[j]])
          }
        }else{ # 沒有遺失
          mu.o = mu[,i]
          Sig.oo = Sigma[,,i]
          if(pjf[j]==pjo[j]){ #F
            y.hat[j,,i] = Ycm[j,]
            y2.hat[,,j,i] = Ycm[j,]%*%t(Ycm[j,])
          }else if(pjf[j]==0){ #C
            mu.cf = mu.o
            Sig.cc.f = Sig.oo
            d0 = cen[j,]
            d0 = d0[which(d0!=0)]
            bound = Ycm[j,] #censored value
            A=rep(-Inf,length(d0))
            B=rep(Inf,length(d0))
            for(t in 1:length(d0)){
              if(d0[t]==1){
                B[t]=bound[t]
              }else{
                A[t]=bound[t]
              }
            } #censored bound
            EX=meanvarTMD(lower=A,upper=B,mu=mu.cf,Sigma=Sig.cc.f,dist="normal")
            yc.hat = EX$mean
            yc2.hat = EX$EYY
            yhat.cent = t(Cj[[j]])%*%yc.hat - mu.o
            y.hat[j,,i] = t(t(Cj[[j]])%*%yc.hat)
            y2.hat[,,j,i] = t(Cj[[j]]) %*% (yc2.hat%*%Cj[[j]])
          }else{ #F+C
            Sig.ff = Fj[[j]]%*%Sig.oo%*%t(Fj[[j]])
            Sig.cf = Cj[[j]]%*%Sig.oo%*%t(Fj[[j]])
            Sig.cc = Cj[[j]]%*%Sig.oo%*%t(Cj[[j]])
            yf = Fj[[j]]%*%Ycm[j,]
            yf.cent = yf - Fj[[j]]%*%mu.o
            mu.cf = Cj[[j]]%*%mu.o + Sig.cf %*% solve(Sig.ff)%*% yf.cent
            Sig.cc.f = Sig.cc - Sig.cf %*% solve(Sig.ff) %*% t(Sig.cf)
            d0 = cen[j,]
            d0 = d0[which(d0!=0)]
            bound = Cj[[j]]%*%Ycm[j,] #censored value
            A=rep(-Inf,length(d0))
            B=rep(Inf,length(d0))
            for(t in 1:length(d0)){
              if(d0[t]==1){
                B[t]=bound[t]
              }else{
                A[t]=bound[t]
              }
            } #censored bound
            EX=meanvarTMD(lower=A,upper=B,mu=mu.cf,Sigma=Sig.cc.f,dist="normal")
            yc.hat = EX$mean
            yc2.hat = EX$EYY
            yhat.cent = t(Fj[[j]])%*%yf+ t(Cj[[j]])%*%yc.hat - mu.o
            y.hat[j,,i] = t(t(Fj[[j]])%*%yf + t(Cj[[j]])%*%yc.hat)
            y2.hat[,,j,i] = t(Fj[[j]])%*%yf %*% (t(yf)%*%Fj[[j]] + t(yc.hat)%*%Cj[[j]]) + t(Cj[[j]]) %*% (yc.hat%*%t(yf)%*%Fj[[j]] + yc2.hat%*%Cj[[j]])
          }
        }
      }
      # CM-step
      mu[,i] = colSums(Z[,i] * y.hat[,,i])/ni[i] #更新mu值
      W[,,i] = (apply(rep(Z[,i], each=p2)*y2.hat[,,,i], 1:2, sum) - colSums(Z[,i]*y.hat[,,i])%*%t(mu[,i]) - 
                  mu[,i]%*%t(colSums(Z[,i]*y.hat[,,i])) + sum(Z[,i])*mu[,i]%*%t(mu[,i]))
      Bj[,,i] = diag(diag(W[,,i]))/det(diag(diag(W[,,i])))^(1/p)
      lambda[i] = det(diag(diag(W[,,i])))^(1/p)/ni[i]
      Sigma[,,i] = lambda[i]*Bj[,,i]
    } 
    
    # new log-likelihood
    w.den = log.like.B(Ycm, cen, na, w, mu, Sigma, Oj, Cj, Fj)
    max.wden = apply(w.den, 1, max)
    w.den = exp(w.den - max.wden)
    indv.den = rowSums(w.den)
    log.indv.den = log(indv.den) + max.wden
    new.loglik = sum(log.indv.den)
    diff.lnL = (new.loglik-old.loglik)
    if(iter%%per == 0) cat('iter = ', iter, ',\t obs.loglik = ', new.loglik, ',\t diff = ', diff.lnL, sep = ' ', '\n')
    if(diff.lnL < tol || iter >= max.iter) break
    old.loglik = new.loglik 
  }
  end = proc.time()[1]
  run.sec = end - begin #總執行時間
  para = list(w=w,mu=mu, Sigma = Sigma)
  est = c(mu, Sigma[vechS])
  Z = w.den / indv.den 
  post.clus = matrix(apply(Z, 1, order),nrow=g)[g,]
  if(length(unique(post.clus))==length(unique(true.clus)))
   {
    CCR=1-classError(true.clus,post.clus)$errorRate
   }
   else {
    CCR=sum(apply(table(post.clus, true.clus),1,max))/n
  }
  ARI=adjustedRandIndex(true.clus, post.clus)
  post.pred = apply(rep(Z, p) * aperm(y.hat, perm = c(1,3,2)), 3, rowSums)
  alpha=g*p+g-1;m=alpha+g*p   #總參數個數
  aic = 2 * m - 2 * new.loglik
  bic = m * log(n) - 2 * new.loglik
  icl = bic - 2*sum(Z*log(Z+1e-300))
  awe = icl + 3*m + m*log(n)  
  model.inf = c(run.sec = run.sec, iter = iter,loglik = new.loglik, aic = aic, bic = bic, icl=icl, awe=awe, ARI=ARI, CCR=CCR, m=m)
  data.inf = list(cen.subj = cen.subj, na.subj=na.subj, cen.rate=num.cen/(n*p), na.rate=num.na/(n*p))
  print(round(rbind(model.inf),4))
  return(list(run.sec = run.sec, iter = iter, model.inf = model.inf, para = para, est=est, data.inf=data.inf, post.pred = post.pred, post.clus = post.clus, ARI=ARI, CCR=CCR))    
}

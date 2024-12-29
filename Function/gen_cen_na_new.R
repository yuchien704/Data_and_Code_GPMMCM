gener.cen.na=function(Data, cen.rate, cen.type, na.rate)
{
  if(length(cen.rate)!=length(cen.type))
    stop('Waring: there is an input error!')
  cen.rate[which(cen.type==0)]=0
  p= ncol(Data)
  n=nrow(Data) #sample size
  r = sum(cen.rate != 0) # number of censored
  idx = which(cen.type != 0) #site of censored
  cen = matrix(0,n,p)
  if(r ==0) { Data;cen
  } else{
    for (k in 1:r) {
      if(cen.type[idx[k]]==1){
        ubd.L = quantile(Data[,idx[k]], prob=cen.rate[idx[k]], na.rm = T) 
        Data[which(Data[,idx[k]]<=ubd.L),idx[k]] = ubd.L
        cen[which(Data[,idx[k]]<=ubd.L),idx[k]]=1
      }else{
        ubd.R = quantile(Data[,idx[k]], prob=1-cen.rate[idx[k]], na.rm = T) 
        Data[which(Data[,idx[k]]>=ubd.R),idx[k]] = ubd.R
        cen[which(Data[,idx[k]]>=ubd.R),idx[k]]=2
      }
    }
  }
  if(na.rate==0){Data;cen
  }else{
    cen.posi = which(t(cen)!=0) # 觀測資料向量的設限位置
    keep.part = sample(p, n, replace = T) + p * 0:(n-1)#每列保留至少一個觀測值
    keep.part1 = unique(c(keep.part,cen.posi))
    num.na = floor(n*p*na.rate)
    na.posi = tabulate(sample((1:(n*p))[-keep.part1], num.na), nbins = n*p)
    na.posi = matrix(na.posi, ncol = p, byrow = T)
    Data[na.posi == 1] = NA
    cen[na.posi == 1] = NA
  }
  return(list(Data = Data, cen = cen))
}

library(mvtnorm)
library(tmvtnorm)
library(mclust)
library(MomTrunc)
Yc = as.matrix(read.table(paste(SPATH, '/Data/source/moesm.txt', sep=''), head=F))
Ycen = read.table(paste(SPATH, '/Data/source/moesm.cen.txt', sep=''), head=F)
cen = cbind(as.numeric(Ycen$V1=='<0.1'), as.numeric(Ycen$V2=='<0.1'), as.numeric(Ycen$V3=='<1'),
            as.numeric(Ycen$V4=='<1' | Ycen$V4=='<0.5'), as.numeric(Ycen$V5=='<1' | Ycen$V5=='<0.5'))
cen.rate = round(colMeans(cen)*100, 2)
Yc_mean = colMeans(Yc)
Yc_sd = apply(Yc,2,sd)

#table2
cname2 = c("Censoring proportion", "Sample mean", "Sample deviation")
rname2 = c('Cu','Pb','Zn','Ca','Mg')
Table2 = cbind(cen.rate,round(Yc_mean,3),round(Yc_sd,3))
rownames(Table2) = rname2
colnames(Table2) = cname2
write.csv(Table2, paste(SPATH,'/Results/Table2.csv',sep=""), row.names = TRUE)

#-----------------------------------------------#
load(paste(SPATH, "/Data/vdeq.RData", sep = ""))

#table3
T13 = (as.matrix(T13)[c(1,4,2,3),])
T23 = (as.matrix(T23)[c(1,4,2,3),])
Table3 = cbind(T13, T23)
rownames(Table3)=paste("(VVE,4)-C" ,c(1:4), sep = "")
colnames(Table3)=c(paste("(VVE,3)-C" ,c(1:3), sep = ""), paste("(VEV,3)-C" ,c(1:3), sep = ""))
write.csv(Table3, paste(SPATH,'/Results/Table3.csv',sep=""), row.names = TRUE)


#table4

W = round(rbind(w0, w1, w2),3)
M1 = t(round(cbind(mu0[,1], mu1[,1], mu2[,1]),3))
M2 = t(round(cbind(mu0[,2], mu1[,2], mu2[,2]),3))
M3 = t(round(cbind(mu0[,3], mu1[,3], mu2[,3]),3))

S01 = round(S0[,,1][lower.tri(S0[,,1])==F],3)
S02 = round(S0[,,2][lower.tri(S0[,,2])==F],3)
S03 = round(S0[,,3][lower.tri(S0[,,3])==F],3)

S11 = round(S1[,,1][lower.tri(S1[,,1])==F],3)
S12 = round(S1[,,2][lower.tri(S1[,,2])==F],3)
S13 = round(S1[,,3][lower.tri(S1[,,3])==F],3)

S21 = round(S2[,,1][lower.tri(S2[,,1])==F],3)
S22 = round(S2[,,2][lower.tri(S2[,,2])==F],3)
S23 = round(S2[,,3][lower.tri(S2[,,3])==F],3)

Sigma1 = rbind(S01, S11, S21)
Sigma2 = rbind(S02, S12, S22)
Sigma3 = rbind(S03, S13, S23)
colnames(Sigma1)=colnames(Sigma2)=colnames(Sigma3)=c("sigma11","sigma21","sigma22","sigma31","sigma32", "sigma33","sigma41","sigma42","sigma43","sigma44","sigma51","sigma52","sigma53","sigma54", "sigma55")

Table4 = data.frame("pi"=W, "mu1"=M1, "mu2"=M2, "mu3"=M3, "S1"=Sigma1, "S2"=Sigma2, "S3"=Sigma3, row.names = c("(VVV,3)","(VVE,3)", "(VEV,3)"))
write.csv(Table4, paste(SPATH,'/Results/Table4.csv',sep=""), row.names = TRUE)



library(mvtnorm)
library(tmvtnorm)
library(mclust)
library(MomTrunc)

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
source(paste(SPATH, '/Function/Both_gmixcm/EEI_GMIXCMB_1.R', sep=''))
source(paste(SPATH, '/Function/Both_gmixcm/EEE_GMIXCMB_1.R', sep=''))

Yc = scale(as.matrix(read.table(paste(SPATH, '/Data/source/moesm.txt', sep=''), head=F)))
Ycen = read.table(paste(SPATH, '/Data/source/moesm.cen.txt', sep=''), head=F)
cen = cbind(as.numeric(Ycen$V1=='<0.1'), as.numeric(Ycen$V2=='<0.1'), as.numeric(Ycen$V3=='<1'),
            as.numeric(Ycen$V4=='<1' | Ycen$V4=='<0.5'), as.numeric(Ycen$V5=='<1' | Ycen$V5=='<0.5'))
cen.rate = round(colMeans(cen)*100, 2)

g=1; clus1=kmeans(Yc, g)$cluster
tol = 1e-5;neq.S = T;max.iter=1000;tol.FG=1e-3;max.iter.FG=10;per=1
fit.vvv1 = VVV.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus1, tol, neq.S = T, max.iter, per)
fit.eee1 = EEE.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus1, tol, neq.S = F, max.iter, per)
fit.eev1 = EEV.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus1, tol, neq.S = T, max.iter, per)
fit.eve1 = EVE.GMIX.MVNCM.B.ECM(Yc, cen, g, clus1, tol, neq.S = T, max.iter, tol.FG, max.iter.FG, per)
fit.evv1 = EVV.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus1, tol, neq.S = T, max.iter, per)
fit.vee1 = VEE.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus1, tol, neq.S = T, max.iter, per)
fit.vev1 = VEV.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus1, tol, neq.S = T, max.iter, per)
fit.vve1 = VVE.GMIX.MVNCM.B.ECM(Yc, cen, g, clus1, tol, neq.S = T, max.iter, tol.FG, max.iter.FG, per)
fit.vei1 = VEI.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus1, tol, neq.S = T, max.iter, per)
fit.eei1 = EEI.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus1, tol, neq.S = F, max.iter, per)
fit.eii1 = EII.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus1, tol, neq.S = F, max.iter, per)
fit.evi1 = EVI.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus1, tol, neq.S = T, max.iter, per)
fit.vii1 = VII.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus1, tol, neq.S = T, max.iter, per)
fit.vvi1 = VVI.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus1, tol, neq.S = T, max.iter, per)
a1=c(fit.vvv1$model.inf[[5]],fit.eee1$model.inf[[5]],fit.eev1$model.inf[[5]],fit.eve1$model.inf[[5]],fit.evv1$model.inf[[5]],
     fit.vee1$model.inf[[5]],fit.vev1$model.inf[[5]],fit.vve1$model.inf[[5]],fit.vei1$model.inf[[5]],fit.eei1$model.inf[[5]],
     fit.vvi1$model.inf[[5]],fit.evi1$model.inf[[5]],fit.vii1$model.inf[[5]],fit.eii1$model.inf[[5]])
names(a1) = paste(c("VVV", "EEE", "EEV","EVE","EVV","VEE","VEV","VVE","VEI","EEI","VVI","EVI","VII","EII"),1,sep = ",")


g=2; clus2=kmeans(Yc, g)$cluster
fit.vvv2 = VVV.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus2, tol, max.iter, per)
fit.eee2 = EEE.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus2, tol, neq.S = F, max.iter, per)
repeat{
  fit.eev2 = try(EEV.GMIX.MVNCM.B.ECM.renew(Yc, Yc, cen, g, tol, neq.S = T, max.iter, per), silent=T)
  if(class(fit.eev2) != "try-error") break
}
fit.eve2 = EVE.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus2, tol, neq.S = T, max.iter, tol.FG, max.iter.FG, per)
fit.evv2 = EVV.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus2, tol, neq.S = T, max.iter, per)
fit.vee2 = VEE.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus2, tol, neq.S = T, max.iter, per)
fit.vev2 = VEV.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus2, tol, neq.S = T, max.iter, per)
fit.vve2 = VVE.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus2, tol, neq.S = T, max.iter, tol.FG, max.iter.FG, per)
fit.vei2 = VEI.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus2, tol, neq.S = T, max.iter, per)
repeat{
  fit.eei2 = try(EEI.GMIX.MVNCM.B.ECM.renew(Yc, Yc, cen, g, tol, neq.S = F, max.iter, per), silent=T)
  if(class(fit.eei2) != "try-error") break
}
fit.eii2 = EII.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus2, tol, neq.S = F, max.iter, per)
fit.evi2 = EVI.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus2, tol, neq.S = T, max.iter, per)
fit.vii2 = VII.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus2, tol, neq.S = T, max.iter, per)
fit.vvi2 = VVI.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus2, tol, neq.S = T, max.iter, per)
a2=c(fit.vvv2$model.inf[[5]],fit.eee2$model.inf[[5]],fit.eev2$model.inf[[5]],fit.eve2$model.inf[[5]],fit.evv2$model.inf[[5]],
     fit.vee2$model.inf[[5]],fit.vev2$model.inf[[5]],fit.vve2$model.inf[[5]],fit.vei2$model.inf[[5]],fit.eei2$model.inf[[5]],
     fit.vvi2$model.inf[[5]],fit.evi2$model.inf[[5]],fit.vii2$model.inf[[5]],fit.eii2$model.inf[[5]])
names(a2) = paste(c("VVV", "EEE", "EEV","EVE","EVV","VEE","VEV","VVE","VEI","EEI","VVI","EVI","VII","EII"),2,sep = ",")


g=3; clus3=kmeans(Yc, g)$cluster
repeat{
  fit.vvv3 = try(VVV.GMIX.MVNCM.B.ECM.renew(Yc, Yc, cen, g, tol, neq.S = T, max.iter, per), silent=T)
  if(class(fit.vvv3) != "try-error") break
}
fit.eee3 = EEE.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus3, tol, neq.S = F, max.iter, per)
fit.eev3 = EEV.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus3, tol, neq.S = T, max.iter, per)
repeat{
  fit.eev3 = try(EEV.GMIX.MVNCM.B.ECM.renew(Yc, Yc, cen, g, tol, neq.S = T, max.iter, per), silent=T)
  if(class(fit.eev3) != "try-error") break
}
fit.eve3 = EVE.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus3, tol, neq.S = T, max.iter, tol.FG, max.iter.FG, per)
fit.evv3 = EVV.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus3, tol, neq.S = T, max.iter, per)
fit.vee3 = VEE.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus3, tol, neq.S = T, max.iter, per)
fit.vev3 = VEV.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus3, tol, neq.S = T, max.iter, per)
fit.vve3 = VVE.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus3, tol, neq.S = T, max.iter, tol.FG, max.iter.FG, per)
fit.vei3 = VEI.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus3, tol, neq.S = T, max.iter, per)
repeat{
  fit.eei3 = try(EEI.GMIX.MVNCM.B.ECM.renew(Yc, Yc, cen, g, tol, neq.S = F, max.iter, per), silent=T)
  if(class(fit.eei3) != "try-error") break
}
fit.eii3 = EII.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus3, tol, neq.S = F, max.iter, per)
fit.evi3 = EVI.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus3, tol, neq.S = T, max.iter, per)
fit.vii3 = VII.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus3, tol, neq.S = T, max.iter, per)
fit.vvi3 = VVI.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus3, tol, neq.S = T, max.iter, per)
a3=c(fit.vvv3$model.inf[[5]],fit.eee3$model.inf[[5]],fit.eev3$model.inf[[5]],fit.eve3$model.inf[[5]],fit.evv3$model.inf[[5]],
     fit.vee3$model.inf[[5]],fit.vev3$model.inf[[5]],fit.vve3$model.inf[[5]],fit.vei3$model.inf[[5]],fit.eei3$model.inf[[5]],
     fit.vvi3$model.inf[[5]],fit.evi3$model.inf[[5]],fit.vii3$model.inf[[5]],fit.eii3$model.inf[[5]])
names(a3) = paste(c("VVV", "EEE", "EEV","EVE","EVV","VEE","VEV","VVE","VEI","EEI","VVI","EVI","VII","EII"),3,sep = ",")


g=4; clus4=kmeans(Yc, g)$cluster
repeat{
  fit.vvv4 = try(VVV.GMIX.MVNCM.B.ECM.renew(Yc, Yc, cen, g, tol, neq.S = T, max.iter, per), silent=T)
  if(class(fit.vvv4) != "try-error") break
}
fit.eee4 = EEE.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus4, tol, neq.S = F, max.iter, per)
repeat{
  fit.eev4 = try(EEV.GMIX.MVNCM.B.ECM.renew(Yc, Yc, cen, g, tol, neq.S = T, max.iter, per), silent=T)
  if(class(fit.eev4) != "try-error") break
}
repeat{
  fit.eve4 = try(EVE.GMIX.MVNCM.B.ECM.renew(Yc, Yc, cen, g, tol, neq.S = T, max.iter, tol.FG, max.iter.FG, per), silent=T)
  if(class(fit.eve4) != "try-error") break
}
repeat{
  fit.evv4 = try(EVV.GMIX.MVNCM.B.ECM.renew(Yc, Yc, cen, g, tol, neq.S = T, max.iter, per), silent=T)
  if(class(fit.evv4) != "try-error") break
}
fit.vee4 = VEE.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus4, tol, neq.S = T, max.iter, per)
repeat{
  fit.vev4 = try(VEV.GMIX.MVNCM.B.ECM.renew(Yc, Yc, cen, g, tol, neq.S = T, max.iter, per), silent=T)
  if(class(fit.vev4) != "try-error") break
}
repeat{
  fit.vve4 = try(VVE.GMIX.MVNCM.B.ECM.renew(Yc, Yc, cen, g, tol, neq.S = T, max.iter, tol.FG, max.iter.FG, per), silent=T)
  if(class(fit.vve4) != "try-error") break
}
fit.vei4 = VEI.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus4, tol, neq.S = T, max.iter, per)
repeat{
  fit.eei4 = try(EEI.GMIX.MVNCM.B.ECM.renew(Yc, Yc, cen, g, tol, neq.S = F, max.iter, per), silent=T)
  if(class(fit.eei4) != "try-error") break
}
fit.eii4 = EII.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus4, tol, neq.S = F, max.iter, per)
fit.evi4 = EVI.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus4, tol, neq.S = T, max.iter, per)
fit.vii4 = VII.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus4, tol, neq.S = T, max.iter, per)
fit.vvi4 = VVI.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus4, tol, neq.S = T, max.iter, per)
a4=c(NA,fit.eee4$model.inf[[5]],NA,fit.eve4$model.inf[[5]],fit.evv4$model.inf[[5]],
     fit.vee4$model.inf[[5]],NA,fit.vve4$model.inf[[5]],fit.vei4$model.inf[[5]],NA,
     fit.vvi4$model.inf[[5]],fit.evi4$model.inf[[5]],fit.vii4$model.inf[[5]],fit.eii4$model.inf[[5]])
names(a4) = paste(c("VVV", "EEE", "EEV","EVE","EVV","VEE","VEV","VVE","VEI","EEI","VVI","EVI","VII","EII"),4,sep = ",")


g=5; clus5=kmeans(Yc, g)$cluster
repeat{
  fit.vvv5 = try(VVV.GMIX.MVNCM.B.ECM.renew(Yc, Yc, cen, g, tol, neq.S = T, max.iter, per), silent=T)
  if(class(fit.vvv5) != "try-error") break
}
repeat{
  fit.eee5 = try(EEE.GMIX.MVNCM.B.ECM.renew(Yc, Yc, cen, g, tol, neq.S = F, max.iter, per), silent=T)
  if(class(fit.eee5) != "try-error") break
}
repeat{
  fit.eev5 = try(EEV.GMIX.MVNCM.B.ECM.renew(Yc, Yc, cen, g, tol, neq.S = T, max.iter, per), silent=T)
  if(class(fit.eev5) != "try-error") break
}
repeat{
  fit.eve5 = try(EVE.GMIX.MVNCM.B.ECM.renew(Yc, Yc, cen, g, tol, neq.S = T, max.iter, tol.FG, max.iter.FG, per), silent=T)
  if(class(fit.eve5) != "try-error") break
}
repeat{
  fit.evv5 = try(EVV.GMIX.MVNCM.B.ECM.renew(Yc, Yc, cen, g, tol, neq.S = T, max.iter, per), silent=T)
  if(class(fit.evv5) != "try-error") break
}
fit.vee5 = VEE.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus5, tol, neq.S = T, max.iter, per)
repeat{
  fit.vev5 = try(VEV.GMIX.MVNCM.B.ECM.renew(Yc, Yc, cen, g, tol, neq.S = T, max.iter, per), silent=T)
  if(class(fit.vev5) != "try-error") break
}
repeat{
  fit.vve5 = try(VVE.GMIX.MVNCM.B.ECM.renew(Yc, Yc, cen, g, tol, neq.S = T, max.iter, tol.FG, max.iter.FG, per), silent=T)
  if(class(fit.vve5) != "try-error") break
}
fit.vei5 = VEI.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus5, tol, neq.S = T, max.iter, per)
repeat{
  fit.eei5 = try(EEI.GMIX.MVNCM.B.ECM.renew(Yc, Yc, cen, g, tol, neq.S = F, max.iter, per), silent=T)
  if(class(fit.eei5) != "try-error") break
}
fit.eii5 = EII.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus5, tol, neq.S = F, max.iter, per)
fit.evi5 = EVI.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus5, tol, neq.S = T, max.iter, per)
fit.vii5 = VII.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus5, tol, neq.S = T, max.iter, per)
fit.vvi5 = VVI.GMIX.MVNCM.B.ECM(Yc,  cen, g, clus5, tol, neq.S = T, max.iter, per)
a5=c(NA,NA,NA,fit.eve5$model.inf[[5]],NA,
     fit.vee5$model.inf[[5]],fit.vev5$model.inf[[5]],NA,fit.vei5$model.inf[[5]],NA,
     fit.vvi5$model.inf[[5]],fit.evi5$model.inf[[5]],fit.vii5$model.inf[[5]],fit.eii5$model.inf[[5]])
names(a5) = paste(c("VVV", "EEE", "EEV","EVE","EVV","VEE","VEV","VVE","VEI","EEI","VVI","EVI","VII","EII"),5,sep = ",")


BIC_VDEQ5 = rbind(a1,a2,a3,a4,a5)
model = c(a1,a2,a3,a4,a5)
best_model = sort(model)[1:5]
colnames(BIC_VDEQ5)=c("VVV", "EEE", "EEV","EVE","EVV","VEE","VEV","VVE","VEI","EEI","VVI","EVI","VII","EII")
rownames(BIC_VDEQ5)=c(1:5)

#-------------------- fit --------------------------#
ysel = c('Cu','Pb','Zn','Ca','Mg')
x1 = c(-2,-2,-1,-1.5,-2)
xli = c(4,4,6,2,4) #colnames x-axis
alim= c(0.8,1.4,1.3,0.8,0.8) #colnames y-axis
lim= c(1,1.8,1.6,1,1) # y-axis max
at = c(-0.2,-0.25,-0.1,-0.1,-0.1)# y-axis min
bx = c(0.3, 0.3, 0.3, 0.2, 0.2)
c = list()
c[[1]] = which(cen[,1]==1)
c[[2]] = which(cen[,2]==1)
c[[3]] = which(cen[,3]==1)
c[[4]] = which(cen[,4]==1)
c[[5]] = which(cen[,5]==1)
nc1 = length(c[[1]])
nc2 = length(c[[2]])
nc3 = length(c[[3]])
nc4 = length(c[[4]])
nc5 = length(c[[5]])
cutoff1 = unique(Yc[c[[1]], 1])
cutoff2 = unique(Yc[c[[2]], 2])
cutoff3 = unique(Yc[c[[3]], 3])
cutoff4 = unique(Yc[c[[4]], 4])
cutoff5 = unique(Yc[c[[5]], 5])
cutof = c(cutoff1,cutoff2,cutoff3,cutoff4[1],cutoff5[1])

fit.VDEQ1 = fit.vve3
p.clus1 = fit.VDEQ1$post.clus
y.pred1 = fit.VDEQ1$post.pred
w1 = fit.VDEQ1$para[[1]]
mu1 = fit.VDEQ1$para[[2]]
S1 = fit.VDEQ1$para[[3]]

fit.VDEQ2 = fit.vev3
p.clus2 = fit.VDEQ2$post.clus
y.pred2 = fit.VDEQ2$post.pred
w2 = fit.VDEQ2$para[[1]]
mu2 = fit.VDEQ2$para[[2]]
S2 = fit.VDEQ2$para[[3]]

fit.VDEQ3 = fit.vve4
p.clus3 = fit.VDEQ3$post.clus

fit.VDEQ0 = fit.vvv3
p.clus00 = fit.VDEQ0$post.clus
y.pred0 = fit.VDEQ0$post.pred
w00 = fit.VDEQ0$para[[1]]
mu00 = fit.VDEQ0$para[[2]]
S00 = fit.VDEQ0$para[[3]]

## --------------------- 1->1, 2->3, 3->2 --------------- ##
p.clus0=p.clus00
p.clus0[which(p.clus00==2)]=3
p.clus0[which(p.clus00==3)]=2
w0 = c(w00[1],w00[3],w00[2])
mu0 = cbind(mu00[,1],mu00[,3],mu00[,2])
S0 = array(cbind(S00[,,1],S00[,,3],S00[,,2]) ,dim=c(5,5,3))
#--------------------------------------------------------- #

T13 = table(p.clus3,p.clus1)
T23 = table(p.clus3,p.clus2)

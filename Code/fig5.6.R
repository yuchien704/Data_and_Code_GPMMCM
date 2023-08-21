library(mvtnorm)
library(tmvtnorm)
library(mclust)
library(MomTrunc)
load(paste(SPATH, "/Data/vdeq.RData", sep = ""))

#--------------------fig5------------------------------------------------#

Col=c("Dark Orchid4","Purple","Orange","Dark Red","Red","Magenta","Dark Green","Violet Red","Royal Blue","Dark Goldenrod","Blue","Cyan","Dark Violet","Green")
Pch=c(1:14)
DD = 0
postscript(paste(SPATH, '/Results/fig5.eps', sep=''), width=6, height=30)
par(mar=c(4,4,2,0.5), oma = c(0,0,4,0),pty='s')
#VVV
plot(BIC_VDEQ5[,1],type='b',ylim=c(min(BIC_VDEQ5[!is.na(BIC_VDEQ5)])+DD ,max(BIC_VDEQ5[!is.na(BIC_VDEQ5)])),
     xlab= expression("The number of components " (italic("g"))),ylab="BIC values",xaxt='n',
     lwd =1,lty=1,pch=1,col=Col[1],main=" The VDEQ data") #VVV - original

lines(BIC_VDEQ5[,2],type='b',lty = 1, pch = 2 ,col=Col[2], lwd =1)#EEE
lines(BIC_VDEQ5[,3],type='b',lty = 1, pch = 3 ,col=Col[3], lwd =1)#EEV
lines(BIC_VDEQ5[,4],type='b',lty = 1, pch = 4 ,col=Col[4], lwd =1)#EVE
lines(BIC_VDEQ5[,5],type='b',lty = 1, pch = 5 ,col=Col[5], lwd =1)#EVV
lines(BIC_VDEQ5[,6],type='b',lty = 1, pch = 6 ,col=Col[6], lwd =1)#VEE
lines(BIC_VDEQ5[,7],type='b',lty = 1, pch = 7 ,col=Col[7], lwd =1)#VEV
lines(BIC_VDEQ5[,8],type='b',lty = 1, pch = 8 ,col=Col[8], lwd =1)#VVE

lines(BIC_VDEQ5[,9],type='b',lty = 3, pch = 9 ,col=Col[9], lwd =1)#VEI
lines(BIC_VDEQ5[,10],type='b',lty = 3, pch = 10 ,col=Col[10], lwd =1)#EEI
lines(BIC_VDEQ5[,11],type='b',lty = 3, pch = 11 ,col=Col[11], lwd =1)#VVI
lines(BIC_VDEQ5[,12],type='b',lty = 3, pch = 12 ,col=Col[12], lwd =1)#EVI

lines(BIC_VDEQ5[,13]+DD,type='b',lty = 4, pch = 13 ,col=Col[13], lwd =1)#VII
lines(BIC_VDEQ5[,14]+DD,type='b',lty = 4, pch = 14 ,col=Col[14], lwd =1)#EII

new_num = c(1,8,5,3,7,2,6,4,10,9,12,11,14,13)
legend("topright", c("VVV", "EEE", "EEV","EVE","EVV","VEE","VEV","VVE","VEI","EEI","VVI","EVI","VII","EII")[new_num],
       col = Col[new_num], pch = c(1:14)[new_num], lty = c(rep(1,8),rep(3,4),4,4)[new_num],ncol=4,cex=0.7, bty='n') #partial
axis(1,at=c(1:5),labels=c(1:5))
text(3+0.35,BIC_VDEQ5[3,8]-10, "1st. VVE,3",cex = 0.7, col = Col[8])
text(3+0.35,BIC_VDEQ5[3,7]+5, "2nd. VEV,3",cex = 0.7, col = Col[7])
text(4+0.35,BIC_VDEQ5[4,8]+2, "3rd. VVE,4",cex = 0.7, col = Col[8])
dev.off()
###==========fig6================================###

yhat0 = fit.VDEQ0$post.pred #VVV
yhat1 = fit.VDEQ1$post.pred #VVE

postscript(paste(SPATH, '/Results/fig6.eps', sep=''), width=20, height=30)

nf = layout(matrix(c(16,1:5,17,6:10,18,11:15),6,3),widths = c(6.5,6.5,6.5), heights =c(2,5,5,5,5,5),T)
par(mar=c(2.3,2.5,0.5,0.5),pty='s',oma=c(0,0,0,0), cex.axis=0.9, las=2,xpd = T)
p=ncol(Yc)

for(i in 1: p){
  plot(0:1,0:1, ylim=c(at[i]-0.1, lim[i]), xlim=c(min(Yc[,i])-0.5, max(Yc[,i])+0.5), 
       type='n', xlab='', ylab='', axes=T, main='', yaxt='n')
  if(i==2){
    abline(v=par("usr")[1]-4.75, h=par("usr")[3:4]+c(-0.655,0.15), lwd=3,col='red')
  }
  axis(2, seq(0,lim[i],length = 6))
  h1 = hist(Yc[,i], breaks=15, plot=F)
  a1 = sum(h1$breaks>=cutof[i])
  a2 = length(h1$counts) + 1
  hist(Yc[,i], breaks=15, prob=T, ylim=c(0, lim[i]), col=c(rep(gray(.6),a2-a1), rep(0,a1)), border=c(rep(gray(.8),a2-a1), rep(1,a1)), xlab='', ylab='Density', add=T)
  text(xli[i], alim[i], ysel[i], font=2, cex=1.2)
  boxplot(Yc[,i], xlab='', ylab='', axes=F, horizontal=T, ylim=c(min(Yc[,i]), max(Yc[,i])), add=T, at=at[i], col=0, cex=0.5, pch=16,  boxwex = bx[i])
}

lim.new = c(0.8,0.6,0.75,0.9,0.8)
alim.new = c(0.64,0.48,0.6,0.72,0.64)
at.new = c(-0.2,-0.15,-0.1,-0.1,-0.1)


#VVV
for(i in 1:p){
  plot(0:1,0:1, ylim=c(at.new[i]-0.1, lim.new[i]), xlim=c(min(c(yhat1[,i],yhat0[,i])), max(c(yhat1[,i],yhat0[,i]))),
       type='n', xlab='', ylab='', axes=T, main='', yaxt='n')
  if(i==2){
    abline(h=par("usr")[3:4]+c(-0.26,0.06), lwd=3,col='red')
  }
  axis(2, seq(0,lim.new[i],length = 6))
  h1 = hist(yhat0[,i], breaks=15, plot=F)
  a1 = sum(h1$breaks>=cutof[i])
  a2 = length(h1$counts) + 1
  hist(yhat0[,i], breaks=15, prob=T, ylim=c(0, lim.new[i]), col=c(rep(4,a2-a1), rep(0,a1)), border=c(rep("cyan",a2-a1), rep(1,a1)), xlab='', ylab='Density', add=T)
  text(xli[i], alim.new[i], ysel[i], font=2, cex=1.2)
  boxplot(yhat0[,i], xlab='', ylab='', axes=F, horizontal=T, ylim=c(min(yhat0[,i]), max(yhat0[,i])), add=T, at=at.new[i], col=0, cex=0.5, pch=16,  boxwex = bx[i])
}
#VVE

for(i in 1:p){
  plot(0:1,0:1, ylim=c(at.new[i]-0.1, lim.new[i]), xlim=c(min(c(yhat1[,i],yhat0[,i])), max(c(yhat1[,i],yhat0[,i]))), 
       type='n', xlab='', ylab='', axes=T, main='', yaxt='n')
  if(i==2){
    abline(v=par("usr")[2]+3, h=par("usr")[3:4]+c(-0.26,0.06), lwd=3,col='red')
  }
  axis(2, seq(0,lim.new[i],length = 6))
  h1 = hist(yhat1[,i], breaks=15, plot=F)
  a1 = sum(h1$breaks>=cutof[i])
  a2 = length(h1$counts) + 1
  hist(yhat1[,i], breaks=15, prob=T, ylim=c(0, lim.new[i]), col=c(rep(3,a2-a1), rep(0,a1)), border=c(rep("light green",a2-a1), rep(1,a1)), xlab='', ylab='Density', add=T)
  text(xli[i], alim.new[i], ysel[i], font=2, cex=1.2)
  boxplot(yhat1[,i], xlab='', ylab='', axes=F, horizontal=T, ylim=c(min(yhat1[,i]), max(yhat1[,i])), add=T, at=at.new[i], col=0, cex=0.5, pch=16,  boxwex = bx[i])
}
plot(1, type='n',axes=F, xlab="", ylab="")
title(main = "(a) Original data", line = -2, cex.main=1.1)
plot(1, type='n',axes=F, xlab="", ylab="")
title(main = "(b) Predicted by (VVV,3)", line = -2, cex.main=1.1)
plot(1, type='n',axes=F, xlab="", ylab="")
title(main = "(c) Predicted by (VVE,3)", line = -2, cex.main=1.1)
dev.off()
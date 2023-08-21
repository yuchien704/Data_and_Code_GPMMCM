load(paste(SPATH, "/Data/faithful.RData", sep = ""))

BIC[,c(13,14)] = BIC[,c(13,14)]/3+1200
BIC[1,c(1:8)] = BIC[1,c(1:8)]-150
BIC[1,c(9:12)] = BIC[1,c(9:12)]-380

BIC[3,12] = BIC[3,12]-15
BIC[3,11] = BIC[3,11]-10
BIC[2,1] = BIC[2,1]-5

#fig3
postscript(paste(SPATH, '/Results/fig3.eps', sep=''), width=6, height=30)
par(mar=c(4,4,2,0.5), oma = c(0, 0, 4, 0),pty='s')
Col=c("Dark Orchid4","Purple","Orange","Dark Red","Red","Magenta","Cyan","Dark Green","Royal Blue","Dark Goldenrod","Blue","Violet Red","Dark Violet","Green")
Pch=c(1:14)
plot(NA, xlim=c(0.5,5),ylim=c(min(BIC[!is.na(BIC)]), max(BIC[!is.na(BIC)])),
     xlab= expression("The number of components " (italic("g"))),ylab="BIC values",axes = F,main="The old faithful data")
for(i in 1:14)
  lines(1:5,BIC[,i], type = "b" ,col=Col[i] ,pch=Pch[i])

new_num = c(1,8,5,3,7,2,6,4,10,9,12,11,14,13)
legend("topright", c("VVV", "EEE", "EEV","EVE","EVV","VEE","VEV","VVE","VEI","EEI","VVI","EVI","VII","EII")[new_num],
       col=Col[new_num],pch=Pch[new_num],ncol=4,cex=0.85,bty ='n')

box(bty="]");lines(c(0,0.5),c(2263,2263), lwd=4, col="white");lines(c(0,0.5),c(1930,1930), lwd=4, col="white")
axis(1)
axis(2, at = c(2050,2000,1950), labels = c(2040,2000,1960), pos = 0.5)
lines(c(0.5,0.5), c(2265, 1930))
lines(c(0.5,0.5), c(2130, 2140), lwd=3, type = "l", col='white')
lines(c(0.5,0.8), c(2130, 2133), type = "l")
lines(c(0.8,0.3), c(2133, 2136), type = "l")
lines(c(0.3,0.5), c(2136, 2140), type = "l")
axis(2, at=c(2240,2200,2160), labels=c(2990,2950,2910), pos = 0.5)
text(3+0.3,BIC[3,12]-5, "1st.  EVI,3",cex = 0.7, col = Col[12])
text(3-0.3,BIC[3,11]-5, "2nd. \n VVI,3",cex = 0.7, col = Col[11])
text(2-0.2,BIC[2,1]-5, "3rd. \n VVV,2",cex = 0.7, col = Col[1])
dev.off()



#==========================================fig4========================#

m=30
tmp1=range(na.omit(Ycm[,1]))
tmp2=range(na.omit(Ycm[,2]))
xx=expand.grid(x1<-seq(tmp1[1],tmp1[2],length=m),x2<-seq(tmp2[1],tmp2[2],length=m))
den=0
for(i in 1:g) den=den+mvtnorm::dmvnorm(xx, mean=fit1$para$mu[,i], sigma=fit1$para$Sigma[,,i])
den=matrix(den,m,m)
p.clus = fit1$post.clus
y.pred = fit1$post.pred

postscript(paste(SPATH, '/Results/fig4.eps', sep=''), width=8, height=10)
mat <- matrix(c(2, 4, 1, 3), 2, 2,byrow=TRUE)
nf <- layout(mat, widths=c(4,2), heights=c(2,4), respect = TRUE)
Col_0=c("pink", "lightgreen", "light blue")
par(mar=c(4, 4, 0, 0), las=0)
plot(Ycm[-c(na.posi,cen_11,cen_21),], pch = 2*p.clus[-c(na.posi,cen_11,cen_21)]+13,xlim = c(min(Y[,1]),max(Y[,1])), ylim = c(min(Y[,2]),max(Y[,2])),cex.lab=1.2,cex=1.5,col = Col_0[p.clus[-c(na.posi,cen_11,cen_21)]])

cm.posi = na.posi[which(rowSums(cen[na.posi,],na.rm = T)!=0)]
for(i in 1:length(cm.posi)){
  if(y.pred[cm.posi[i],1]>cen_11_value){
    y.pred[cm.posi[i],1]=cen_11_value
  }else if(y.pred[cm.posi[i],2]<cen_21_value){
    y.pred[cm.posi[i],2]=cen_21_value
  }
}
Col_1=c("red", "darkgreen", "blue")
points(y.pred[na.posi,], xlim = c(min(Y[,1]),cen_11_value), ylim = c(cen_21_value,max(Y[,2])), pch = 4,cex=2, col = Col_1[p.clus[na.posi]])
points(Ycm[c(cen_11,cen_21), ], pch = 20,cex=2, col = "grey")
abline(h = cen_21_value, v = cen_11_value, lty = 2)
contour(x1, x2, den, labcex = 0.3, add=T, col=gray(.6), lwd=0.005, drawlabels = T, nlevels=20)

#hist. of x
par(mar=c(0, 4, 4, 0))
xcmhist=hist(x_cm,plot=FALSE,breaks = seq(min(x),max(x),length=21))
barplot(c(xcmhist$counts),axes=FALSE, space=0, col = c(rep(gray(.95),17),rep(gray(.7),17)))

#hist. of y
ycmhist=hist(y_cm,plot=FALSE,breaks = seq(min(y),max(y),length=21))
par(mar=c(4, 0, 0, 1))
barplot(c(ycmhist$counts),axes=FALSE, space=0, horiz=TRUE, col = c(rep(gray(.7),4),rep(gray(.95),16)))
plot.new()
legend("bottom",c("observed values","censored values", "imputed missing values in group 1", "imputed missing values in group 2", "imputed missing values in group 3")
       , pch=c(15,15,4, 4, 4), col = c(gray(.95),gray(.7), Col_1), cex=0.8, bty = 'n')
legend("center",c("Group 1","Group 2", "Group 3"), 
       col = c(Col_0),pch = c(15,17,19), cex=0.8, bty = 'n',pt.cex=1.2, ncol=3)
dev.off()
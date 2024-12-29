library(mvtnorm)
#library(mixtools)
#library(mixture)
library(mclust)
library(MomTrunc)
source(paste(SPATH, "/Function/gen_cen_na_new.R", sep = ""))

gen.tmix = function(ni, mu, Sigma, nu)                # generate Y directly
{
  g = ncol(mu); 
  Y = list(g)
  for(i in 1: g) Y[[i]] = t(mu[,i]+t(rmvt(ni[i], Sigma[,,i], df=nu[i])))
  y = NULL
  for(i in 1: g) y = rbind(y, Y[[i]])
  Y.gen = y
  return(Y.gen) 
}

set.seed(0123)
n=300; p=3; g=3; 
ni=c(n/3, n/3, n/3)
w=c(1/3,1/3,1/3)
# mu=matrix(c(-3, -4, -3, 2, 0, -1, 5, 2, 6), 3, 3)
mu = matrix(c(-9, 7, -4, 0, 1, 2, 8, -5, 6), 3, 3)
Sigma = array(NA, dim = c(3, 3, 3))
lambda1 = 16; lambda2 = 4; lambda3 = 9
theta = pi/3
D = matrix(c(1,0,0,0,cos(theta),sin(theta),0,-sin(theta),cos(theta)),3,3)
A1 = diag(c(1, 2, 0.5))
A2 = diag(c(1, 3, 1/3))
A3 = diag(c(1, 4, 1/4))
Sigma[, , 1] = lambda1 * D %*% A1 %*% t(D)
Sigma[, , 2] = lambda2 * D %*% A2 %*% t(D)
Sigma[, , 3] = lambda3 * D %*% A3 %*% t(D)
nu=c(4, 7, 10)

Y = gen.tmix(ni , mu, Sigma, nu) 
true_clus = rep(c(3,2,1), each = ni[1])

# ----------- generate the censored and missing values --------------#
Data = gener.cen.na(Y, cen.type = c(0, 1, 2),cen.rate=c(0, 0.1, 0.1),na.rate = 0.1)
Ycm = as.matrix(Data$Data)
cen = Data$cen
cen_1 = Ycm[which(cen[,2]==1),2][1]
cen_2 = Ycm[which(cen[,3]==2),3][1]
na.posi = which(rowSums(is.na(Ycm[,2:3]))!=0)

#--------- The scatter plot with variable contains left/right-censored values -------#
postscript(paste(SPATH, '/Results/Fig.S2(b).eps', sep=''), width=10, height=10) # break = 20
mat <- matrix(c(2, 4, 1, 3), 2, 2,byrow=TRUE)
nf <- layout(mat, widths=c(4,2), heights=c(2,4), respect = TRUE)
par(mar=c(4, 4, 0, 0), las=0)
plot(Ycm[,2:3], xlab='Variable 2', ylab='Variable 3',xlim=c(min(Y[,2]),max(Y[,2])), ylim=c(min(Y[,3]),max(Y[,3])), col=true_clus+1, pch=true_clus-1)
points(Y[na.posi,2:3], xlim=c(min(Y[,2]),max(Y[,2])), ylim=c(min(Y[,3]),max(Y[,3])), pch = 4)
abline(h=cen_2, v=cen_1, col='grey', lty = 2, lwd=2)

par(mar=c(0, 9, 1, 0))
x1 = list(x11 = Ycm[true_clus==1, 2], x12 = Ycm[true_clus==2, 2], x13 = Ycm[true_clus==3, 2])
h3x = lapply(x1, hist, breaks = seq(min(Ycm[,2], na.rm = T), max(Ycm[,2], na.rm = T), length=20), plot=F)
t3 = rbind(h3x[[1]]$counts, h3x[[2]]$counts, h3x[[3]]$counts)
rownames(t3)=names(h3x)
colnames(t3)=h3x[[1]]$mids
a3 = which(h3x$x11$breaks>=cen_1)[-20]
barplot(t3[1,a3], ylim=c(0,max(t3)+5),col=2, border="pink", xaxt='n', las=1, space=0, xaxt='n', yaxt='n')
barplot(t3[2,a3], ylim=c(0,max(t3)+5),col=3, border="lightgreen",  xaxt='n', las=1, space=0, add=T, xaxt='n', yaxt='n')
barplot(t3[3,a3], ylim=c(0,max(t3)+5),col=0, border="blue",xaxt='n', las=1, space=0, add=T, xaxt='n', yaxt='n')

par(mar=c(4, 0, 5, 1))
y1 = list(y11 = Ycm[true_clus==1, 3], y12 = Ycm[true_clus==2, 3], y13 = Ycm[true_clus==3, 3])
h3y = lapply(y1, hist, breaks = seq(min(Ycm[,3], na.rm = T), max(Ycm[,3], na.rm = T), length=20), plot=F)
t3y = rbind(h3y[[1]]$counts, h3y[[2]]$counts, h3y[[3]]$counts)
rownames(t3y)=names(h3y)
colnames(t3y)=h3y[[1]]$mids
a4 = which(h3y$y11$breaks<=cen_2)[-20]
barplot(t3y[1,a4], xlim=c(0,max(t3y)+5),  horiz=T, col=2, border="pink", xaxt='n', las=1, space=0, xaxt='n', yaxt='n')
barplot(t3y[2,a4], xlim=c(0,max(t3y)+5),  horiz=T, col=3, border="lightgreen", xaxt='n', las=1, space=0, add=T, xaxt='n', yaxt='n')
barplot(t3y[3,a4], xlim=c(0,max(t3y)+5),  horiz=T, col=0, border="blue", xaxt='n', las=1, space=0, add=T, xaxt='n', yaxt='n')

plot.new()
legend("bottomleft",c("Group 1","Group 2","Group 3", "missing values"), col = c(2,3,4,1), pch = c(0,1,2,4), cex=0.9, bty='n')
title('(b) Well separated mixture samples', outer=T, line=-1)
dev.off()


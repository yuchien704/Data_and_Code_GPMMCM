library(mvtnorm)
library(mixtools)
library(mixture)
library(mclust)
library(MomTrunc)
source(paste(SPATH, "/Function/gen_cen_na.R", sep = ""))

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

set.seed(1234)
n=450; p=2; G=3; 
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

Y = gen.nmix(ni , mu, Sigma) 
Data = gener.cen.na(Y, cen.type = c(1, 2),cen.rate=c(0.1,0.1),na.rate = 0.1)
Ycm = as.matrix(Data$Data)
cen = Data$cen
cen_1 = Ycm[which(cen[,1]==1),1][1]
cen_2 = Ycm[which(cen[,2]==2),2][1]
na.posi = which(rowSums(is.na(Ycm))==1)
clus = rep(c(1,2,3), each = 150)

postscript(paste(SPATH, '/Results/fig6.eps', sep=''), width=6, height=30)
mat <- matrix(c(2, 4, 1, 3), 2, 2,byrow=TRUE)
nf <- layout(mat, widths=c(4,2), heights=c(2,4), respect = TRUE)

par(mar=c(4, 4, 0, 0), las=0)
plot(Ycm,  xlab = 'Variable 1', ylab = 'Variable 1',xlim=c(min(Y[,1]),max(Y[,1])), ylim=c(min(Y[,2]),max(Y[,2])), col = clus+1, pch = clus-1)
points(Y[na.posi,], xlim=c(min(Y[,1]),max(Y[,1])), ylim=c(min(Y[,2]),max(Y[,2])), pch = 4)
abline(h = cen_2, v = cen_1, col='grey', lty = 2, lwd=2)

par(mar=c(0, 10, 0.5, 0))
x1 = list(x11 = Ycm[clus==1, 1], x12 = Ycm[clus==2, 1], x13 = Ycm[clus==3, 1])
h3x = lapply(x1, hist, breaks = seq(min(Ycm[,1], na.rm = T), max(Ycm[,1], na.rm = T), length=30), plot=F)
t3 = rbind(h3x[[1]]$counts, h3x[[2]]$counts, h3x[[3]]$counts)
rownames(t3)=names(h3x)
colnames(t3)=h3x[[1]]$mids
a3 = which(h3x$x11$breaks>=cen_1)[-30]
barplot(t3[1,a3], ylim=c(0,max(t3)),col=2, border="pink", xaxt='n', las=1, space=0, xaxt='n', yaxt='n')
barplot(t3[2,a3], ylim=c(0,max(t3)),col=3, border="lightgreen",  xaxt='n', las=1, space=0, add=T, xaxt='n', yaxt='n')
barplot(t3[3,a3], ylim=c(0,max(t3)),col=0, border="blue",xaxt='n', las=1, space=0, add=T, xaxt='n', yaxt='n')

par(mar=c(4, 0, 4, 1))
y1 = list(y11 = Ycm[clus==1, 2], y12 = Ycm[clus==2, 2], y13 = Ycm[clus==3, 2])
h3y = lapply(y1, hist, breaks = seq(min(Ycm[,2], na.rm = T), max(Ycm[,2], na.rm = T), length=30), plot=F)
t3y = rbind(h3y[[1]]$counts, h3y[[2]]$counts, h3y[[3]]$counts)
rownames(t3y)=names(h3y)
colnames(t3y)=h3y[[1]]$mids
a4 = which(h3y$y11$breaks<=cen_2)[-30]
barplot(t3y[1,a4], xlim=c(0,max(t3y)),  horiz=T, col=2, border="pink", xaxt='n', las=1, space=0, xaxt='n', yaxt='n')
barplot(t3y[2,a4], xlim=c(0,max(t3y)),  horiz=T, col=3, border="lightgreen", xaxt='n', las=1, space=0, add=T, xaxt='n', yaxt='n')
barplot(t3y[3,a4], xlim=c(0,max(t3y)),  horiz=T, col=0, border="blue", xaxt='n', las=1, space=0, add=T, xaxt='n', yaxt='n')

plot.new()
legend("topleft",c("Group 1","Group 2","Group 3", "missing values"), col = c(2,3,4,1), pch = c(0,1,2,4), cex=0.9, bty='n')
dev.off()
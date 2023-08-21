library(mvtnorm)
library(tmvtnorm)
library(mclust)
library(MomTrunc)
library(mixtools)
library("RColorBrewer")
source(paste(SPATH, "/Function/gen_cen_na.R", sep = ""))

vech.posi=function(dim) cbind(rep(1:dim, 1:dim), unlist(mapply(':', 1, 1:dim)))

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

my.ellipse = function(mu, Sigma){
  ellipse(mu[,1],Sigma[,,1])
  ellipse(mu[,2],Sigma[,,2])
  ellipse(mu[,3],Sigma[,,3])
  points(mu[1,1],mu[2,1], pch=20)
  points(mu[1,2],mu[2,2], pch=20)
  points(mu[1,3],mu[2,3], pch=20)
} 

set.seed(13)
n=1800; p=2; g=3
Sigma = array(NA, dim = c(p, p, g))
w=c(1/3,1/3.1/3)
ni=c(n/3, n/3, n/3)

postscript(paste(SPATH, '/Results/fig1.eps', sep=''), width=6, height=30)

mat = matrix(c(17,17,17,17,1:8,18,18,18,18,9:12,19,19,19,19,13:16),ncol = 4, byrow = T)
nf = layout(mat ,widths = c(3,3,3,3), heights =c(1,3,3,1,3,1,3),T)
par(mar=c(2,2,2,2), oma = c(0.5,0.5,0.5,0.5),pty='s',cex.axis=0.9, las=2)
#====================================#

#function_VVV
mu = matrix(c(-2, -5, 7, 1, -8, 6), p, g)
Sigma[, , 1] = matrix(c(3,-1,-1,1), p, p)
Sigma[, , 2] = matrix(c(1,0.5,0.5,0.5), p, p)
Sigma[, , 3] = matrix(c(2,-0.25,-0.25,2), p, p)

Y = gen.nmix(ni, mu, Sigma)
Data = gener.cen.na(Y, cen.type = c(0, 0),cen.rate=c(0,0),na.rate = 0)
Ycm = as.matrix(Data$Data)
cen = Data$cen
plot(Ycm,xlab= "",ylab="", main= "VVV",col=rep(brewer.pal(3, "Pastel1"),each=600),cex=0.6)
my.ellipse(mu,Sigma)
#==================(2EVV)=====================
mu = matrix(c(-2, -5, 7, 1, -8, 4), p, g)
Sigma = array(NA, dim = c(p, p, g))
lambda= 1 

theta1 = pi/4
theta2 = -pi/9
theta3 = pi
D1 = matrix(c(cos(theta1),-sin(theta1),sin(theta1),cos(theta1)),p,p)
D2 = matrix(c(cos(theta2),-sin(theta2),sin(theta2),cos(theta2)),p,p)
D3 = matrix(c(cos(theta3),-sin(theta3),sin(theta3),cos(theta3)),p,p)

A1 = diag(c(3,1/3))
A2 = diag(c(6,1/6))
A3 = diag(c(4,1/4))

Sigma[, , 1] = lambda*D1%*%A1%*%t(D1)
Sigma[, , 2] = lambda*D2%*%A2%*%t(D2)
Sigma[, , 3] = lambda*D3%*%A3%*%t(D3)
Y = gen.nmix(ni , mu, Sigma) 
Data = gener.cen.na(Y, cen.type = c(0, 0),cen.rate=c(0,0),na.rate = 0)
Ycm = as.matrix(Data$Data)
cen = Data$cen
plot(Ycm,xlab= "",ylab="", main= "EVV",col=rep(brewer.pal(3, "Pastel1"),each=600),cex=0.6)
my.ellipse(mu,Sigma)
#==================(3VEV)=====================
mu = matrix(c(-6, -5, 7, 1, -8, 6), p, g)
Sigma = array(NA, dim = c(p, p, g))

lambda1 = 1
lambda2 = 6
lambda3 = 3
theta1 = pi/4
theta2 = pi/3
theta3 = -pi/4

D1 = matrix(c(cos(theta1),-sin(theta1),sin(theta1),cos(theta1)),p,p)
D2 = matrix(c(cos(theta2),-sin(theta2),sin(theta2),cos(theta2)),p,p)
D3 = matrix(c(cos(theta3),-sin(theta3),sin(theta3),cos(theta3)),p,p)
A = matrix(c(3,0,0,1/3), p , p)

Sigma[, , 1] = lambda1*D1%*%A%*%t(D1)
Sigma[, , 2] = lambda2*D2%*%A%*%t(D2)
Sigma[, , 3] = lambda3*D3%*%A%*%t(D3)

Y = gen.nmix(ni , mu, Sigma) 
Data = gener.cen.na(Y, cen.type = c(0, 0),cen.rate=c(0,0),na.rate = 0)
Ycm = as.matrix(Data$Data)
cen = Data$cen
plot(Ycm,xlab= "",ylab="", main= "VEV",col=rep(brewer.pal(3, "Pastel1"),each=600),cex=0.6)
my.ellipse(mu,Sigma)
#==================(4VVE)=====================
mu = matrix(c(-6, -5, 7, 1, -8, 6), p, g)
Sigma = array(NA, dim = c(p, p, g))

lambda1 = 1
lambda2 = 6
lambda3 = 3

theta1 = pi/4
D = matrix(c(cos(theta1),-sin(theta1),sin(theta1),cos(theta1)),p,p)

A1 = diag(c(2,1/2))
A2 = diag(c(9,1/9))
A3 = diag(c(4,1/4))

Sigma[, , 1] = lambda1*D%*%A1%*%t(D)
Sigma[, , 2] = lambda2*D%*%A2%*%t(D)
Sigma[, , 3] = lambda3*D%*%A3%*%t(D)
Y = gen.nmix(ni , mu, Sigma) 
Data = gener.cen.na(Y, cen.type = c(0, 0),cen.rate=c(0,0),na.rate = 0)
Ycm = as.matrix(Data$Data)
cen = Data$cen
plot(Ycm,xlab= "",ylab="", main= "VVE",col=rep(brewer.pal(3, "Pastel1"),each=600),cex=0.6)
my.ellipse(mu,Sigma)
#==================(5EEV)=====================
mu = matrix(c(-2, -5, 7, 1, -8, 6), p, g)
Sigma = array(NA, dim = c(p, p, g))

lambda=1
theta1 = pi/4
theta2 = pi/3
theta3 = -pi/4

D1 = matrix(c(cos(theta1),-sin(theta1),sin(theta1),cos(theta1)),p,p)
D2 = matrix(c(cos(theta2),-sin(theta2),sin(theta2),cos(theta2)),p,p)
D3 = matrix(c(cos(theta3),-sin(theta3),sin(theta3),cos(theta3)),p,p)
A = matrix(c(4,0,0,1/4), p , p)

Sigma[, , 1] = lambda*D1%*%A%*%t(D1)
Sigma[, , 2] = lambda*D2%*%A%*%t(D2)
Sigma[, , 3] = lambda*D3%*%A%*%t(D3)
Y = gen.nmix(ni , mu, Sigma) 
Data = gener.cen.na(Y, cen.type = c(0, 0),cen.rate=c(0,0),na.rate = 0)
Ycm = as.matrix(Data$Data)
cen = Data$cen
plot(Ycm,xlab= "",ylab="", main= "EEV",col=rep(brewer.pal(3, "Pastel1"),each=600),cex=0.6)
my.ellipse(mu,Sigma)
#==================(6EVE)=====================
mu = matrix(c(-6, -5, 7, 1, -8, 6), p, g)
Sigma = array(NA, dim = c(p, p, g))
lambda= rep(1, g) 
lambda= 1
theta1 = pi/4
D = matrix(c(cos(theta1),-sin(theta1),sin(theta1),cos(theta1)),p,p)

A1 = diag(c(2,1/2))
A2 = diag(c(9,1/9))
A3 = diag(c(4,1/4))

Sigma[, , 1] = lambda*D%*%A1%*%t(D)
Sigma[, , 2] = lambda*D%*%A2%*%t(D)
Sigma[, , 3] = lambda*D%*%A3%*%t(D)
Y = gen.nmix(ni , mu, Sigma) 
Data = gener.cen.na(Y, cen.type = c(0, 0),cen.rate=c(0,0),na.rate = 0)
Ycm = as.matrix(Data$Data)
cen = Data$cen
plot(Ycm,xlab= "",ylab="", main= "EVE",col=rep(brewer.pal(3, "Pastel1"),each=600),cex=0.6)
my.ellipse(mu,Sigma)

#==================(7VEE)=====================
mu = matrix(c(-6, -5, 7, 1, -8, 6), p, g)
Sigma = array(NA, dim = c(p, p, g))
lambda1=1
lambda2=6
lambda3=3
theta = pi/4
D = matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),p,p)
A = matrix(c(4,0,0,1), p , p)

Sigma[, , 1] = lambda1*D%*%A%*%t(D)
Sigma[, , 2] = lambda2*D%*%A%*%t(D)
Sigma[, , 3] = lambda3*D%*%A%*%t(D)

Y = gen.nmix(ni , mu, Sigma) 
Data = gener.cen.na(Y, cen.type = c(0, 0),cen.rate=c(0,0),na.rate = 0)
Ycm = as.matrix(Data$Data)
cen = Data$cen
plot(Ycm,xlab= "",ylab="", main= "VEE",col=rep(brewer.pal(3, "Pastel1"),each=600),cex=0.6)
my.ellipse(mu,Sigma)
#==================(8EEE)=====================
mu = matrix(c(-6, -5, 7, 1, -8, 6), p, g)
Sigma = array(NA, dim = c(p, p, g))
Sigma[, , 1] = matrix(c(3,-1,-1,1), p, p)
Sigma[, , 3] = Sigma[, , 2] = Sigma[, , 1] 

Y = gen.nmix(ni , mu, Sigma) 
head(Y)
true.clus=rep(1:g, ni) 
Data = gener.cen.na(Y, cen.type = c(0, 0),cen.rate=c(0,0),na.rate = 0)
Ycm = as.matrix(Data$Data)
cen = Data$cen
plot(Ycm,xlab= "",ylab="", main= "EEE",col=rep(brewer.pal(3, "Pastel1"),each=600),cex=0.6)
my.ellipse(mu,Sigma)

#==================(9VVI)=====================
mu = matrix(c(-2, -10, 7, 1, -8, 6), p, g)
Sigma = array(NA, dim = c(p, p, g))

lambda1=1
lambda2=3
lambda3=6

B1 = diag(c(9,1))
B2 = diag(c(1,9))
B3 = diag(c(2,2))
Sigma[, , 1] = lambda1*B1 
Sigma[, , 2] = lambda2*B2 
Sigma[, , 3] = lambda3*B3 
Y = gen.nmix(ni , mu, Sigma) 
Data = gener.cen.na(Y, cen.type = c(0, 0),cen.rate=c(0,0),na.rate = 0)
Ycm = as.matrix(Data$Data)
cen = Data$cen
plot(Ycm,xlab= "",ylab="", main= "VVI",col=rep(brewer.pal(3, "Pastel1"),each=600),cex=0.6)
my.ellipse(mu,Sigma)

#==================(10EVI)=====================
mu = matrix(c(-2, -5, 7, 1, -8, 6), p, g)
Sigma = array(NA, dim = c(p, p, g))

lambda=2

B1 = diag(c(9,1))
B2 = diag(c(1,9))
B3 = diag(c(3,3))
Sigma[, , 1] = lambda*B1 
Sigma[, , 2] = lambda*B2 
Sigma[, , 3] = lambda*B3 
Y = gen.nmix(ni , mu, Sigma) 
Data = gener.cen.na(Y, cen.type = c(0, 0),cen.rate=c(0,0),na.rate = 0)
Ycm = as.matrix(Data$Data)
cen = Data$cen

plot(Ycm,xlab= "",ylab="", main= "EVI",col=rep(brewer.pal(3, "Pastel1"),each=600),cex=0.6)
my.ellipse(mu,Sigma)

#==================(11VEI)=====================
mu = matrix(c(-2, -5, 7, 1, -8, 6), p, g)
Sigma = array(NA, dim = c(p, p, g))

lambda1=1
lambda2=3
lambda3=6
B = diag(c(1,6))

Sigma[, , 1] = lambda1*B 
Sigma[, , 2] = lambda2*B 
Sigma[, , 3] = lambda3*B
Y = gen.nmix(ni , mu, Sigma) 
Data = gener.cen.na(Y, cen.type = c(0, 0),cen.rate=c(0,0),na.rate = 0)
Ycm = as.matrix(Data$Data)
cen = Data$cen
plot(Ycm,xlab= "",ylab="", main= "VEI",col=rep(brewer.pal(3, "Pastel1"),each=600),cex=0.6)
my.ellipse(mu,Sigma)
#==================(12EEI)=====================

mu = matrix(c(-2, -5, 7, 1, -8, 6), p, g)
Sigma = array(NA, dim = c(p, p, g))
lambda=2
B = diag(c(1,6))

Sigma[, , 1] = lambda*B
Sigma[, , 3] = Sigma[, , 2] = Sigma[, , 1] 
Y = gen.nmix(ni , mu, Sigma) 
Data = gener.cen.na(Y, cen.type = c(0, 0),cen.rate=c(0,0),na.rate = 0)
Ycm = as.matrix(Data$Data)
cen = Data$cen

plot(Ycm,xlab= "",ylab="", main= "EEI",col=rep(brewer.pal(3, "Pastel1"),each=600),cex=0.6)
my.ellipse(mu,Sigma)
#==================(13VII)=====================
plot.new()
mu = matrix(c(-2, -5, 7, 1, -8, 6), p, g)
Sigma = array(NA, dim = c(p, p, g))

lambda1=1
lambda2=4
lambda3=9
I= diag(rep(1, p))

Sigma[, , 1] = lambda1*I
Sigma[, , 2] = lambda2*I
Sigma[, , 3] = lambda3*I
Y = gen.nmix(ni , mu, Sigma) 
true.clus=rep(1:g, ni) 
Data = gener.cen.na(Y, cen.type = c(0, 0),cen.rate=c(0,0),na.rate = 0)
Ycm = as.matrix(Data$Data)
cen = Data$cen
plot(Ycm,xlab= "",ylab="", main= "VII",col=rep(brewer.pal(3, "Pastel1"),each=600),cex=0.6)
my.ellipse(mu,Sigma)
#==================(14EII)=====================

mu = matrix(c(-2, -5, 7, 1, -8, 6), p, g)
Sigma = array(NA, dim = c(p, p, g))
lambda=2

I= diag(rep(1, p))

Sigma[, , 1] = lambda*I
Sigma[, , 3] =Sigma[, , 2] = Sigma[, , 1]
Y = gen.nmix(ni , mu, Sigma) 
head(Y)
true.clus=rep(1:g, ni) 
Data = gener.cen.na(Y, cen.type = c(0, 0),cen.rate=c(0,0),na.rate = 0)
Ycm = as.matrix(Data$Data)
cen = Data$cen
plot(Ycm,xlab= "",ylab="", main= "EII",col=rep(brewer.pal(3, "Pastel1"),each=600),cex=0.6)
my.ellipse(mu,Sigma)

title(main = "(a) Ellipsoidal family", outer = T, line = -4,cex.main=1.5)
title(main = "(b) Diagonal family", outer = T,line = -30.5,cex.main=1.5)
title(main = "(c) Spherical family", outer = T,line = -45,cex.main=1.5)
dev.off()
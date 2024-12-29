library(mvtnorm)
library(tmvtnorm)
library(mclust)
library(MomTrunc)
load(paste(SPATH, "/Data/hawks_new.RData", sep = ""))

bic = c(fit_EEE$model.inf[[5]],fit_VEE$model.inf[[5]],fit_EEV$model.inf[[5]],fit_EVE$model.inf[[5]],
        fit_EVV$model.inf[[5]],fit_VEV$model.inf[[5]],fit_VVE$model.inf[[5]],fit_VVV$model.inf[[5]],
        fit_EEI$model.inf[[5]],fit_VEI$model.inf[[5]],fit_EVI$model.inf[[5]],fit_VVI$model.inf[[5]],
        fit_EII$model.inf[[5]],fit_VII$model.inf[[5]])

names(bic) = c("EEE","VEE","EEV","EVE","EVV","VEV","VVE","VVV","EEI","VEI","EVI","VVI","VII","EII")
Best_hawks = sort(bic)[1:3]

bic[(13:14)] = bic[(13:14)]-5000

postscript(paste(SPATH, '/Results/fig5.eps', sep=''), width=8, height=8)
mat <- matrix(1:2, 1, 2)
nf = layout(mat, widths=c(6,3), heights=2)
par(mar=c(4,4,2,0), oma = c(0,0,4,1),pty='s')
plot(bic, type = 'b',lwd = 2,axes=FALSE, xlab = "structure", ylab = "BIC values")

points(8, bic[8], col = 2, pch = 15,cex = 1.5)
points(6, bic[6], col = 3, pch = 16,cex = 1.5)
points(7, bic[7], col = 4, pch = 17,cex = 1.5)
text(8+0.8,bic[8]-50, "1st \n VVV,3",cex = 0.8, col = 2)
text(6+0.8,bic[6]-50, "2nd \n VEV,3",cex = 0.8, col = 3)
text(7+0.8,bic[7]-50, "3rd \n VVE,3",cex = 0.8, col = 4)
axis(1, 1:14, c("EEE","VEE","EEV","EVE","EVV","VEV","VVE","VVV","EEI","VEI","EVI","VVI","VII","EII"),cex = 0.3, las = 2)

box(bty="]");lines(c(0,0.6),c(44400,44400), lwd=4, col="white");lines(c(0,0.6),c(32550,32550), lwd=4, col="white")
axis(2, at = c(40000,38000,36000,34000), pos = 0.6)
lines(c(0.6,0.6), c(44500,32500))
lines(c(0.6,0.6), c(42000, 41600), lwd=3, type = "l", col='white')
lines(c(0.6,0.9), c(41600, 41725), type = "l")
lines(c(0.9,0.4), c(41725, 41850), type = "l")
lines(c(0.4,0.6), c(41850, 42000), type = "l")
axis(2, at=c(44000,43000), labels=c(49000,48000), pos = 0.6)
plot.new()
legend("left", legend = c("1st VVV,3","2nd VEV,3","3rd VVE,3"), text.col = c(2,3,4),pch = c(15,16,17),col=c(2,3,4))
dev.off()

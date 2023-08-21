load(paste(SPATH, "/Data/faithful_01.RData", sep = ""))
m=50
tmp01=range(Y[,1])+c(-0.5,0.5)
tmp02=range(Y[,2])+c(-0.5,0.5)

xx0=expand.grid(x01<-seq(tmp01[1],tmp01[2],length=m),x02<-seq(tmp02[1],tmp02[2],length=m))
den0=0
for(i in 1:g) den0=den0+fit0$para$w[i]*mvtnorm::dmvnorm(xx0, mean=fit0$para$mu[,i], sigma=fit0$para$Sigma[,,i])
den0=matrix(den0,m,m)

postscript(paste(SPATH, '/Results/fig2.eps', sep=''), width=6, height=30)

nbcol = 10
jet.colors1 <- colorRampPalette( c("#F0F0F0","#BEBEBE") )  
color1 <- jet.colors1(nbcol)
jet.colors2 <- colorRampPalette( c("#BEBEBE","#7B7B7B") ) 
color2 <- jet.colors2(nbcol)
jet.colors3 <- colorRampPalette( c("#7B7B7B","#3C3C3C") ) 
color3 <- jet.colors3(nbcol)
color<-matrix(c(color1, color2, color3), ncol=1)

x=x01; y=x02; z=den0
zfacet <- z[-1, -1] + z[-1, -length(y)] + z[-length(x), -1] + z[-length(x), -length(y)]
facetcol <- cut(zfacet, breaks=nbcol*3)
trans <- persp(x,y,z, zlim=c(0, 0.06), box=T, theta=20, lwd=.5, xlab = "waiting", 
               ylab = "eruptions", zlab = "density", expand=.8, col=color[facetcol], phi=30, border=T)
clines <- contourLines(x, y, z)
lapply(clines, function(contour){ lines(trans3d(contour$x, contour$y, 0.06, trans), col = 'grey') })

for(i in 1:n)
{
  points(trans3d(Y[i,1],Y[i,2], 0.06, trans), col = c(2:4)[p.clus0[i]],pch=2*p.clus0[i]+13,cex=0.7)
}
dev.off()

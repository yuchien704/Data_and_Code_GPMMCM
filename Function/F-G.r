FG.alg= function(p, g, S, tol.FG, max.iter.FG)
{
#begin = proc.time()[1]
R = array(NA, dim=c(p,p,g))
B = diag(p); B.old = diag(p);phi.old = 1
TT = array(NA,dim=c(2,2,g))
iter.FG=0
repeat
{
iter.FG=iter.FG+1; Q.old=diag(2)
    for(dd in 2:p)
    {
        for(ee in 1:(dd-1))
        {
            Q.new=diag(2)
            if(ee < dd)
            {
            iter.Q=0
                repeat
                {
                  iter.Q=iter.Q+1
                  M = matrix(0, 2, 2)
                  d1= d2 = rep(NA,g)
                  for(i in 1:g)
                  {
                    H = B[,c(dd,ee)]
                    TT[,,i] = t(H)%*%S[,,i]%*%H
                    d1[i] = t(Q.new[,1])%*%TT[,,i]%*%Q.new[,1]
                    d2[i] = t(Q.new[,2])%*%TT[,,i]%*%Q.new[,2]
                    M = M + (d1[i]-d2[i])/(d1[i]*d2[i])*TT[,,i]
                  }
                  normala = t(Q.new[,1])%*%M%*%Q.new[,2]
                  Q.new = eigen(M)$vector
                  B[,c(dd,ee)] = H%*%Q.new
                  #cat('iter.Q=',iter.Q,'\t Q=', Q.new, '\n')
                  if(sqrt(sum((Q.old-Q.new)^2)) < tol.FG | iter.Q == max.iter.FG)break
                  Q.old=Q.new
                }
            }
        }
    }
    Crit = sum((B-B.old)^2)
    B.old = B
    SSLF = 0
    phi.new = 1
    for(i in 1:g)
    {
      R[,,i] = t(B) %*% S[,,i] %*% B
      phi.new = phi.new * det(diag(diag(R[,,i]))) / det(R[,,i])
      SSLF = SSLF + sum((R[,,i] - diag(diag(R[,,i])))^2)
    }
    #cat('iter.F=',iter.FG,'\tCrit =', Crit, '\tphi =', phi.new, '\tSSLF =', SSLF, '\n')
    if(abs(phi.old-phi.new)< tol.FG | iter.FG == max.iter.FG) break
    phi.old=phi.new
}

La = matrix(NA, p, g)
for(i in 1:g)
{La[,i] = diag(t(B)%*%S[,,i]%*%B)}

hS = array(NA, dim=c(p,p,g))
for(i in 1:g)
{hS[,,i] = B%*%diag(La[,i],p)%*%t(B)}
end = proc.time()[1]
#cat("It took", end - begin, "seconds.\n")
return(list(phi.new=phi.new,B=B,phi.new=phi.new,R=R,hS=hS))
}
#S  Sample Covariance Matrices
#hS  MLE's of Population Covariance Matrices
#B  Coefficients of Common Principal Components
#R  The simultaneously diagonalizable matrices
 
#iris
#class.label = unique(iris[,5])
#p = 4; g = 3; 
#S= array(NA, dim=c(p,p,g))
#for(i in 1:g)
#S[,,i]=round(100*cov(iris[iris[,5]==class.label[i],-5]),4) } 
#FG.iris=FG.alg(p=p,g=g,S, tol.FG = 1e-6)
#FG.iris

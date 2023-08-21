MM1 = function(p, g, S, A, D, Sigma ,tol.MM, max.iter.MM)
{
  begin = proc.time()[1]
  iter.MM = 0
  g.old = 0
  
  repeat
  {
    iter.MM = iter.MM+1
    C.new = 0
    F1 = array(NA, dim=c(p,p,g))
    for (i in 1:g) {
      F1[,,i] = solve(A[,,i])%*%t(D)%*%S[,,i]-eigen(S[,,i])$val[1]*solve(A[,,i])%*%t(D)
      C.new = C.new + sum(diag(S[,,i]%*%D%*%solve(A[,,i])%*%t(D)))
    }
    F1.sum = apply(F1, 1:2, sum)
    P = svd(F1.sum)$u
    R = svd(F1.sum)$v
    D = R%*%t(P)
    g.new = C.new
    diff = (g.new-g.old)
    if(diff < tol.MM || iter.MM >= max.iter.MM) break
    g.old = g.new
  }
  end = proc.time()[1]
  return(list(D = D, iter = iter.MM, time = end-begin, diff = diff))
}

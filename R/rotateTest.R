# k=2
# p=20
# n=100
# Rinv=matrix(rnorm(p*n),n,p)
# H=matrix(rnorm(p*k*p*k),p*k)
# Rk=kronecker(Rinv,diag(k))
# library(microbenchmark)
# microbenchmark(rotateHessian(Rinv,H,k),times=100)
# #microbenchmark(t(Rk)%*%H%*%Rk,times=1000)
# #microbenchmark(solve(H+diag(p*k)))
# microbenchmark(fullMultiply(Rk,H),times=100)
# 

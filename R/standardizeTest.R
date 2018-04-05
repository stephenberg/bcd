# library(glmnet)
# p=5
# set.seed(33)
# beta=rnorm(p)
# n=400
# X0=matrix(rnorm(n*p),n,p)+4
# y=X0%*%beta+rnorm(n)+4
# probs=1/(1+exp(-y))
# y=rbinom(n,1,probs)
# colNorms=apply(X0*X0,2,mean)
# colMeans=apply(X0,2,mean)
# colVars=(colNorms-colMeans^2)
# scaleMat=diag(1/sqrt(colVars))
# 
# X1=X0%*%scaleMat
# 
# # colNorms=apply(X*X,2,mean)
# # colMeans=apply(X,2,mean)
# # colVars=(colNorms-colMeans^2)
# # scaleMat=diag(1/sqrt(colVars))
# fit0=glmnet(x=X0,y=y,standardize = TRUE,thresh = 10^-12,family="binomial")
# fit1=glmnet(x=X1,y=y,standardize = FALSE,thresh=10^-12,family="binomial")
# beta0=as.matrix(fit0$beta)
# beta1=((scaleMat)%*%as.matrix(fit1$beta))
# plot(beta0,beta1)
# 

# library(glmnet)
# library(microbenchmark)
# set.seed(33)
# p=10
# n=100
# k=1
# X=matrix(rnorm(n*p),ncol=p)+100
# #X[,1]=rep(1,n)
# beta=matrix(rnorm(p*k),nrow=p)
# beta[4:p,]=0
# y=X%*%beta+rnorm(n*k)
# #grpreg(X,y,group=c(1,2,3,4,5,6,7,8,9,10))
# 
# groups=as.list(1:p)
# groupLengths=sapply(groups,length)
# groupStart=c(0,cumsum(groupLengths[1:(length(groups)-1)]))
# groupEnd=cumsum(groupLengths)-1
# groupColumns=unlist(groups)-1
# 
# t=multiResponseGaussianDense(nGroups_ = p,groupStart,groupEnd,groupColumns,10^-12,rep(1,p),10000,y,X,1,k,100)
# fitglmnet=glmnet(x=X,y=y,intercept=FALSE,standardize=TRUE,thresh=10^-9)
# 
# fitglmnet$beta
# fitglmnet$lambda
# 
# t[[length(t)]][1]
# fitglmnet$lambda[1]
# 
# lambda=rep(1,p)
# for (i in 1:p){
#   xi=X[,i]
#   xi=xi
#   xi=xi/((sum(xi^2))^0.5)
#   lambda[i]=((t(xi)%*%y)^2)^0.5/n
# }
# 
# #Test
# lambda2=rep(1,p)
# for (i in 1:p){
#   xi=X[,i]
#   xi=xi-mean(xi)
#   vari=((sum(xi^2))^0.5)
#   lambda2[i]=1/vari*((t(X[,i])%*%y)^2)^0.5/n
# }
# 
# t[[length(t)]][1]
# fitglmnet$lambda[1]
# max(lambda)
# max(lambda2)
# #t=multiResponseGaussianDense(nGroups_ = p,groupStart,groupEnd,groupColumns,10^-9,c(0,rep(1,p-1)),10000,y,X,1,k,length(fitglmnet$lambda))
# 
# #
# # microbenchmark(multiResponseGaussianDense(nGroups_ = p,groupStart,groupEnd,groupColumns,10^-9,rep(1,p),10000,y,X,1,k,length(fitglmnet$lambda))
# # ,times=1)
# # microbenchmark(glmnet(x=X,y=y,family="mgaussian",intercept=FALSE,standardize = TRUE,thresh = 10^-9),times=1)

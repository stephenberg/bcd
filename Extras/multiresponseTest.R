# library(grpreg)
# library(grplasso)
# library(microbenchmark)
# 
# set.seed(33)
# tol=10^-8
# tolGrp=10^-10
# k=1
# n=50
# p=100
# X=matrix(rnorm(n*p),n,p)
# X[,1]=rep(1,n)
# beta=rnorm(p)
# groups=list(1:1,2:5,6:15,16:20,21:22,23:25,26:30)
# groups2=as.list(31:100)
# groups=c(groups,groups2)
# #X[,2]=c(rep(1,n/2),rep(0,n/2))
# #X[,3]=c(rep(0,n/2),rep(1,n/2))
# nGroups=length(groups)
# groupLengths=sapply(groups,length)
# groupStart=c(0,cumsum(groupLengths[1:(length(groups)-1)]))
# groupEnd=cumsum(groupLengths)-1
# groupColumns=unlist(groups)-1
# beta=rnorm(length(groupColumns)*k)
# beta[21:22]=0
# beta[26:p]=0
# penaltyFactor=c(0,rep(1,length(groups)-1))
# sampleWeights=rep(1,n)
# XSparse=Matrix(X,sparse=TRUE)
# eta=X%*%beta
# response=matrix(X%*%beta+rnorm(n))
# 
# 
# fit=bcdFit(X=X,y=response,k=dim(response)[2],family="gaussian",groups,penaltyFactor = penaltyFactor,tolerance = 10^-10)
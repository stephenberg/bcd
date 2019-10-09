library(grpreg)
library(grplasso)
library(microbenchmark)
library(glmnet)

set.seed(33)
tol=10^-8

k=1
n=100
p=30
X=matrix(rnorm(n*p),n,p)
X[,1]=rep(1,n)
beta=rnorm(p)/5

#the grouping:
groups=list(1:1,2:5,6:15,16:20,21:22,23:25,26:30)
nGroups=length(groups)
groupLengths=sapply(groups,length)
groupStart=c(0,cumsum(groupLengths[1:(length(groups)-1)]))
groupEnd=cumsum(groupLengths)-1
groupColumns=unlist(groups)-1
beta[21:22]=0
beta[26:30]=0
penaltyFactor=c(0,1,1,1,1,1,1)
sampleWeights=rep(1,n)
XSparse=Matrix(X,sparse=TRUE)
eta=X%*%beta
mu=exp(X%*%beta)
response=rpois(n,mu)
z=response
offset=rep(0,n)
response=matrix(response,ncol=1)

X=X
X[,1]=rep(1,n)


#fitting with this package, dense matrix, offset=0
fitDense0=bcdFit(X=X,y=response,k=1,family="poisson",groups=groups,penaltyFactor = penaltyFactor,offset=rep(0,n),lambdaMinRatio = 0.05)
# #using sample weights
fitDense1=bcdFit(X=X,y=response,k=1,family="poisson",groups=groups,penaltyFactor = penaltyFactor,sampleWeights = rep(1,n),offset=rep(0,n),lambdaMinRatio = 0.05)

# #fitting with this package, sparse matrix, offset=0
fitSparse0=bcdFit(X=Matrix(X,sparse=TRUE),y=response,k=1,family="poisson",groups=groups,penaltyFactor = penaltyFactor,offset=rep(0,n),lambdaMinRatio = 0.05)
# #using sample weights
fitSparse1=bcdFit(X=Matrix(X,sparse=TRUE),y=response,k=1,family="poisson",groups=groups,penaltyFactor = penaltyFactor,sampleWeights = rep(1,n),offset=rep(0,n),lambdaMinRatio = 0.05)

#fitting with grpreg for comparison purposes
fitgrpreg=grpreg(X[,2:dim(X)[2]],response,group=c(rep(1,4),rep(2,10),rep(3,5),rep(4,2),rep(5,3),rep(6,5)),family="poisson",lambda.min=0.05,eps=10^-12,lambda=fitDense0$lambda*sqrt(n))
ng=dim(fitgrpreg$beta)[2]
max(abs(unlist(fitDense0$beta[1:ng])-c(fitgrpreg$beta)))
max(abs(unlist(fitDense1$beta[1:ng])-c(fitgrpreg$beta)))
max(abs(unlist(fitSparse0$beta[1:ng])-c(fitgrpreg$beta)))
max(abs(unlist(fitSparse1$beta[1:ng])-c(fitgrpreg$beta)))

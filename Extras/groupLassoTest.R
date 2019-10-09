library(grpreg)
library(grplasso)
library(microbenchmark)
library(msgl)
set.seed(33)
tol=10^-8


k=1
n=300
p=60

XSparse=replicate(p,rbinom(n,1,0.05))
X=matrix(rnorm(n*p),n,p)
X[,1]=rep(1,n)
XSparse[,1]=rep(1,n)
#X=X+10
X[,1]=rep(1,n)

##########
#X=XSparse
XSparse=X
XSparse=Matrix(XSparse,sparse=TRUE)
##########

beta=rnorm(p)/3
groups=list(1,2:5,6:15,16:20,21:22,23:25,26:30)
#groups=c(groups,as.list(31:p))
groups=as.list(1:p)
#groups=list(1,2:5,6:15)

#X[,2]=c(rep(1,n/2),rep(0,n/2))
#X[,3]=c(rep(0,n/2),rep(1,n/2))
nGroups=length(groups)
groupLengths=sapply(groups,length)
groupStart=c(0,cumsum(groupLengths[1:(length(groups)-1)]))
groupEnd=cumsum(groupLengths)-1
groupColumns=unlist(groups)-1
beta=rnorm(length(groupColumns)*k)
penaltyFactor=c(0,rep(1,length(groups)-1))
sampleWeights=rep(1,n)
eta=X%*%beta
response=eta+rnorm(n)
grouping=c(rep(1,4),rep(2,10),rep(3,5),rep(4,2),rep(5,3),rep(6,5))
grouping=c(grouping,7:(length(groups)-1))
grouping=(1:(p-1))
grpregFit=NULL
microbenchmark(grpregFit<-grpreg(X = X[,2:dim(X)[2]],y=response,eps = 10^-12,group=grouping,lambda.min = 0.0001),times=1)
nLambda=dim(grpregFit$beta)[2]-1
test=NULL

testDense=NULL
testSparse=NULL
lambdaMinRatio=grpregFit$lambda[nLambda]/grpregFit$lambda[1]
testDense<-fit_bcd(X=X,
                  y=response,
                  penaltyFactor = penaltyFactor,
                  groups=groups)
microbenchmark(testDense<-fit_bcd(X=X,y=response,penaltyFactor=penaltyFactor,groups=groups,tolerance = 10^-12,lambdaMinRatio=lambdaMinRatio,nLambda = nLambda,devTol = 0.99),times=1)
microbenchmark(testSparse<-fit_bcd(X=Matrix(X,sparse=TRUE),y=response,penaltyFactor=penaltyFactor,groups=groups,tolerance = 10^-12,lambdaMinRatio=lambdaMinRatio,nLambda = nLambda),times=1)

betaMatFromList<-function(test,n){
  res=NULL
  for (i in 1:(n)){
    res=cbind(res,test[[i]])
  }
  return(res)
}
#
library(glmnet)
fitglmnet=NULL
microbenchmark(fitglmnet<-glmnet(x=XSparse[,2:p],y=response,lambda = grpregFit$lambda[1:nLambda],thresh = 10^-24),times=1)

grpregBeta=grpregFit$beta[,1:nLambda]
myBeta=betaMatFromList(testDense$beta,nLambda)
glmnetBeta=as.matrix(coef(fitglmnet))[,1:nLambda]
#msglBeta=betaMatFromList(beta_msgl,nLambda)
diffsGrpReg=myBeta-grpregBeta
diffsGlmnet=myBeta-glmnetBeta
diffsGrpNet=glmnetBeta-grpregBeta
print(max(abs(diffsGrpReg[])))
print(max(abs(diffsGlmnet[])))
print(max(abs(diffsGrpNet)))

grpregFit$lambda[1]/testDense$lambda[1]
sqrt(n)


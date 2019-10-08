# library(grpreg)
# library(grplasso)
# library(microbenchmark)
# 
# set.seed(33)
# tol=10^-8
# 
# k=1
# n=30
# p=30
# X=matrix(rnorm(n*p),n,p)+1
# X[,1]=rep(1,n)
# beta=rnorm(p)
# groups=list(1,2:5,6:15,16:20,21:22,23:25,26:30)
# #groups=list(1,2:5,6:15)
# 
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
# penaltyFactor=c(0,1,1,1,1,1,1)
# sampleWeights=rep(1,n)
# XSparse=Matrix(X,sparse=TRUE)
# eta=X%*%beta
# probs=exp(eta)/(1+exp(eta))
# z=rbinom(n,1,probs)
# response=cbind(1-z,z)
# 
# 
# fitreg=grpreg(X[,2:dim(X)[2]],z,group=c(rep(1,4),rep(2,10),rep(3,5),rep(4,2),rep(5,3),rep(6,5)),family="binomial",eps=10^-12)
# print(dim(fitreg$beta))
# nLambda=dim(fitreg$beta)[2]
# #
# #
# # test2=NULL
# # for (i in 1:length(test)){
# #   test2=cbind(test2,test[[i]][,2]-test[[i]][,1])
# # }
# # c(test2[,1:nLambda])-c(fitreg$beta[,1:nLambda])
# # fitreg$lambda/test[1:length(fitreg$lambda)]
# 
# 
# 
# u0=rep(1,n)
# Proj=diag(n)-1/n*u0%*%t(u0)
# u1=qr.Q(qr(Proj%*%X[,2:5]))
# u2=qr.Q(qr(Proj%*%X[,6:15]))
# 
# means=t(replicate(n,apply(response,2,mean)))
# 
# 
# resids=Proj%*%(response-means)*2
# means0=Proj%*%means
# v1=t(u1)%*%(response-means0)*2
# v2=t(u2)%*%(response-means0)*2
# norm1=sum(v1^2)^0.5
# norm2=sum(v2^2)^0.5
# n1=norm1/sqrt(length(v1))/n
# n2=norm2/sqrt(length(v2))/n
# 
# test=multinomialDense(p_=p,n_=n,k_=2,nGroups=length(groups),groupStart_ = groupStart
#                       ,groupEnd_ = groupEnd,groupColumns_ = groupColumns
#                       ,outerTol_ = 10^-9,innerTol_ = 10^-9,10000,10000,penaltyFactor_ = penaltyFactor
#                       ,response_ = response,X_ = X,nLambda_=nLambda)
# testSparse=multinomialSparse(p_=p,n_=n,k_=2,nGroups=length(groups),groupStart_ = groupStart
#                       ,groupEnd_ = groupEnd,groupColumns_ = groupColumns
#                       ,outerTol_ = 10^-9,innerTol_ = 10^-9,10000,10000,penaltyFactor_ = penaltyFactor
#                       ,response_ = response,X_ = Matrix(X,sparse=TRUE),nLambda_=nLambda)
# test2=NULL
# test3=NULL
# for (i in 1:length(test)){
#   test2=cbind(test2,test[[i]][,2]-test[[i]][,1])
#   test3=cbind(test3,testSparse[[i]][,2]-testSparse[[i]][,1])
# }
# 
# # microbenchmark(multinomialDense(p_=p,n_=n,k_=2,nGroups=length(groups),groupStart_ = groupStart
# # ,groupEnd_ = groupEnd,groupColumns_ = groupColumns
# # ,outerTol_ = 10^-5,innerTol_ = 10^-9,10000,10000,penaltyFactor_ = penaltyFactor
# # ,response_ = response,X_ = X,nLambda_=nLambda),times=5)
# # microbenchmark(grpreg(X[,2:dim(X)[2]],z,group=c(rep(1,4),rep(2,10),rep(3,5),rep(4,2),rep(5,3),rep(6,5)),family="binomial",eps=10^-12),times=5)
# max(abs((c(test2[,])-c(fitreg$beta[,]))))
# max(abs(test3-test2))
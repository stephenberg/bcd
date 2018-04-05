# library(grpreg)
# library(grplasso)
# library(microbenchmark)
# library(msgl)
# library(glmnet)
# set.seed(33)
# tol=10^-8
# 
# 
# config=msgl.algorithm.config(tolerance_penalized_main_equation_loop = 1e-10,
#                              tolerance_penalized_inner_loop_alpha = 1e-06,
#                              tolerance_penalized_inner_loop_beta = 0.0001,
#                              tolerance_penalized_middel_loop_alpha = 0.00001,
#                              tolerance_penalized_outer_loop_alpha = 0.00001,
#                              tolerance_penalized_outer_loop_beta = 0,
#                              tolerance_penalized_outer_loop_gamma = 1e-05,
#                              use_bound_optimization = TRUE,
#                              use_stepsize_optimization_in_penalizeed_loop = TRUE,
#                              stepsize_opt_penalized_initial_t = 1, stepsize_opt_penalized_a = 0.1,
#                              stepsize_opt_penalized_b = 0.1, max_iterations_outer = 1e+05,
#                              inner_loop_convergence_limit = 1e+05, verbose = TRUE)
# 
# k=1
# n=500
# p=500
# X=matrix(rnorm(n*p),n,p)
# X[,1]=rep(1,n)
# beta=rnorm(p)/3
# groups=list(1:1,2:5,6:15,16:20,21:22,23:25,26:30)
# groups=c(groups,as.list(31:p))
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
# penaltyFactor=c(0,rep(1,length(groups)-1))
# sampleWeights=rep(1,n)
# XSparse=Matrix(X,sparse=TRUE)
# eta=X%*%beta
# probs=exp(eta)/(1+exp(eta))
# z=rbinom(n,1,probs)
# response=cbind(1-z,z)
# grouping=c(rep(1,4),rep(2,10),rep(3,5),rep(4,2),rep(5,3),rep(6,5))
# grouping=c(grouping,7:(length(groups)-1))
# 
# 
# # u=rep(1,n)
# # projMat=diag(n)-u%*%t(u)/n
# # RinvList=list()
# # X2=X
# # for (i in 1:(length(groups)-1)){
# #   Xi=projMat%*%X[,groups[[i+1]]]
# #   RinvList[[i]]=solve(qr.R(qr(Xi)))
# #   X2[,groups[[i+1]]]=X[,groups[[i+1]]]%*%RinvList[[i]]
# # }
# 
# #microbenchmark(msglFit<-fit(X2[,2:dim(X2)[2]],factor(response[,2]),grouping=grouping,standardize=FALSE,alpha=0,lambda=0.5,d=100,intercept=TRUE,algorithm.config = config),times=2)
# grpregFit=NULL
# microbenchmark(grpregFit<-grpreg(X = X[,2:dim(X)[2]],y=response[,2],eps = 10^-15,group=grouping,family="binomial",lambda.min = 0.05),times=1)
# nLambda=dim(grpregFit$beta)[2]-1
# lambdaMinRatio=grpregFit$lambda[nLambda]/grpregFit$lambda[1]
# testDense=NULL
# testSparse=NULL
# microbenchmark(testDense<-bcdFit(X=X,y=response,k=2,family="multinomial",groups=groups,
#                             penaltyFactor = penaltyFactor,lambdaMinRatio = lambdaMinRatio,nLambda=nLambda,tolerance = 10^-13),times=1)
# microbenchmark(testSparse<-bcdFit(X=Matrix(X,sparse=TRUE),y=response,k=2,family="multinomial",groups=groups,
#                                  penaltyFactor = penaltyFactor,lambdaMinRatio = lambdaMinRatio,nLambda=nLambda,tolerance = 10^-13),times=1)
# 
# # microbenchmark(testDense<-bcdFit(X=X,y=response,k=2,family="multinomial",groups=groups,
# #                                  penaltyFactor = penaltyFactor,lambdaMinRatio = 0.0001,nLambda=100,tolerance = 10^-5),times=1)
# # microbenchmark(testSparse<-bcdFit(X=Matrix(X,sparse=TRUE),y=response,k=2,family="multinomial",groups=groups,
# #                                   penaltyFactor = penaltyFactor,lambdaMinRatio = 0.0001,nLambda=100,tolerance = 10^-5),times=1)
# 
# #fitglmnet=glmnet(x=X[,2:dim(X)[2]],y=response[,2],family="binomial",lambda.min.ratio = 0.0001,thresh = 10^-15)
# betaMatFromList<-function(test,n){
#   res=NULL
#   for (i in 1:(n)){
#     res=cbind(res,test$beta[[i]][,2]-test$beta[[i]][,1])
#   }
#   return(res)
# }
# #
# # #
# grpregBeta=grpregFit$beta[,1:nLambda]
# myBeta=betaMatFromList(testDense,length(testSparse$beta))
# #glmnetBeta=as.matrix(coef(fitglmnet)[,1:58])
# # # #msglBeta=betaMatFromList(beta_msgl,nLambda)
# diffsGrpReg=myBeta-grpregBeta
# print(max(abs(diffsGrpReg)))

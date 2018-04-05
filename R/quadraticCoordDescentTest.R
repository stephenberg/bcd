# set.seed(33)
# k=1
# n=5000
# p=15
# X=matrix(rnorm(n*p),n,p)
# beta=rnorm(p)
# groups=list(1:10,11:15)
# groupLengths=sapply(groups,length)
# groupStart=c(0,cumsum(groupLengths[1:(length(groups)-1)]))
# groupEnd=cumsum(groupLengths)-1
# groupColumns=unlist(groups)-1
# beta=rnorm(length(groupColumns)*k)
# penaltyFactor=rep(1,2)
# sampleWeights=rep(1,n)
# response=X%*%beta+rnorm(n)
# 
# # Eigen::VectorXd constructorTest(int k_, Eigen::VectorXi groupStart_, Eigen::VectorXi groupEnd_
# #                                 ,Eigen::VectorXi groupColumns_, Eigen::VectorXd beta_, Eigen::VectorXd penaltyFactor_
# #                                 , Eigen::VectorXd response_, Eigen::VectorXd sampleWeights_, Eigen::VectorXd residuals_
# #                                 ,Eigen::MatrixXd X_) 
# # 
# #RinvList=constructorTest(k,groupStart,groupEnd,groupColumns,rep(0,15),penaltyFactor,response,sampleWeights,X)
# #RinvListSparse=constructorTestSparse(k,groupStart,groupEnd,groupColumns,beta,penaltyFactor,response,Matrix::Matrix(X,sparse=TRUE))
# # # 
# # # RinvList2=list()
# # # for (i in 1:length(groups)){
# # #   RinvList2[[i]]=solve(qr.R(qr(X[,groups[[i]]])))
# # # }
# # 
# # q=qr.Q(qr(X))
# # 
# # Q=NULL
# # for (i in 1:6){
# #   Q=cbind(Q,qr.Q(qr(X[,groups[[i]]])))
# # }
# # #beta[c(46:54)]=0
# # #beta[c(1:45,55:72)]=0
# # #beta[46:72]=0
# # #beta[1:55]=0
# # # 
# # # kronecker(Q,diag(3))%*%beta-RinvList
# # # Q2=kronecker(Q,diag(3))
# # # microbenchmark(Q2%*%beta)
# # # Q3=t(Q2)
# # # resp3=rep(response,3)
# # # microbenchmark(Q3%*%resp3)
# # # microbenchmark(constructorTest(k,groupStart,groupEnd,groupColumns,beta,penaltyFactor,response,sampleWeights,residuals,X)
# # # )
# # # microbenchmark(constructorTestSparse(k,groupStart,groupEnd,groupColumns,beta,penaltyFactor,response,sampleWeights,residuals,Matrix::Matrix(X,sparse=TRUE))
# # # )
# # # microbenchmark(1*X)
library(bcd)
library(testthat)
library(glmnet)
library(grpreg)
data("exampleData")

test_that("Linear regression", {
  fitLinear=fit_bcd(X=X,y=y_gaussian,family="gaussian",groups=grouping,penaltyFactor=penaltyFactor)
  fitGrpreg=grpreg::grpreg(X=X[,-1],y=y_gaussian,family="gaussian",group=c(rep(1,9),rep(2,20),rep(3,20)),eps = 10^-16)
  betaBCD=matrix(unlist(fitLinear$beta),ncol=length(fitLinear$beta))
  expect_equal(betaBCD,matrix(fitGrpreg$beta,ncol=100),tol=10^-8)
})

test_that("Logistic regression vs. grpreg", {
  fitLogistic=fit_bcd(X=X,y=y_binary,family="logistic",groups=grouping,penaltyFactor=penaltyFactor)
  fitGrpreg=grpreg(X=X[,-1],y=y_binary,family="binomial",group=c(rep(1,9),rep(2,20),rep(3,20)),eps = 10^-15)
  betaBCD=matrix(0,length(fitLogistic$beta[[1]][,1]),ncol(fitGrpreg$beta))
  for (i in 1:ncol(fitGrpreg$beta)){
    betaBCD[,i]=fitLogistic$beta[[i]][,2]-fitLogistic$beta[[i]][,1]
  }
  expect_equal(betaBCD[,1:ncol(fitGrpreg$beta)],matrix(fitGrpreg$beta,ncol=ncol(fitGrpreg$beta)),tol=10^-8)
})

test_that("Poisson regression vs. glmnet", {
  #grpreg seems to be off by a bit
  fitPoisson=fit_bcd(X=X,y=y_count,family="poisson",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)))
  fitGlmnet=glmnet(x=X[,-1],y=y_count,family="poisson",thresh=10^-30)
  fitGlmnet=glmnet(x=X[,-1],y=y_count,family="poisson",lambda=fitPoisson$lambda*fitGlmnet$lambda[1]/fitPoisson$lambda[1],thresh=10^-30)
  for (i in 1:length(fitPoisson$beta)){
    b1=fitPoisson$beta[[i]]
    b2=coefficients(fitGlmnet)[,i]
    expect_equal(matrix(b1,ncol=1),matrix(b2,ncol=1),tol=10^-10)
  }
})

test_that("Multiresponse linear vs. glmnet", {
  fitMultiresponse=fit_bcd(X=X,y=y_multiresponse,family="gaussian",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)))
  fitGlmnet=glmnet(x=X[,-1],y=y_multiresponse,family="mgaussian",thresh=10^-30)
  fitGlmnet=glmnet(x=X[,-1],y=y_multiresponse,family="mgaussian",lambda=fitMultiresponse$lambda*fitGlmnet$lambda[1]/fitMultiresponse$lambda[1],thresh=10^-30)
  for (i in 1:length(fitMultiresponse$beta)){
    b1=fitMultiresponse$beta[[i]]
    b2=NULL
    b2temp=coefficients(fitGlmnet)
    for (j in 1:3){
      b2=cbind(b2,b2temp[[j]][,i])
    }
    expect_equal(matrix(b1,ncol=ncol(b1)),matrix(b2,ncol=ncol(b2)),tol=10^-10)
  }
})

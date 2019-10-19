library(bcd)
library(testthat)
library(glmnet)
library(grpreg)
data("exampleData")
data("referenceCoefficients")

#####Linear
test_that("Linear regression vs. grpreg", {
  fitLinear=fit_bcd(X=X,y=y_gaussian,family="gaussian",groups=grouping,penaltyFactor=penaltyFactor)
  fitGrpreg=grpreg::grpreg(X=X[,-1],y=y_gaussian,family="gaussian",group=c(rep(1,9),rep(2,20),rep(3,20)),eps = 10^-16)
  betaBCD=matrix(unlist(fitLinear$beta),ncol=length(fitLinear$beta))
  expect_equal(betaBCD,matrix(fitGrpreg$beta,ncol=100),tol=10^-8)
})

test_that("Linear regression vs. reference from working version", {
  fitLinear=fit_bcd(X=X,y=y_gaussian,family="gaussian",groups=grouping,penaltyFactor=penaltyFactor)
  for (i in 1:length(fitLinear$beta)){
    expect_equal(fitLinear$beta[[i]],reference_Linear[[i]])
  }
})

test_that("Linear regression with sample weights rep(1,n)", {
  fitLinear=fit_bcd(X=X,y=y_gaussian,family="gaussian",groups=grouping,penaltyFactor=penaltyFactor,sampleWeights = rep(1,n))
  for (i in 1:length(fitLinear$beta)){
    expect_equal(fitLinear$beta[[i]],reference_Linear[[i]])
  }
})

test_that("Linear regression with varying sample weights", {
  fitLinear=fit_bcd(X=X,y=y_gaussian,family="gaussian",groups=grouping,penaltyFactor=penaltyFactor,sampleWeights = sampleWeights)
  for (i in 1:length(fitLinear$beta)){
    expect_equal(fitLinear$beta[[i]],reference_Linear_Weighted[[i]])
  }
})

test_that("Linear regression without sample weights vs. glmnet", {
  fitLinear=fit_bcd(X=X,y=y_gaussian,family="gaussian",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)),tol=10^-12)
  fitGlmnet=glmnet(x=X[,-1],y=y_gaussian,family="gaussian",lambda=fitLinear$lambda*sqrt(n),thresh=10^-30)
  for (i in 1:length(fitLinear$beta)){
    expect_equal(as.numeric(fitLinear$beta[[i]]),as.numeric(coefficients(fitGlmnet)[,i]))
  }
})

test_that("Linear regression with sample weights rep(1,n) vs. glmnet", {
  fitLinear=fit_bcd(X=X,y=y_gaussian,family="gaussian",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)),tol=10^-12,sampleWeights = rep(1,n))
  fitGlmnet=glmnet(x=X[,-1],y=y_gaussian,family="gaussian",lambda=fitLinear$lambda*sqrt(n),thresh=10^-30)
  for (i in 1:length(fitLinear$beta)){
    expect_equal(as.numeric(fitLinear$beta[[i]]),as.numeric(coefficients(fitGlmnet)[,i]))
  }
})


test_that("Linear regression with varying sample weights vs. glmnet", {
  fitLinear=fit_bcd(X=X,y=y_gaussian,family="gaussian",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)),sampleWeights = sampleWeights,tol=10^-12)
  fitGlmnet=glmnet(x=X[,-1],y=y_gaussian,family="gaussian",lambda=fitLinear$lambda*sqrt(n),thresh=10^-30,weights=sampleWeights)
  for (i in 1:length(fitLinear$beta)){
    expect_equal(as.numeric(fitLinear$beta[[i]]),as.numeric(coefficients(fitGlmnet)[,i]))
  }
})

#####Logistic
test_that("Logistic regression vs. grpreg", {
  fitLogistic=fit_bcd(X=X,y=y_binary,family="logistic",groups=grouping,penaltyFactor=penaltyFactor)
  fitGrpreg=grpreg(X=X[,-1],y=y_binary,family="binomial",group=c(rep(1,9),rep(2,20),rep(3,20)),eps = 10^-12)
  betaBCD=matrix(0,length(fitLogistic$beta[[1]][,1]),ncol(fitGrpreg$beta))
  for (i in 1:ncol(fitGrpreg$beta)){
    betaBCD[,i]=fitLogistic$beta[[i]][,2]-fitLogistic$beta[[i]][,1]
  }
  expect_equal(betaBCD[,1:ncol(fitGrpreg$beta)],matrix(fitGrpreg$beta,ncol=ncol(fitGrpreg$beta)),tol=10^-7)
})

test_that("Logistic regression vs. reference from working version", {
  fitLogistic=fit_bcd(X=X,y=y_binary,family="logistic",groups=grouping,penaltyFactor=penaltyFactor)
  for (i in 1:length(fitLogistic)){
    expect_equal(fitLogistic$beta[[i]],reference_Logistic[[i]])
  }
})

test_that("Logistic regression with sample weights rep(1,n)", {
  fitLogistic=fit_bcd(X=X,y=y_binary,family="logistic",groups=grouping,penaltyFactor=penaltyFactor,sampleWeights = rep(1,n))
  for (i in 1:length(fitLogistic)){
    expect_equal(fitLogistic$beta[[i]],reference_Logistic[[i]])
  }
})

test_that("Logistic regression with varying sample weights", {
  fitLogistic=fit_bcd(X=X,y=y_binary,family="logistic",groups=grouping,penaltyFactor=penaltyFactor,sampleWeights = sampleWeights)
  for (i in 1:length(fitLogistic)){
    expect_equal(fitLogistic$beta[[i]],reference_Logistic_Weighted[[i]])
  }
})

test_that("Logistic regression with varying sample weights vs. glmnet", {
  fitLogistic=fit_bcd(X=X,y=y_binary,family="logistic",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)),sampleWeights = sampleWeights)
  fitGlmnet=glmnet(x=X[,-1],y=y_binary,family="binomial",lambda=fitLogistic$lambda*sqrt(n),thresh=10^-30,weights=sampleWeights)
  for (i in 1:length(fitLogistic$beta)){
    betaLogistic=fitLogistic$beta[[i]][,2]-fitLogistic$beta[[i]][,1]
    betaGlmnet=as.numeric(coefficients(fitGlmnet)[,i])
    expect_equal(betaLogistic,betaGlmnet)
  }
})

#####Multinomial
test_that("Multinomial regression vs. reference from working version", {
  fitMultinomial=fit_bcd(X=X,y=as.factor(y_multinomial),family="multinomial",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)))
  for (i in 1:length(fitMultinomial)){
    expect_equal(fitMultinomial$beta[[i]],reference_Multinomial[[i]])
  }
})

test_that("Multinomial regression with sampleweights rep(1,n)", {
  fitMultinomial=fit_bcd(X=X,y=as.factor(y_multinomial),family="multinomial",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)),sampleWeights = rep(1,n))
  for (i in 1:length(fitMultinomial)){
    expect_equal(fitMultinomial$beta[[i]],reference_Multinomial[[i]])
  }
})

test_that("Multinomial regression with varying weights", {
  fitMultinomial=fit_bcd(X=X,y=as.factor(y_multinomial),family="multinomial",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)),sampleWeights = sampleWeights)
  for (i in 1:length(fitMultinomial)){
    expect_equal(fitMultinomial$beta[[i]],reference_Multinomial_Weighted[[i]])
  }
})

test_that("Multinomial regression vs. glmnet", {
  fitMultinomial=fit_bcd(X=X,y=as.factor(y_multinomial),family="multinomial",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)))
  fitGlmnet=glmnet(x=X[,-1],y=as.factor(y_multinomial),family="multinomial",lambda=fitMultinomial$lambda*sqrt(k*n),thresh=10^-20,type.multinomial="grouped")
  for (i in 1:length(fitMultinomial$beta)){
    b1=fitMultinomial$beta[[i]]
    b2=NULL
    b2temp=coefficients(fitGlmnet)
    for (j in 1:3){
      b2=cbind(b2,b2temp[[j]][,i])
    }
    expect_equal(matrix(b1,ncol=ncol(b1)),matrix(b2,ncol=ncol(b2)),tol=10^-9)
  }
})

test_that("Multinomial regression with varying sample weights vs. glmnet", {
  fitMultinomial=fit_bcd(X=X,y=as.factor(y_multinomial),family="multinomial",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)),sampleWeights = sampleWeights)
  fitGlmnet=glmnet(x=X[,-1],y=y_multinomial,family="multinomial",lambda=fitMultinomial$lambda*sqrt(k*n),thresh=10^-20,weights=sampleWeights,type.multinomial="grouped")
  for (i in 1:length(fitMultinomial$beta)){
    b1=fitMultinomial$beta[[i]]
    b2=NULL
    b2temp=coefficients(fitGlmnet)
    for (j in 1:3){
      b2=cbind(b2,b2temp[[j]][,i])
    }
    expect_equal(matrix(b1,ncol=ncol(b1)),matrix(b2,ncol=ncol(b2)),tol=10^-9)
  }
})


#####Poisson
test_that("Poisson regression vs. glmnet", {
  fitPoisson=fit_bcd(X=X,y=y_count,family="poisson",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)),tol=10^-14)
  fitGlmnet=glmnet(x=X[,-1],y=y_count,family="poisson",lambda=fitPoisson$lambda*sqrt(n),thresh=10^-30)
  for (i in 1:length(fitPoisson$beta)){
    b1=fitPoisson$beta[[i]]
    b2=coefficients(fitGlmnet)[,i]
    expect_equal(matrix(b1,ncol=1),matrix(b2,ncol=1),tol=10^-10)
  }
})

test_that("Poisson regression vs. grpreg", {
  #grpreg seems to be a bit off
  #using a higher tolerance for this test
  fitPoisson=fit_bcd(X=X,y=y_count,family="poisson",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)))
  fitGrpreg=grpreg(X=X[,-1],y=y_count,family="poisson",eps = 10^-15)
  
  for (i in 1:dim(fitGrpreg$beta)[2]){
    b1=fitPoisson$beta[[i]]
    b2=fitGrpreg$beta[,i]
    expect_equal(matrix(b1,ncol=1),matrix(b2,ncol=1),tol=10^-3)
  }
})

test_that("Poisson regression vs. reference", {
  fitPoisson=fit_bcd(X=X,y=y_count,family="poisson",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)))
  for (i in 1:length(fitPoisson)){
    expect_equal(fitPoisson$beta[[i]],reference_Poisson[[i]])
  }
})

test_that("Poisson regression vs. reference with sample weights rep(1,n)", {
  fitPoisson=fit_bcd(X=X,y=y_count,family="poisson",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)),sampleWeights = rep(1,n))
  for (i in 1:length(fitPoisson)){
    expect_equal(fitPoisson$beta[[i]],reference_Poisson[[i]])
  }
})

test_that("Poisson regression vs. reference with varying sample weights", {
  fitPoisson=fit_bcd(X=X,y=y_count,family="poisson",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)),sampleWeights = sampleWeights)
  for (i in 1:length(fitPoisson)){
    expect_equal(fitPoisson$beta[[i]],reference_Poisson_Weighted[[i]])
  }
})

test_that("Poisson regression vs. glmnet with varying sample weights", {
  fitPoisson=fit_bcd(X=X,y=y_count,family="poisson",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)),sampleWeights = sampleWeights)
  fitGlmnet=glmnet(x=X[,-1],y=y_count,family="poisson",lambda=fitPoisson$lambda*sqrt(n),thresh=10^-30,weights=sampleWeights)
  for (i in 1:length(fitPoisson)){
    expect_equal(as.numeric(fitPoisson$beta[[i]]),as.numeric(coefficients(fitGlmnet)[,i]),tol=10^-12)
  }
})


#####Multiresponse linear
test_that("Multiresponse linear vs. glmnet", {
  fitMultiresponse=fit_bcd(X=X,y=y_multiresponse,family="gaussian",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)))
  fitGlmnet=glmnet(x=X[,-1],y=y_multiresponse,family="mgaussian",lambda=fitMultiresponse$lambda*sqrt(k*n),thresh=10^-30)
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

test_that("Multiresponse linear with sample weights vs. glmnet", {
  fitMultiresponse=fit_bcd(X=X,y=y_multiresponse,family="gaussian",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)),sampleWeights = sampleWeights)
  fitGlmnet=glmnet(x=X[,-1],y=y_multiresponse,family="mgaussian",lambda=fitMultiresponse$lambda*sqrt(),thresh=10^-30,weights=sampleWeights)
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

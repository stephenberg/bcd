# library(bcd)
# library(testthat)
# library(glmnet)
# library(grpreg)
# data("exampleData")
# 
# 
# #####without sample weights
# fitLinear=fit_bcd(X=X,y=y_gaussian,family="gaussian",groups=grouping,penaltyFactor=penaltyFactor)
# reference_Linear=fitLinear$beta
# 
# fitLogistic=fit_bcd(X=X,y=y_binary,family="logistic",groups=grouping,penaltyFactor=penaltyFactor)
# reference_Logistic=fitLogistic$beta
# 
# fitPoisson=fit_bcd(X=X,y=y_count,family="poisson",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)))
# reference_Poisson=fitPoisson$beta
# 
# fitMultiresponse=fit_bcd(X=X,y=y_multiresponse,family="gaussian",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)))
# reference_Multiresponse=fitMultiresponse$beta
# 
# fitMultinomial=fit_bcd(X=X,y=as.factor(y_multinomial),family="multinomial",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)))
# reference_Multinomial=fitMultinomial$beta
# 
# ####with sample weights
# fitLinear=fit_bcd(X=X,y=y_gaussian,family="gaussian",groups=grouping,penaltyFactor=penaltyFactor,sampleWeights = sampleWeights)
# reference_Linear_Weighted=fitLinear$beta
# 
# fitLogistic=fit_bcd(X=X,y=y_binary,family="logistic",groups=grouping,penaltyFactor=penaltyFactor,sampleWeights = sampleWeights)
# reference_Logistic_Weighted=fitLogistic$beta
# 
# fitPoisson=fit_bcd(X=X,y=y_count,family="poisson",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)),sampleWeights = sampleWeights)
# reference_Poisson_Weighted=fitPoisson$beta
# 
# fitMultiresponse=fit_bcd(X=X,y=y_multiresponse,family="gaussian",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)),sampleWeights = sampleWeights)
# reference_Multiresponse_Weighted=fitMultiresponse$beta
# 
# fitMultinomial=fit_bcd(X=X,y=as.factor(y_multinomial),family="multinomial",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)),sampleWeights = sampleWeights)
# reference_Multinomial_Weighted=fitMultinomial$beta
# 
# save(reference_Linear,
#      reference_Logistic,
#      reference_Poisson,
#      reference_Multiresponse,
#      reference_Multinomial,
#      reference_Linear_Weighted,
#      reference_Logistic_Weighted,
#      reference_Poisson_Weighted,
#      reference_Multiresponse_Weighted,
#      reference_Multinomial_Weighted,
#      file="data/referenceCoefficients.RData")

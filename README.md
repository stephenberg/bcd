## bcd: block coordinate descent for group lasso problems

### Description

The **bcd** package uses block coordinate descent to fit common GLMs with a group lasso penalty. The model fitting implementation is in C++, and the package is used in R with wrapper functions from Rcpp. Run help("fit_bcd") in R for descriptions of the function arguments and for some examples of using the package.

Features:

* tested against **glmnet** and **grpreg**
* fits linear, logistic, multinomial, multiresponse Gaussian, and Poisson models
* supports sparse and dense design matrices
* supports multiple group penalization options:
  * sparse group lasso penalization
  * overlapping group lasso penalization
  * models with 1 or more than 1 unpenalized coefficient
  * models with 0 unpenalized coefficients, including models where the intercept is penalized or the model does not contain an intercept
  * varying weights for each sample
  
### Installation

Install using **devtools** package:

```
devtools::install_github("stephenberg/bcd")
```

### Example

Linear regression with simulated data

```
library(bcd)
data(exampleData)
fit_linear=fit_bcd(X=X,y=y_gaussian,family="gaussian",groups=grouping,penaltyFactor=penaltyFactor)
```

### Example with overlapping groups

A more complicated logistic regression example with overlapping groups:

```
data(exampleData)
grouping=list(as.integer(1),2:10,11:30,11:50)
fit_overlap=fit_bcd(X=X,y=y_binary,family="logistic",groups=grouping,penaltyFactor=penaltyFactor)
```

### Comparison to results from glmnet and grpreg

The tests/ folder of the repository contains a substantial collection of test functions for comparing **bcd**, **glmnet**, and **grpreg** under many settings. Below, we show the setting and output for two of the linear regression test examples.

The tests use the example data that comes with the package:

```
library(bcd)
data(exampleData)
```

Testing linear regression vs. **grpreg**, with 3 groups and an unpenalized intercept:

```
fitBCD=fit_bcd(X=X,y=y_gaussian,family="gaussian",groups=grouping,penaltyFactor=penaltyFactor)
fitGrpreg=grpreg::grpreg(X=X[,-1],y=y_gaussian,family="gaussian",group=c(rep(1,9),rep(2,20),rep(3,20)),eps = 10^-16)
betaBCD=matrix(unlist(fitBCD$beta),ncol=length(fitBCD$beta))
max(abs(betaBCD-fitGrpreg$beta))
```
[1] 4.468648e-15

Testing linear regression vs. **glmnet**, with lasso (*l*<sub>1</sub>) penalty only

```
fitBCD=fit_bcd(X=X,y=y_gaussian,family="gaussian",groups=as.list(1:50),penaltyFactor=c(0,rep(1,49)),tol=10^-12)
fitGlmnet=glmnet(x=X[,-1],y=y_gaussian,family="gaussian",lambda=fitBCD$lambda*sqrt(n),thresh=10^-30)
betaBCD=matrix(unlist(fitBCD$beta),ncol=length(fitBCD$beta))
max(abs(betaBCD-coefficients(fitGlmnet)))
```
[1] 2.250804e-15

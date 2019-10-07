#' Fit model by block coordinate descent
#'
#' Fits a model by block coordinate descent
#'
#' @param X design matrix
#' @param y response vector (or matrix)
#' @param k number of response categories (multinomial) or dimension of response matrix (Gaussian)
#' @param family type of regression: one of "gaussian", "multinomial", "poisson", or "logistic"
#' @param groups list containing the group index of each column
#' @param penaltyFactor vector containing penalty level for each group
#' @param lambda penalty vector for model fitting
#' @param lambdaMinRatio value to use for lambdaMin/lambdaMax (defaults to 0.0001)
#' @param nLambda number of lambda values to use
#' @param sampleWeights different weights for each sample
#' @param maxit maximum number of coordinate descent iterations (default 100000)
#' @param tolerance solution tolerance (default 10^-12)
#' @param offset offset to use for Poisson regression
#' @param eigenValueTolerance tolerance for deciding whether columns in a group are linearly independent (default 10^-9)
#' @param rescale force columns to have norm 1 after orthogonalization (defaults to TRUE)
#' @param useDevTol use deviance tolerance
#' @param devTol tolerance for deviance 
#' 
#' @return model fit
#' 
#' @examples
#' #no example
#' 
#'
#' @export
bcdFit<-function(X,
                 y,
                 k=1,
                 family = c("gaussian", "multinomial", "poisson", "logistic"),
                 groups,
                 penaltyFactor,
                 lambda=NULL,
                 lambdaMinRatio=0.0001,
                 nLambda=100,
                 sampleWeights=NULL,
                 maxit=100000,
                 tolerance=10^-12,
                 offset=NULL,
                 eigenValueTolerance=10^-9,
                 rescale=TRUE,
                 useDevTol=TRUE,
                 devTol=0.9
){

  family <- match.arg(family)
  
  if ((!is.list(groups)) | (length(groups)<1)){
    stop("groups must be a list containing vectors of integers")
  }
  if (!is.integer(unlist(groups))){
    stop("the entries in each group must be integers between 1 and the number of columns in X")
  }

  sparse=FALSE
  if (class(X)!="matrix" & class(X)!="dgeMatrix"){
    if (class(X)=="dgCMatrix"){
      sparse=TRUE
    }
    if (class(X)!="dgCMatrix"){
      stop("X must be of type matrix, dgCMatrix, or dgeMatrix")
    }
  }

  nGroups=length(groups)
  groupLengths=sapply(groups,length)
  groupStart=c(0,cumsum(groupLengths[1:(length(groups)-1)]))
  groupEnd=cumsum(groupLengths)-1
  groupColumns=unlist(groups)-1

  useLambda=FALSE

  if (!is.null(lambda)){
    useLambda=TRUE
    lambda=lambda
    if (nLambda<1){
      stop("the lambda path must have length at least 1.")
    }
    if (any(lambda<=0)){
      stop("All lambda values must be (strictly) >0")
    }
  }
  if (is.null(lambda)){
    lambda=c(-1)
  }

  useWeights=FALSE
  if (!is.null(sampleWeights)){
    useWeights=TRUE
    if (!is.vector(sampleWeights)){
      stop("sampleWeights, if supplied, must be a vector with length the number of rows in X")
    }
    if (any(sampleWeights<0)){
      stop("all entries in sampleWeights must be >=0.")
    }
  }
  if (is.null(sampleWeights)){
    sampleWeights=c(1)
  }


  if (!is.vector(penaltyFactor)){
    stop("a penalty factor vector must be supplied with length nGroups. The code will internally multiply
         the entries by sqrt(group size) for each group.")
  }
  if (!is.numeric(penaltyFactor)){
    stop("a numeric penalty factor vector must be supplied with length nGroups. The code will internally multiply
         the entries by sqrt(group size) for each group.")
  }
  if (length(penaltyFactor)!=nGroups){
      stop("a penalty factor vector must be supplied with length nGroups. The code will internally multiply
           the entries by sqrt(group size) for each group.")
  }
  if (any(penaltyFactor<0)){
    stop("all entries in penaltyFactor must be >=0.")
  }
  unpenalized=which(penaltyFactor==0)
  if (length(unpenalized)>1){
    stop("at most one group may be unpenalized. Please combine all unpenalized groups into a single group.")
  }


  if (family=="gaussian"){
    if (is.vector(y)){
      if (length(y)==dim(X)[1]){
        y=matrix(y)
      }
      if (length(y)!=dim(X)[1]){
        stop("y must have the same length as the number of rows in X")
      }
    }
    if (is.matrix(y)){
      if (dim(y)[1]!=dim(X)[1]){
        stop("y must have the same number of rows as X")
      }
      if (dim(y)[2]!=k){
        stop("for the multiresponse Gaussian, y must have k columns")
      }
    }
    result=NULL
    if (!sparse){
      result=(multiResponseGaussianDense(nGroups,
                                         groupStart,
                                         groupEnd,
                                         groupColumns,
                                         tolerance,
                                         penaltyFactor,
                                         maxit,
                                         y,
                                         X,
                                         1,
                                         k,
                                         nLambda,
                                         lambdaMinRatio,
                                         eigenValueTolerance,
                                         rescale,
                                         useLambda,
                                         lambda,
                                         useWeights,
                                         sampleWeights,
                                         useDevTol,
                                         devTol))
    }
    if (sparse){
      result=(multiResponseGaussianSparse(nGroups,
                                          groupStart,
                                          groupEnd,
                                          groupColumns,
                                          tolerance,
                                          penaltyFactor,
                                          maxit,
                                          y,
                                          X,
                                          1,
                                          k,
                                          nLambda,
                                          lambdaMinRatio,
                                          eigenValueTolerance,
                                          rescale,
                                          useLambda,
                                          lambda,
                                          useWeights,
                                          sampleWeights,
                                          useDevTol,
                                          devTol))
    }
  }
  if (family=="multinomial" | family=="logistic"){
    response=NULL
    if (is.factor(y)){
      if (length(unique(y))!=k){
        stop("number of levels in y must equal k.")
      }
      if (length(y)!=(dim(X))[1]){
        stop("length of y must equal number of rows in X")
      }
      response=matrix(0,dim(X)[1],k)
      whichNonZero=as.numeric(y)
      for (i in 1:dim(X)[1]){
        response[i,whichNonZero[i]]=1
      }
      if (any(apply(response,2,sum)==0)){
        stop("all response categories must have >0 observations.")
      }
    }
    if (is.vector(y)){
      if (k!=2){
        stop("if the response is given as a non-factor-type vector, k must equal 2.")
      }
      if (any(y<0)|any(y>1)){
        stop("all responses must be between 0 and 1")
      }
      response=c(1-y,y)
    }
    if (is.matrix(y)){
      if (dim(y)[2]!=k){
        stop("the dimension of y not equal to k")
      }
      if (dim(y)[1]!=dim(X)[1]){
        stop("the number of rows in y must equal the number of rows in X.")
      }
      if (any(y<0)|any(y>1)){
        stop("all responses must be between 0 and 1")
      }
      if (any(abs(apply(y,1,sum)-1)>10^-10)){
        stop("all rows in the response must add to 1.")
      }
      if (any(apply(y,2,sum)==0)){
        stop("all response categories must have >0 observations.")
      }
      response=y
    }
    if (!sparse){
      result=(multinomialDense(nGroups,
                               groupStart,
                               groupEnd,
                               groupColumns,
                               tolerance,
                               penaltyFactor,
                               maxit,
                               response,
                               X,
                               k,
                               nLambda,
                               lambdaMinRatio,
                               eigenValueTolerance,
                               rescale,
                               useLambda,
                               lambda,
                               useWeights,
                               sampleWeights,
                               useDevTol,
                               devTol))
    }
    else{
      result=(multinomialSparse(nGroups,
                                groupStart,
                                groupEnd,
                                groupColumns,
                                tolerance,
                                penaltyFactor,
                                maxit,
                                response,
                                X,
                                k,
                                nLambda,
                                lambdaMinRatio,
                                eigenValueTolerance,
                                rescale,
                                useLambda,
                                lambda,
                                useWeights,
                                sampleWeights,
                                useDevTol,
                                devTol))
    }
  }
  if (family=="poisson"){
    if (k!=1){
      stop("k must equal 1 for poisson regression")
    }
    if (!is.vector(y)){
      if (!is.matrix(y)){
        stop("y must be a length n vector or a n-by-1 matrix.")
      }
    }
    if (!is.numeric(y)){
      stop("y must be numeric.")
    }
    if (any(y)<0){
      stop("Responses must be positive for Poisson regression.")
    }
    response=matrix(y,ncol=1)
    if (dim(response)[1]!=dim(X)[1]){
      stop("number of responses must be equal to the number of rows in X.")
    }
    if (is.null(offset)){
      offset=rep(0,dim(X)[1])
    }
    if (!is.null(offset)){
      if (!is.vector(offset)){
        stop("offset must be a numeric vector of length n.")
      }
      if (!is.numeric(offset)){
        stop("offset must be a numeric vector of length n.")
      }
    }
    offset=matrix(offset,ncol=1)
    if (!sparse){
      result=(poissonDense(nGroups,
                           groupStart,
                           groupEnd,
                           groupColumns,
                           tolerance,
                           penaltyFactor,
                           maxit,
                           response,
                           X,
                           nLambda,
                           lambdaMinRatio,
                           eigenValueTolerance,
                           rescale,
                           useLambda,
                           lambda,
                           useWeights,
                           sampleWeights,
                           offset,
                           useDevTol,
                           devTol))
    }
    if (sparse){
      result=(poissonSparse(nGroups,
                            groupStart,
                            groupEnd,
                            groupColumns,
                            tolerance,
                            penaltyFactor,
                            maxit,
                            response,
                            X,
                            nLambda,
                            lambdaMinRatio,
                            eigenValueTolerance,
                            rescale,
                            useLambda,
                            lambda,
                            useWeights,
                            sampleWeights,
                            offset,
                            useDevTol,
                            devTol))
    }
  }
  if (!(family%in%c("gaussian","logistic","multinomial","poisson"))){
    stop("family must be one of gaussian, logistic, multinomial, or poisson.")
  }
  return(result)
  }
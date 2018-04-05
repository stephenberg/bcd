multinomialCV<-function(X,response,k,nFolds=10,folds,groups,lambda,penaltyFactor,tolerance=10^-9,useDevTol=FALSE){
  betaList=list()
  devMat=matrix(0,length(lambda),nFolds)
  for (foldInd in 1:nFolds){
    print(foldInd)
    Xholdout=X[-folds[[foldInd]],]
    responseHoldout=response[-folds[[foldInd]],]
    betaList[[foldInd]]=bcdFit(X=Xholdout,y=responseHoldout,k=dim(response)[2],family = "multinomial",groups=groups,penaltyFactor = penaltyFactor,lambda=lambda,tolerance=tolerance,useDevTol = useDevTol)$beta
    for (lamInd in 1:length(lambda)){
      devMat[lamInd,foldInd]=(multinomialDeviance(X[folds[[foldInd]],],response[folds[[foldInd]],],betaList[[foldInd]][[lamInd]]))
    }
  }
  return(devMat)
}

multinomialDeviance<-function(X,response,beta){
  n=dim(X)[1]
  k=dim(response)[2]
  linPred=X%*%beta
  deviance=0
  for (nInd in 1:n){
    linPred[nInd,]=linPred[nInd,]-max(linPred[nInd,])
    mu=exp(linPred[nInd,])
    mu=mu/sum(mu)
    deviance=deviance+sum(response[nInd,]*log(mu))
  }
  return(-2*deviance)
}
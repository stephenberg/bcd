# multinomial<-function(X
#                       ,y
#                       ,k=-1
#                       ,grouping=list()
#                       ,penaltyFactor=numeric(0)
#                       ,innerTolerance=10^-9
#                       ,innerMaxit=10000
#                       ,outerTolerance=10^-6
#                       ,outerMaxit=10000
#                       ,eigenValueTolerance=10^-9
#                       ,lambdaMinRatio=10^-4
#                       ,lambda=NULL
#                       ,nLambda=100
#                       ,scale=TRUE){
#   
#   #check matrix type
#   matrixType=""
#   if (class(X)=="matrix" | class(X)=="dgeMatrix"){
#     matrixType="matrix"
#   }else if (class(X)=="dgCMatrix"){
#     matrixType="sparse matrix"
#   }else{
#       stop("Error: X must be a sparse or dense matrix")
#   }
#     
#   if (k<=0){
#     stop("Error: please supply a positive number k of classification categories")
#   }
#   n=dim(X)[1]
#   response=NULL
#   whichCategory=NULL
#   if (class(y)=="factor"){
#     whichCategory=as.numeric(y)
#     if (k!=sum(table(y)>0)){
#       stop(paste("Error: k=",k," but there are ",sum(table(y)>0)," levels with count>0 in the supplied response"))
#     }
#     response=matrix(0,n,k)
#     for (i in 1:n){
#       response[i,whichCategory(i)]=1
#     }
#     response=c(t(response))
#   }else{
#     if (!is.numeric(y)){
#       stop("Error: response must be either a factor or a length (nk) vector of doubles between 0 and 1")
#     }
#     if (length(y)!=(n*k)){
#       stop("Error: response must be either a length (n) factor or a length (nk) vector of doubles between 0 and 1")
#     }
#     if (min(y)<0 | max(y)>1){
#       stop("Error: response must be either a length(n) factor or a length (nk) vector of doubles between 0 and 1")
#     }
#     y=matrix(y,ncol=n)
#     if (!all(apply(y,2,sum)==1)){
#       stop("Error: responses must add to 1 for each observation")
#     }
#     response=c(y)
#   }
#   if (length(grouping)==0){
#     stop("Error: no grouping supplied")
#   }
#   if ((length(penaltyFactor)!=0) & (length(penaltyFactor)!=length(grouping))){
#     stop("Error: penaltyFactor must either be the same length as the grouping vector or else have length 0")
#   }
#   if (length(penaltyFactor)==0){
#     print("Penalty factor vector not supplied-all groups will be penalized")
#     penaltyFactor=rep(1,length(grouping))
#   }
#   if (sum(penaltyFactor==0)>1){
#     stop("Error: only 1 unpenalized group is allowed: please move all unpenalized columns to this group")
#   }
#   
#   givenLambda=FALSE;
#   if (!is.null(lambda)){
#     givenLambda=TRUE
#   }else{
#     lambda=1
#   }
#   
#   #compute groupings
#   groupLengths=sapply(grouping,length)
#   groupStart=c(0,cumsum(groupLengths[1:(length(grouping)-1)]))
#   groupEnd=cumsum(groupLengths)-1
#   groupColumns=unlist(grouping)-1
#   nGroups=length(groupStart)
#   
#   if (matrixType=="matrix"){
#     # Eigen::MatrixXd X_
#     # ,Eigen::VectorXd response_
#     # ,Eigen::VectorXd penaltyFactor_
#     # ,double innerTolerance_
#     # ,int innerIterationsMax_
#     # ,double outerTolerance_
#     # ,int outerIterationsMax_
#     # ,double eigenValueTolerance_
#     # ,double lambdaMinRatio_
#     # ,int nLambda_
#     # ,bool givenLambda_
#     # ,Eigen::VectorXd lambda_
#     # ,bool scale_
#     # ,Eigen::VectorXi groupStart_
#     # ,Eigen::VectorXi groupEnd_
#     # ,Eigen::VectorXi groupColumns_
#     # ,int nGroups_
#     # ,int k_
#     return(multinomialFitDense(X
#                                 ,response
#                                 ,penaltyFactor
#                                 ,innerTolerance
#                                 ,innerMaxit
#                                 ,outerTolerance
#                                 ,outerMaxit
#                                 ,eigenValueTolerance
#                                 ,lambdaMinRatio
#                                 ,nLambda
#                                 ,givenLambda
#                                 ,lambda
#                                 ,scale
#                                 ,groupStart
#                                 ,groupEnd
#                                 ,groupColumns
#                                 ,nGroups
#                                 ,k))
#   }
#   else{
#     return(multinomialFitSparse(X
#                                ,response
#                                ,penaltyFactor
#                                ,innerTolerance
#                                ,innerMaxit
#                                ,outerTolerance
#                                ,outerMaxit
#                                ,eigenValueTolerance
#                                ,lambdaMinRatio
#                                ,nLambda
#                                ,givenLambda
#                                ,lambda
#                                ,scale
#                                ,groupStart
#                                ,groupEnd
#                                ,groupColumns
#                                ,nGroups
#                                ,k))
#   }
# }
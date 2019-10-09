# ###script to generate the example data
# ###uncomment to run
# 
# 
# set.seed(33)
# 
# n=1000
# p=50
# k=3
# X=matrix(rnorm(n*p),ncol=p)
# beta=rnorm(p)/sqrt(p)
# beta_multiresponse=matrix(rnorm(p*k),ncol=k)/(sqrt(p*k))
# 
# beta[11:30]=0
# beta_multiresponse[11:30,]=0
# 
# y_gaussian=rnorm(X%*%beta)+rnorm(n)
# y_count=rpois(n,lambda=exp(X%*%beta))
# y_binary=rbinom(prob=1/(1+exp(-X%*%beta)),size=1,n=n)
# y_multinomial=rep(0,n)
# for (i in 1:n){
#   y_multinomial[i]=which(rmultinom(exp(X[i,,drop=FALSE]%*%beta_multiresponse)/sum(exp(X[i,,drop=FALSE]%*%beta_multiresponse)),n=1,size=1)!=0)
# }
# y_multiresponse=X%*%beta_multiresponse+rnorm(n*k)
# 
# grouping=list(1:1,c(2:10),11:30,31:50)
# penaltyFactor=c(0,1,1,1)
# save(n,p,k,X,beta,beta_multiresponse,y_gaussian,y_count,y_binary,y_multinomial,y_multiresponse,penaltyFactor,grouping,file = "data/exampleData.RData")

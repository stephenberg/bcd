n=100
p=20
X=matrix(rnorm(n*p),ncol=20)

beta=rnorm(20)

y=X%*%beta+rnorm(100)


library(bcd)


# bcdFit<-function(X,
#                  y,
#                  k=1,
#                  family = c("gaussian", "multinomial", "poisson", "logistic"),
#                  groups,
#                  penaltyFactor,
#                  lambda=NULL,
#                  lambdaMinRatio=0.0001,
#                  nLambda=100,
#                  sampleWeights=NULL,
#                  maxit=100000,
#                  tolerance=10^-12,
#                  offset=NULL,
#                  eigenValueTolerance=10^-9,
#                  rescale=TRUE,
#                  useDevTol=TRUE,
#                  devTol=0.9
# )

fit=bcdFit(X=X,
           y=y,
           k=1,
           family="gaussian",
           groups=as.list(1:p),
           penaltyFactor=c(0,rep(1,19)),
           lambda=NULL,
           lambdaMinRatio = 0.0001,
           nLambda = 100,
           sampleWeights = NULL,
           maxit=1000,
           tolerance=10^-8,
           offset=NULL,
           eigenValueTolerance = 10^-9,
           rescale=TRUE,
           devTol = 0.9)
              

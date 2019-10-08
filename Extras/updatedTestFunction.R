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
#                  eigenValueTolerance=10^-8,
#                  rescale=TRUE,
#                  useDevTol=TRUE,
#                  devTol=0.9
# )

fit=fit_bcd(X=X,
           y=y,
           family="gaussian",
           groups=as.list(1:p),
           penaltyFactor=c(0,rep(1,19)),
           lambda=NULL,
           lambdaMinRatio = 0.0001,
           nLambda = 100,
           sampleWeights = NULL,
           maxit=1000,
           eigenValueTolerance = 10^-12,
           tolerance=10^-8,
           offset=NULL,
           rescale=TRUE,
           devTol = 0.9)
              
library(grpregOverlap)
data(pathway.dat)
X <- pathway.dat$expression
group <- pathway.dat$pathways
y <- pathway.dat$mutation
fit_grpreg <- grpregOverlap(X, y, group, penalty = 'grLasso', family = 'binomial')

groupNums=list()
for (i in 1:length(group)){
  groupNums[[i]]=which(fit_grpreg$incidence.mat[i,]!=0)
}
groupNums[[length(groupNums)+1]]=as.integer(1)
X_mod=X[,unique(unlist(groupNums))]
X_mod=cbind(rep(1,50),X_mod)
penaltyFactor=c(rep(1,length(groupNums)-1),0)
bcdfit<- fit_bcd(X = X_mod,
                 y=as.factor(y),
                 family="logistic",
                 groups=groupNums,
                 penaltyFactor = penaltyFactor,
                 eigenValueTolerance=10^-9)
plot(fit_grpreg$loss)
points(bcdfit$deviance,add=TRUE,col="blue")


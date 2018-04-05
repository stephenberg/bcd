ez1z2<-function(alpha,gamma,delta,beta){
  return(beta)
}

ez1ez2<-function(alpha,gamma,delta,beta){
  return((delta+beta)*(gamma+beta))
}

eez1z2ez2z1<-function(alpha,gamma,delta,beta){
  val=0
  val=val+alpha*(delta/(alpha+delta))*(gamma/(alpha+gamma))
  val=val+gamma*(beta/(gamma+beta))*(gamma/(alpha+gamma))
  val=val+delta*(delta/(delta+alpha))*(beta/(delta+beta))
  val=val+beta*(beta/(gamma+beta))*(beta/(beta+delta))
  return(val)
}

all<-function(alpha,gamma,delta,beta){
  return(c(ez1ez2(alpha,gamma,delta,beta),eez1z2ez2z1(alpha,gamma,delta,beta),ez1z2(alpha,gamma,delta,beta)))
}

n=1000000
values=matrix(0,n,3)
for (i in 1:n){
  probs=runif(4)
  probs=probs/sum(probs)
  alpha=probs[1]
  beta=probs[2]
  gamma=probs[3]
  delta=probs[4]
  values[i,]=all(alpha,beta,gamma,delta)
}


any((values[,2]==apply(values,1,min))|(values[,2]==apply(values,1,max)))


prob<-function(z1,z2,eta){
  normalizer=2*exp(eta)+2*exp(-eta)
  return(exp(eta*z1*z2)/normalizer)
}

z1z2<-function(z1,z2,eta){
  return(z1*z2)
}
z1fun<-function(z1,z2,eta){
  return(z1)
}

condExpProd<-function(z1,z2,eta){
  p1=exp(eta)/(exp(-eta)+exp(eta))
  p0=1-p1
  e1=z2*p1+-z2*p0
  e2=z1*p1+-z1*p0
  return(e1*e2)
}

condExp<-function(z1,z2,eta){
  p1=exp(eta)/(exp(-eta)+exp(eta))
  p0=1-p1
  e1=z2*p1+-z2*p0
  return(e1)
}

expectation<-function(eta,fun){
  val=0
  for (z1 in c(-1,1)){
    for (z2 in c(-1,1)){
      val=val+prob(z1,z2,eta)*fun(z1,z2,eta)
    }
  }
  return(val)
}
expectation(0.3,z1fun)
expectation(0.3,condExp)
expectation(0.3,z1z2)
expectation(0.3,condExpProd)

eta=seq(-1,1,0.01)
condExpProdExp=expectation(eta,condExpProd)
prodExp=expectation(eta,z1z2)
plot(eta,prodExp)
points(eta,condExpProdExp,pch=19)



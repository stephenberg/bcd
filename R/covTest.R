expit<-function(alpha){
  return(1/(1+exp(-alpha)))
}
logit<-function(p){
  return(log(p/(1-p)))
}

numerator<-function(h,gamma,z){
  return(exp(sum(h*z)+gamma*z[1]*z[2]+gamma*z[1]*z[3]+gamma*z[2]*z[4]+gamma*z[3]*z[4]))
}


y_z<-function(y,z,alpha,beta){
  pyzMat=cbind(c(1-alpha,alpha),c(1-beta,beta))
  prob=1
  col=0
  row=0
  for (i in 1:4){
    if (z[i]==-1){
      col=1
    }
    else{
      col=2
    }
    if (y[i]==-1){
      row=1
    }
    else{
      row=2
    }
    prob=prob*pyzMat[row,col]
  }
  return(prob)
}
normalize<-function(h,gamma){
  constant=0;
  for (z1 in c(-1,1)){
    for (z2 in c(-1,1)){
      for (z3 in c(-1,1)){
        for (z4 in c(-1,1)){
          z=c(z1,z2,z3,z4)
          constant=constant+numerator(h,gamma,z)
        }
      }
    }
  }
  return(constant)
}

prob<-function(h,gamma,z){
  return(numerator(h,gamma,z)/normalize(h,gamma))
}



covariance<-function(gamma,alpha,beta,gamma1){
  m0=c(0,0)
  m1=c(0,0)
  m2=c(0,0)
  m1m2=matrix(0,2,2)
  m1m1=matrix(0,2,2)
  m2m2=matrix(0,2,2)
  m0m0=m1m1
  m0m2=m1m1
  h=c(0,0)
  gamma=gamma
  print(gamma)
  print(gamma1)
  gamma_grad=0
  val=c(0,0)
  val2=matrix(0,2,2)
  for (z1 in c(-1,1)){
    for (z2 in c(-1,1)){
      for (z3 in c(-1,1)){
        for (z4 in c(-1,1)){
          
          z=c(z1,z2,z3,z4)
          p1234real=prob(h,gamma,z)   
          
          for (y1 in c(-1,1)){
            for (y2 in c(-1,1)){
              for (y3 in c(-1,1)){
                for (y4 in c(-1,1)){
                  y=c(y1,y2,y3,y4)    
                  py_z=y_z(y,z,alpha,beta)
                  p12341234=p1234real*py_z
                  val=val+p12341234*e_y(pseudoGrad,gamma1,y,alpha,beta)
                  t_etz=t_stat(z,gamma)-et_z(c(0,0,0,0),gamma,z)
                  etstat_y=e_y(t_stat,gamma,y,alpha,beta)
                  val2=val2+p12341234*t_etz%*%t(etstat_y)
                }
              }       
            }
          }
        }
      }
    }
  }
  return(list(val=val,val2=val2))
}

condH<-function(y,alpha,beta){
  h=(1-y)/2*log((1-beta)/(1-alpha))+(1+y)/2*log(beta/alpha)
  h=h/2
  return(h)
}

et_z<-function(h,gamma,z){
  val=c(0,0)
  val[1]=tanh(z[1]*gamma+h[1])+tanh(z[2]*gamma+h[2])+tanh(z[3]*gamma+h[3])+tanh(z[4]*gamma+h[4])
  val[2]=z[1]*tanh((z[1]+z[4])*gamma+h[2])+z[2]*tanh((z[2]+z[3])*gamma+h[1])
  val[2]=val[2]+z[1]*tanh((z[1]+z[4])*gamma+h[3])+z[3]*tanh((z[2]+z[3])*gamma+h[1])
  val[2]=val[2]+z[4]*tanh((z[1]+z[4])*gamma+h[2])+z[2]*tanh((z[2]+z[3])*gamma+h[4])
  val[2]=val[2]+z[4]*tanh((z[1]+z[4])*gamma+h[3])+z[3]*tanh((z[2]+z[3])*gamma+h[4])
  val[2]=val[2]/2
  return(val)
}

e_y<-function(fun,gamma,y,alpha,beta){
  val=c(0,0)
  

  h=condH(y,alpha,beta)
  for (z1 in c(-1,1)){
    for (z2 in c(-1,1)){
      for (z3 in c(-1,1)){
        for (z4 in c(-1,1)){
          z=c(z1,z2,z3,z4)
          prob1234=prob(h=h,gamma=gamma,z)
          val=val+fun(z,gamma)*prob1234
        }
      }
    }
  }      
  return(val)
}

t_stat<-function(z,gamma){
  return(c(z[1]+z[2]+z[3]+z[4],z[1]*z[2]+z[2]*z[4]+z[3]*z[1]+z[3]*z[4]))
}

pseudoGrad<-function(z,gamma){
  return(t_stat(z,gamma)-et_z(c(0,0,0,0),gamma,z))
}

alpha=0.1
beta=0.9
gamma=0.4
v=covariance(gamma,alpha,beta,0.35)
e_gradient<-function(gamma1,alpha,beta,gamma){
  return(covariance(gamma,alpha,beta,gamma1)$val[2])
}
numDeriv::grad(e_gradient,0.4,method="Richardson",method.args=list(),side=NA,alpha,beta,0.4)
-v$val2[2,2]


etaSeq=seq(0,1,by=0.02)
eGrads=1:length(etaSeq)
for (i in 1:length(etaSeq)){
  eGrads[i]=e_gradient(etaSeq[i],alpha,beta,gamma)
}
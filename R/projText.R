set.seed(33)
n=100
proj=matrix(0,n,2)
for (i in 1:n){
  x=rnorm(3)
  p1=matrix(rnorm(6),3,2)
  p2=matrix(rnorm(6),3,2)
  p3=matrix(rnorm(6),3,2)
  p1=qr.Q(qr(p1))
  p1=p1%*%t(p1)
  p2=qr.Q(qr(p2))
  p2=p2%*%t(p2)
  p3=qr.Q(qr(p3))
  p3=p3%*%t(p3)
  proj[i,1]=t(x)%*%p1%*%x
  p21=t(x)%*%p2%*%p1%*%x
  p21=p21+t(x)%*%p3%*%p1%*%x
  proj[i,2]=p21/2
  
}
proj

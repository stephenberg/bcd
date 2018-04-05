# X=cbind(c(rep(1,25),rep(0,75)),c(rep(0,25),rep(1,10),rep(0,65)),c(rep(0,35),rep(1,65)))
# b=rnorm(100)
# q1=qr.Q(qr(X))
# q2=q1[,c(1,3)]
# q3=q1[,c(2,3)]
# norm<-function(x){sum(x^2)^0.5}
# norm(q1%*%t(q1)%*%b)
# norm(q2%*%t(q2)%*%b)
# norm(q3%*%t(q3)%*%b)
# 
# J=matrix(1,100,1)
# Iproj=diag(100)-J%*%t(J)/100
# 
# b=b-mean(b)
# norm(q2%*%t(q2)%*%b)
# norm(q3%*%t(q3)%*%b)
# q4=(qr.Q(qr(Iproj%*%X)))
# norm(q4%*%t(q4)%*%b)
# #norm(q4[,2:3]%*%t(q4[,2:3])%*%b)
# #norm(q4[,1:2]%*%t(q4[,1:2])%*%b)
# #norm(q4[,c(1,3)]%*%t(q4[,c(1,3)])%*%b)
# 
# Xt=X
# Xt[,1]=rep(1,100)
# q4=(qr.Q(qr(Iproj%*%X[,2:3])))
# norm(q4%*%t(q4)%*%b)
# norm(Xt%*%solve(t(Xt)%*%Xt)%*%t(Xt)%*%b)
# norm(X%*%solve(t(X)%*%X)%*%t(X)%*%b)
# norm(Iproj%*%Xt%*%solve(t(Xt)%*%Xt)%*%t(Xt)%*%b)
# norm(Iproj%*%X%*%solve(t(X)%*%X)%*%t(X)%*%b)
# norm(t(q4)%*%b)
# 
# #norm(q4[,2:3]%*%t(q4[,2:3])%*%b)
# #norm(q4[,1:2]%*%t(q4[,1:2])%*%b)
# #norm(q4[,c(1,3)]%*%t(q4[,c(1,3)])%*%b)
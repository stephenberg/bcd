# # probs=1:10
# # 
# # maxEig=1:n;
# # maxEig[1]=0
# # hess=matrix(0,10,10)
# # for (i in 2:100000){
# #   probs=runif(10)
# #   probs=probs/sum(probs)
# #   hess=diag(probs)-probs%*%t(probs)
# #   maxEig[i]=max(maxEig[i-1],max(eigen(hess)$values))
# # }
# 
# n=2000
# eigList=1:80
# for (i in 2:80){
#   probs=runif(i)
#   maxEig=0
#   for (j in 1:n){
#       probs=runif(i)
#       probs=probs/sum(probs)
#       hess=diag(probs)-probs%*%t(probs)
#       maxEig=max(maxEig,max(eigen(hess)$values))
#   }
#   eigList[i]=maxEig
# }
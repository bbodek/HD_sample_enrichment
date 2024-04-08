library(MASS)


sigma<-matrix(c(1,0.5,0.5,1),nrow=2,byrow = T)
R<-chol(sigma)


# mvn.sim is a function taking inputs:
#   n - the number of simulated values to return
#   R - the cholesky decomposition of a covariance matrix for multivariate normal dist
#   mu - mean vector for multivariate normal dist
# The function returns a matrix of n simulated values from a multivariate normal 
# distribution with mean vector mu and covariance R%*%t(R)
mvn.sim<-function(n,R,mu){
  X = t(R)%*%matrix(rnorm(n*2),2)+matrix(rep(mu,n),nrow=2,byrow=F)
  X = t(X)
  }

lme.hyp.test<-function(data,f,)




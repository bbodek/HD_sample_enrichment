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

# extract.lme.parameters is a function taking inputs:
#   model - a linear mixed effectds model with 1 random slope, 1 random intercept,
#           1 fixed slope, and 1 fixed intercept
#   name - the desired name of the model to use to reference estimates
# The function returns a dataframe containing the estimates of the
# fixed effect parameters, variance of the random effects, residual variance
# and covariance of the random effects
extract.lme.parameters<-function(model,name){
  coef<-fixef(model)
  fixed.int<-coef[[1]]
  fixed.slope<-coef[[2]]
  vc<-as.data.frame(VarCorr(model))
  var.rndm.int<-vc[1,4]
  var.rndm.slope<-vc[2,4]
  covar.rndm<-vc[3,4]
  var.rndm.int<-vc[4,4]
  data.frame(model = name,fixed.int = fixed.int, fixed.slope = fixed.slope, 
             var.rndm.int = var.rndm.int,var.rndm.slope = var.rndm.slope,
             covar.rndm = covar.rndm)
  }






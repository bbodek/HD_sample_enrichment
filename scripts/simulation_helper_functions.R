library(tidyverse)

# mvn.sim is a function taking inputs:
#   n - the number of simulated values to return
#   R - the cholesky decomposition of a covariance matrix for multivariate normal dist
#   mu - mean vector for multivariate normal dist
# The function returns a matrix of n simulated values from a multivariate normal 
# distribution with mean vector mu and covariance R%*%t(R)
mvn.sim<-function(n,R,mu){
  X = t(R)%*%matrix(rnorm(n*5),nrow = 5)+matrix(rep(mu,n),nrow=5,byrow=F)
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
  var.res<-vc[4,4]
  data.frame(model = name,fixed.int = fixed.int, fixed.slope = fixed.slope, 
             var.rndm.int = var.rndm.int,var.rndm.slope = var.rndm.slope,
             covar.rndm = covar.rndm,var.res = var.res)
}

# lme.y.dist is a function taking inputs:
#   Z - a matrix of covariates for a single individual
#   (including a column of 1's for the intercept)
#   beta - a list containing the model coefficients from an LME
#   D - the covariance matrix of the random effects
#   sigma.2 - the estimated residual variance
# The function returns the resulting expected Y vector mu and the associate 
#  covariance matrix cov.Y
lme.y.dist<-function(Z,beta, D, sigma.2){
  # calculate mean vector 
  mu<-Z%*%beta
  # create identity matrix of size n = # of rows in Z
  I<-diag(nrow(Z))
  cov.Y<-Z%*%D%*%t(Z)+sigma.2*I
  return(list(mu = mu, cov.Y = cov.Y))
}

# sim.param.inputs is a function taking inputs:
#   parameter.df - a dataframe containing the parameter estimates from LME models
#   model - the name of the model to simulate from
#   treatment.flag - indicates if want to simulate from the control (=0)
#      or treatment group (=1)
# The function returns the resulting mu vector and R matrix 
# (the cholsky decomposition of the covariance matrix which is used as input into 
# the function to simulate from a multivariate normal distribution)
sim.param.inputs<-function(parameter.df,subcohort,treatment.flag = 0){
  # extract fixed effects
  beta<-unlist(parameter.df%>%filter(model == subcohort)%>%
                 dplyr::select(fixed.int,fixed.slope))
  # simulate treatment effect of 30% if treatment.flag == 1 by multiplying 
  # slope by 0.7
  if(treatment.flag==1){
    beta[2] = 0.7*beta[2]
  }
  # create matrix of covariates
  Z<-matrix(c(1,0,1,0.5,1,1,1,1.5,1,2),nrow=5,byrow=T)
  # extract variance and covariance of random effects
  sigma11<-unlist(parameter.df%>%filter(model == subcohort)%>%
                    dplyr::select(var.rndm.int))
  sigma22<-unlist(parameter.df%>%filter(model == subcohort)%>%
                    dplyr::select(var.rndm.slope))
  sigma12<-unlist(parameter.df%>%filter(model == subcohort)%>%
                    dplyr::select(covar.rndm))
  sigma<-unlist(parameter.df%>%filter(model == subcohort)%>%dplyr::select(var.res))
  D<-matrix(c(sigma11,sigma12,sigma12,sigma22),nrow=2,byrow=T)
  # extract residual variance
  sigma.2<-unlist(parameter.df%>%filter(model == subcohort)%>%
                    dplyr::select(var.res))
  # calculate mu and cov.Y
  params<-lme.y.dist(Z,beta,D,sigma.2)
  mu<-params$mu
  R<-chol(params$cov.Y)
  return(list(mu = mu, R = R))
}


# hyp.test is a function taking inputs:
#   parameter.df - a dataframe containing the parameter estimates from LME models
#   model - the name of the model to simulate from
# The function returns the resulting mu vector and R matrix 
# (the cholsky decomposition of the covariance matrix which is used as input into 
# the function to simulate from a multivariate normal distribution)



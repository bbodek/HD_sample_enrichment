library(tidyverse)
library(lmerTest)

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

# sim.hypothesis.test is a function taking inputs:
#   n - total number of participants
#   n.treat - number of participants assigned to the treatment group
#   R - the cholesky decomposition of the covariance matrix
#   mu.treat - the mean vector for the treated group
#   mu.control - the mean vector for the untreated group
# The function fits simulates data for the treatment and control groups,
# fit an LME model with random slope and intercept and conducts a
# hypothesis test on the parameter of the treatment by time interaction turn
# It returns a 1 if the resulting hypothesis test successfully rejects 
# the null and a 0 if it fails to reject the null
sim.hypothesis.test<-function(n.control,n.treat,R,mu.control,mu.treat){
  n<-n.treat+n.control
  # create simulated data for treated group under alternative hypothesis
  treated<-mvn.sim(n.treat,R,mu.treat)
  # create simulated data for control group
  control<-mvn.sim(n.control,R,mu.control)
  # turn into a concatenated vector
  outcome.vec<-c(t(treated),t(control))
  # create vector of subjids
  subjid<-rep(1:n,each=5)
  # create vector of timepoints at each observation
  vis.yr<-rep(c(0,0.5,1,1.5,2),n)
  # vector of treatment flags to indicate group membership
  trt.flg<-c(rep(1,n.treat*5),rep(0,n.control*5))
  # fit model and return p value for coefficient relating to treatment effect
  model<-lmer(outcome.vec~1 + vis.yr + vis.yr:trt.flg + (1+vis.yr|subjid))
  p.val<-summary(model)$coefficients["vis.yr:trt.flg", "Pr(>|t|)"]
  # return 1 if p value < 0.05, 0 otherwise
  return(ifelse(p.val<=0.05,1,0))
}

# approx.sample.calc is a function taking inputs:
#   parameter.df- the dataframe of parameter estimates
#   param.model - the name of the model from which to grab parameter estimates 
#           for sample size calculations
#   effect.size - the desired effect size as a rate of change in slope (between 0 and 1)
#   trt.prop - the proportion of participants recieving treatment
#   power - the desired power (between 0 and 1)
#   alpha - the desired type I error rate 
#   num.visits - the total number of study visits
#   study.duration - study duration in years
# The function returns an approximate sample size necessary to reach the desired
# power while maintaining the desired type I error rate.
# See "Applied Longitudinal Analysis" by Fitzmaurice, Laird, and Ware pg 585
# for more details
approx.sample.calc<-function(parameter.df, param.model,effect.size,trt.prop=0.5,
                             power = 0.9,alpha = 0.05,
                             num.visits,study.duration){
  # grab the residual variance
  sigma2.epsilon<-unlist(parameter.df%>%
                           filter(model == param.model)%>%
                           select(var.res))
  # grab the variance of the random slope
  g.22<-unlist(parameter.df%>%
                 filter(model == param.model)%>%
                 select(var.rndm.slope))
  # calculate the effect size as a change in slope
  delta<-unlist(parameter.df%>%
                  filter(model == param.model)%>%
                  select(fixed.slope))*effect.size
  # calculate the variance of the individual slope estimate
  time.sum<-study.duration^2*num.visits*(num.visits+1)/(12*(num.visits-1))
  sigma.2.beta<-sigma2.epsilon/time.sum+g.22
  # calculate z values for power and type I error
  z.alpha<-qnorm(1-alpha/2)
  z.power<-qnorm(power)
  # calculate approximate sample size
  N<-sigma.2.beta*(z.alpha+z.power)^2/(trt.prop*(1-trt.prop)*delta^2)
  N  
}


simulate<-function(n.sim,model,lower.n,upper.n, increment){
  column_names <- c("n", "Power", "Alpha")
  # Initialize an empty dataframe with predefined column names but no rows
  summary_df <- data.frame(n = 0, power = 0, alpha = 0)
  # gather parameters needed for the simulation
  R<-sim.param.inputs(params.df,model,treatment.flag = 0)$R
  mu.control<-sim.param.inputs(params.df,model,treatment.flag = 0)$mu
  mu.treat<-sim.param.inputs(params.df,model,treatment.flag = 1)$mu
  # repeat for sequence of sample sizes
  for(n in seq(lower.n,upper.n, by = increment)){
    # run n.sim simulations using parallelization
    res <- foreach(i=1:n.sim, .combine = rbind,.packages = c("lmerTest"),
                   .export=c("sim.hypothesis.test","mvn.sim","extract.lme.parameters",
                             "lme.y.dist","sim.param.inputs","approx.sample.calc"))%dorng% {
      # select number of treated and control participants by assigning treatment 
      # with 50% probability
      n.treat<-rbinom(1,n,0.5)
      n.control<-n-n.treat
      # simulate hypothesis test results under alternative hypothesis
      result.alt<-sim.hypothesis.test(n.control = n.control,n.treat = n.treat,R = R,
                                      mu.control = mu.control,mu.treat = mu.treat)
      # simulated hypothesis test results under the null hypothesis
      result.null<-sim.hypothesis.test(n.control = n.control,n.treat = n.treat,R = R,
                                       mu.control = mu.control,mu.treat = mu.control)
      # create data frame of results
      data.frame(sample_size = n, iteration = i, 
                 hypothesis.test.null = result.null,
                 hypothesis.test.alt = result.alt)
    }
    # turn all results into dataframe
    result_df<-as.data.frame(res)
    power<-mean(result_df$hypothesis.test.alt)
    alpha<-mean(result_df$hypothesis.test.null)
    print(n)
    print(power)
    print(alpha)
    summary_df<-rbind(summary_df,c(n,power,alpha))
  }
  summary_df
}








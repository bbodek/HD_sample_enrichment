###########################
##### 05 - Simulation #####
###########################

## This script runs simulations to estimate the required sample size to reach the
## desired power for each enrichment strategy

library(tidyverse)
library(doParallel)
library(doRNG)

## set up cluster for parallelization of the simulation
# This detects the total number of cores available on the machine
nworkers <- detectCores() 
# use all available cores
cl <- makePSOCKcluster(nworkers)
registerDoParallel(cl)
source("scripts/simulation_helper_functions.R")

# set random seed for reproducibility
set.seed(11)

# read in dataframe of parameters (created in "03 - Parameter Estimation.R")
# these parameters will be used to generate the simulated data
params.df<-read.csv("./data/parameters_estimates.csv")
# set number of simulated datasets generated for each power estimate
n.sim<-10

# calculate the approximate sample sizes needed to reach 90% power for each
# recruitment strategy using an approximate sample size formula
# these will be used as initial values for the simulation
pin.sapprox<-unlist(approx.sample.calc(params.df,"PIN",effect.size = 0.3,
                   num.visits = 5,study.duration = 2,power=0.9))
unenriched.sapprox<-unlist(approx.sample.calc(params.df,"Unenriched",effect.size = 0.3,
                   num.visits = 5,study.duration = 2,power=0.9))
cap.sapprox<-unlist(approx.sample.calc(params.df,"CAP",effect.size = 0.3,
                   num.visits = 5,study.duration = 2,power=0.9))
print(paste0("Approximate sample size for Unenriched cohort is ",
             round(unenriched.sapprox)))
print(paste0("Approximate sample size for PIN cohort is ",
             round(pin.sapprox)))
print(paste0("Approximate sample size for CAP cohort is ",
             round(cap.sapprox)))


### Sample size calculation for the unenriched cohort ###
print('UNENRICHED')
print('***********************************')

unenriched_results<-simulate(n.sim=3600,model="Unenriched",
                             lower.n = round(unenriched.sapprox)-50,
                             upper.n = round(unenriched.sapprox)+50, 
                             increment=10)
write.csv(file="data/unenriched_results.csv",unenriched_results,row.names=F)

### Sample size calculation for the CAP cohort ###
print('CAP')
print('***********************************')
tic()
cap_results<-simulate(n.sim=3600,model="CAP",
                             lower.n = round(cap.sapprox)-50,
                             upper.n = round(cap.sapprox)+50, 
                             increment=10)
write.csv(file="data/cap_results.csv",cap_results,row.names=F)
toc()
### Sample size calculation for the PIN cohort ###
print('PIN')
print('***********************************')
tic()
pin_results<-simulate(n.sim=3600,model="PIN",
                      lower.n = round(pin.sapprox)-50,
                      upper.n = round(pin.sapprox)+50, 
                      increment=10)
write.csv(file = "data/pin_results.csv",pin_results,row.names=F)
toc()
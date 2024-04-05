#####################################
##### 03 - Parameter Estimation #####
#####################################

## This script fits linear mixed-effects models for each sub-cohort to obtain 
## parameter estimates to be used in the power simulation

library(lme4)
# define linear mixed effects model formula to be used
f<- "cuhdrs ~ visdy + (1 | subjid)+ (visdy |subjid)"

lmer(formula = f,data=unenriched.subco)

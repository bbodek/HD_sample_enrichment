#####################################
##### 03 - Parameter Estimation #####
#####################################

## This script fits linear mixed-effects models for each sub-cohort to obtain 
## parameter estimates to be used in the power simulation

library(lme4)
# define linear mixed effects model formula to be used
f<- " ~ experienced.flag.long +num.nests+ (1 | DNA_ID)+ (1 |Beach), data = df.dispersion"

head(baseline_df)

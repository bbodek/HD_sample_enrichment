#####################################
##### 03 - Parameter Estimation #####
#####################################

## This script fits linear mixed-effects models for each sub-cohort to obtain 
## parameter estimates to be used in the power simulation

library(lme4)
library(lmerTest)
library(tidyverse)
source("scripts/simulation_helper_functions.R")

# define linear mixed effects model formula to be used
f<- "cuhdrs ~ 1+ vis.yr + (1 + vis.yr| subjid)"

# select only subjects who have at least 2 study visits and print number of participants
# for each sub cohort
capenriched_2viz <- capenriched_df %>%
  group_by(subjid) %>%
  filter(max(visit) > 1 )%>%ungroup()
nrow(capenriched_2viz%>%dplyr::select(subjid)%>%unique())

pinenriched_2viz <- pinenriched_df %>%
  group_by(subjid) %>%
  filter(max(visit) > 1 )
nrow(pinenriched_2viz%>%dplyr::select(subjid)%>%unique())

unenriched_2viz <- unenriched_df %>%
  group_by(subjid) %>%
  filter(max(visit) > 1)
nrow(unenriched_2viz%>%dplyr::select(subjid)%>%unique())

# create variable vis.yr as visdy/365 to go from days to years
capenriched_2viz$vis.yr<-capenriched_2viz$visdy/365
pinenriched_2viz$vis.yr<-pinenriched_2viz$visdy/365
unenriched_2viz$vis.yr<-unenriched_2viz$visdy/365

# fit linear mixed-effects model for each subcohort
cap.model<-lmer(formula=f,data = capenriched_2viz)
pin.model<-lmer(formula=f,data = pinenriched_2viz)
unenriched.model<-lmer(formula=f,data = unenriched_2viz)

### Assess goodness of fit for each model
library(HLMdiag)
## cap
summary(cap.model)
plot(cap.model)
infl = hlm_influence(cap.model, level = 'subjid')
plot(infl$cooksd)
df.sens<-capenriched_2viz%>%filter(!subjid %in% unlist(infl[infl$cooksd>0.004,'subjid'] ))
cap.model.sens<-lmer(formula=f,data=df.sens)
summary(cap.model.sens)
summary(cap.model)
## cap
summary(pin.model)
plot(pin.model)
infl = hlm_influence(pin.model, level = 'subjid')
plot(infl$cooksd)
df.sens<-pinenriched_2viz%>%filter(!subjid %in% unlist(infl[infl$cooksd>0.00s,'subjid'] ))
pin.model.sens<-lmer(formula=f,data=df.sens)
summary(pin.model.sens)
summary(pin.model)
## unenriched
summary(unenriched.model)
plot(unenriched.model)
infl = hlm_influence(unenriched.model, level = 'subjid')
plot(infl$cooksd)
df.sens<-unenriched_2viz%>%filter(!subjid %in% unlist(infl[infl$cooksd>0.004,'subjid'] ))
unenriched.model.sens<-lmer(formula=f,data=df.sens)
summary(unenriched.model.sens)
summary(unenriched.model)

# extract parameter estimates from each model
# and create unified dataframe for all parameter estimates
cap.params<-extract.lme.parameters(cap.model,"CAP")
pin.params<-extract.lme.parameters(pin.model,"PIN")
unenriched.params<-extract.lme.parameters(unenriched.model,"Unenriched")
params.df<-rbind(cap.params,pin.params,unenriched.params)
# write parameter dataframe to a csv
write.csv(params.df,"./data/parameters_estimates.csv",row.names = F)


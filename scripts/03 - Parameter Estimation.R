#####################################
##### 03 - Parameter Estimation #####
#####################################

## This script fits linear mixed-effects models for each sub-cohort to obtain 
## parameter estimates to be used in the power simulation

library(lme4)
library(lmerTest)
# define linear mixed effects model formula to be used
f<- "cuhdrs ~ 1+ visdy + (1 + visdy| subjid)"


summary(capenriched_df%>%filter(visit==2 & baseline_TFCS>11 & baseline_TMS<20)%>%mutate(delta = baseline_cuhdrs-cuhdrs)%>%select(delta))
unenriched_df%>%filter(baseline_cuhdrs<5)
summary(unenriched_df$baseline_TFCS)


# select only subjects who have at least 2 study visits
capenriched_2_df <- capenriched_df %>%
  group_by(subjid) %>%
  filter(max(visit) > 1 & baseline_TFCS>11 & baseline_TMS<20 & !is.na(baseline_cuhdrs))


test<-selected_df%>%
  filter(baseline_PIN>0.4)s

#capenriched_2_df$cuhdrs<-capenriched_2_df$cuhdrs/(max(capenriched_2_df$baselins_cuhdrs)-min(capenriched_2_df$baseline_cuhdrs))
capenriched_2_df$visdy<-capenriched_2_df$visdy/365

model<-lmer(formula = f,data=capenriched_2_df, REML = T,control = lmerControl(optimizer ="bobyqa",optCtrl=list(maxfun=2e5)))
summary(model)
library(ggplot2)
ggplot(unenriched_df, aes(x = visdy, y = cuhdrs, group = subjid, color = as.factor(subjid))) +
  geom_line() +
  labs(x = "visdy", y = "cuhdrs") +
  theme_minimal()

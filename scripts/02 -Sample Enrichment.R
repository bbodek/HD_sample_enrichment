##################################
##### 02 - Sample Enrichment #####
##################################

## This script calculates the thresholds to be used as sample enrichment critera
## based off of median CAP and PIN values in the eligible participants
## within the ENROLL-HD dataset.
## Then, 3 sub-cohorts are created, one for each enrichment strategy
## (PIN-enrichment, CAP-enrichment, and Unenriched)


library(tidyverse)
print("Calculating Thresholds for PIN and CAP")

# create dataset of only baseline observations
baseline_df = selected_df%>%filter(visit==1)
nrow(baseline_df)

#  print quantiles of CAP at baseline
q.cap<-quantile(baseline_df$baseline_CAP, probs = c(0.25, 0.5, 0.75))
# select threshold as median value
cap.threshold<-as.numeric(q.cap[2])
print("CAP Threshold: ")
print(cap.threshold)

#  print quantiles of PIN at baseline
q.pin<-quantile(baseline_df$baseline_PIN, probs = c(0.25, 0.5, 0.75))
# select threshold as the median value
pin.threshold<-as.numeric(q.pin[2])
print("PIN Threshold: ")
print(pin.threshold)


hist(baseline_df$TMS)
hist(baseline_df$SDMT)
hist(baseline_df$SWR)
hist(baseline_df$TFCS)

# create unenriched sub-cohort
unenriched_df<-selected_df
# create CAP-enriched sub-cohort
capenriched_df<-selected_df%>%
  filter(baseline_CAP>cap.threshold)
# create PIN-enriched sub-cohort
pinenriched_df<-selected_df%>%
  filter(baseline_PIN>pin.threshold)

# total number of eligible participants
print("Total number of eligible participants:")
print(unenriched_df%>%dplyr::select(subjid)%>%unique()%>%nrow())


# number of eligible participants with baseline cap above threshold
print("Number of CAP eligible participants:")
print(capenriched_df%>%dplyr::select(subjid)%>%unique()%>%nrow())
# eligible participants with baseline PIN abovse threshold
print("Number of PIN eligible participants")
print(pinenriched_df%>%dplyr::select(subjid)%>%unique()%>%nrow())



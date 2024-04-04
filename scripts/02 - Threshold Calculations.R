#######################################
##### 02 - Threshold Calculations #####
#######################################

## This script calculates the thresholds to be used as sample enrichment critera
## based off of median CAP and PIN values in the eligible participants
## within the ENROLL-HD dataset

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

# total number of eligible participants
print("Total number of eligible participants:")
print(baseline_df%>%nrow())


# number of eligible participants with baseline cap above threshold
print("Number of CAP eligible participants:")
print(baseline_df%>%filter(baseline_CAP>q_cap[2])%>%nrow())
# eligible participants with baseline PIN abovse threshold
print("Number of PIN eligible participants")
print(baseline_df%>%filter(baseline_PIN>q_pin[2])%>%nrow())



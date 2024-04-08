################################
##### 01 - ENROLL dataprep #####
################################

## This script does initial processing of the ENROLL-HD dataset
## and applies initial inclusion/exclusion critera to create the dataset
## used as input for the simulation analysis


# clear workspace
rm(list = ls())
library(tidyverse)


print("Processing ENROLL-HD Dataset")

# # General participant related data, some of 
# # which may be annually updated. This
# # includes: Demographics, HDCC 
# # (HD clinical characteristics), CAG, Mortality
data_profile = read.csv("G:/projects/enroll-hd/DATA-00000464 - Enroll-HD PDS5-R2 Study Dataset/Enroll-HD-Periodic-Dataset-2020-10-R2/profile.csv", sep = "\t") %>%
  dplyr::select(subjid, sex, race, caghigh) %>%
  mutate(sex = ifelse(sex == "f",1,0), 
         # race = ifelse(race > 3 | race == 2,4,race))
         race = ifelse(race ==1,1,0))


# data about ENROLL
data_enroll = read.csv("G:/projects/enroll-hd/DATA-00000464 - Enroll-HD PDS5-R2 Study Dataset/Enroll-HD-Periodic-Dataset-2020-10-R2/enroll.csv", sep = "\t") %>% 
  dplyr::select(subjid, seq, motscore, sdmt1,tfcscore, age, diagconf,visdy,swrt1)

# combine the datasets 
data_enroll = data_enroll %>%
  left_join(data_profile, by=c("subjid")) %>%
  # rename variables to match PREDICT HD
  rename(visit = seq, TMS = motscore, TFCS = tfcscore,SDMT = sdmt1, cag = caghigh, DCL = diagconf,SWR = swrt1) %>%
  mutate(age = as.numeric(age))%>%
  mutate(cag = as.numeric(cag))

# delete unused datasets
rm(data_profile)
data_enroll$subjid  %>% unique() %>% length() # 21116 unique participants
data_enroll$subjid%>% length() # 67000 observations

################################################################################
## begin applying inclusion and exclusion criteria

# sort data by appointment time
data_enroll = data_enroll[with(data_enroll, order(subjid, age)), ]

# create data with baseline DCL, TMS, and SDMT
data_baseline = data_enroll %>%
  # get first observation for each subject
  group_by(subjid) %>%
  slice_head(n = 1) %>%
  # rename X to baseline_X
  rename(baseline_age = age,
         baseline_DCL = DCL,
         baseline_TMS = TMS,
         baseline_SDMT = SDMT,
         baseline_CAG = cag,
         baseline_TFCS = TFCS,
         baseline_SWR = SWR) %>%
  dplyr::select(subjid, contains("baseline"))

# join baseline data to model dataset by SUBJID
data_enroll = data_baseline %>%
  full_join(data_enroll, ., by = "subjid")
rm(data_baseline)

# create age at diagnosis (diag_age)
data_diag = data_enroll %>%
  # motor onset is said to occur when DCL = 4 for the first time,
  # so calculate age at diagnosis for all subjects who experience onset
  filter(DCL == 4) %>%
  group_by(subjid) %>%
  filter(visit == min(visit)) %>%
  summarise(X = as.numeric(age))

# join diagnosis age data to model dataset by SUBJID
data_enroll = data_diag %>%
  left_join(data_enroll, ., by = "subjid")
rm(data_diag)

# create event time variable and delta
data_enroll = data_enroll %>%
  # delta is 1 if diagnosis is observed, 0 otherwise
  mutate(delta = ifelse(test = is.na(X), yes = 0, no = 1)) %>%
  # event time (T) is time from baseline to diagnosis
  mutate(age = as.numeric(age)) %>%
  mutate(AX = age-X)

# for censored subjects, get age at last visit
data_censored = data_enroll %>%
  # filter: subjects who do not experience motor onset
  filter(delta == 0) %>%
  # calculate age at last visit
  group_by(subjid) %>%
  slice_tail() %>%
  mutate(C = as.numeric(age)) %>%
  dplyr::select(subjid, C) 

# join censoring time to model dataset by SUBJID
data_enroll = data_censored %>%
  left_join(data_enroll, ., by = "subjid")
rm(data_censored)

# create censoring time variable
data_enroll = data_enroll %>%
  # censoring time (C) is time from baseline to dropout
  mutate(AC = age - C,
         W = ifelse(delta == 1, X, C),
         AW = ifelse(delta == 1, AX, AC))

# create CAP score variables at baseline and for all visits
data_enroll = data_enroll %>%
  mutate(CAP = age*(cag-33.66),baseline_CAP = baseline_age * (cag-33.66))

# create PIN score variables at baseline and for all visits
data_enroll = data_enroll %>%
  mutate(PIN = (51*TMS-34*SDMT+7*CAP-883)/1044,baseline_PIN = (51*baseline_TMS-34*baseline_SDMT+7*baseline_CAP-883)/1044)

# calculate cUHDRS at baseline and for all visits
cuhdrs_calc<-function(TFC,TMS,SDMT,SWR){
  cuhdrs<-((TFC-10.4)/1.9-(TMS-29.7)/14.9+(SDMT-28.3)/11.3+(SWR-66.1)/20.1)+10
  return(cuhdrs)
}

data_enroll = data_enroll %>% rowwise()%>%
  mutate(baseline_cuhdrs = cuhdrs_calc(baseline_TFCS,baseline_TMS,baseline_SDMT,baseline_SWR),cuhdrs = cuhdrs_calc(TFCS,TMS,SDMT,SWR))

## APPLY EXCLUSION CRITERIA FROM LONG ET AL (2017)

# number of subjects before applying exclusion criteria
data_enroll$subjid %>% unique() %>% length()

## filter: age >= 18
data_age = data_enroll %>%
  filter(visit == 1) %>%
  filter(baseline_age != "<18") %>%
  dplyr::select(subjid)

selected_df = data_enroll

selected_df = selected_df %>%
  filter(subjid %in% data_age$subjid)
rm(data_age)
selected_df$subjid %>% unique() %>% length() #21086


# filter: CAG >= 36 at baseline 
good_subject_id = selected_df %>% 
  filter(visit==1) %>%
  filter(cag >= 36) %>%
  dplyr::select(subjid)
selected_df = selected_df %>%
  subset(subjid %in% good_subject_id$subjid) 
rm(good_subject_id)
selected_df$subjid %>% unique() %>% length() #16079


# filter: DCL < 4
good_subject_id = selected_df %>% 
  filter(visit==1) %>%
  # mutate(baseline_DCL = factor(baseline_DCL)) %>%
  # summary()
  filter(baseline_DCL < 4) %>%
  dplyr::select(subjid)
selected_df = selected_df %>%
  subset(subjid %in% good_subject_id$subjid) 
rm(good_subject_id)
selected_df$subjid %>% unique() %>% length() # 5717

# filter: baseline TMS != NA
`%!in%` = Negate(`%in%`)

bad_subject_id = selected_df %>% 
  subset(is.na(TMS) | is.na(SDMT) | is.na(age)) %>%
  dplyr::select(subjid)
selected_df = selected_df %>%
  subset(subjid %!in% bad_subject_id$subjid) 
rm(bad_subject_id)
selected_df$subjid %>% unique() %>% length() #3719

# filter those with an error code for SDMT
bad_subject_id = selected_df %>% 
  subset(SDMT == 9998) %>%
  dplyr::select(subjid)
selected_df = selected_df %>%
  subset(subjid %!in% bad_subject_id$subjid) 
rm(bad_subject_id)
selected_df$subjid %>% unique() %>% length() #3718

# filter those with an error code for SWR
bad_subject_id = selected_df %>% 
  subset(SWR == 9996 | SWR == 9997 | SWR == 9998) %>%
  dplyr::select(subjid)
selected_df = selected_df %>%
  subset(subjid %!in% bad_subject_id$subjid) 
rm(bad_subject_id)
selected_df$subjid %>% unique() %>% length() #3701


# filter those with an error code for SDMT
bad_subject_id = data_enroll %>% 
  subset(SDMT == 9998 | SDMT == 9997 | SDMT == 9996) %>%
  dplyr::select(subjid)
selected_df = selected_df %>%
  subset(subjid %!in% bad_subject_id$subjid) 
rm(bad_subject_id)
selected_df$subjid %>% unique() %>% length() #3701

# filter those without a baseline cuhdrs
bad_subject_id = data_enroll %>% 
  subset(is.na(baseline_cuhdrs)) %>%
  dplyr::select(subjid)
selected_df = selected_df %>%
  subset(subjid %!in% bad_subject_id$subjid) 
rm(bad_subject_id)
selected_df$subjid %>% unique() %>% length() #3701

# filter those with TMS>20 or TFC<11 at baseline
bad_subject_id = data_enroll %>% 
  subset(baseline_TMS>20 | baseline_TFCS<11) %>%
  dplyr::select(subjid)
selected_df = selected_df %>%
  subset(subjid %!in% bad_subject_id$subjid) 
rm(bad_subject_id)
selected_df$subjid %>% unique() %>% length() #3701

# number of subjects after applying exclusion criteria (3701)
selected_df$subjid %>% unique() %>% length() 



## end of inclusion and exclusion criteria
################################################################################

################################################################################
## calculate basic study statistics 
#(number of participants, number of visits, etc)


# calculate proportion observed
selected_df %>%
  # get first observation for each subject
  group_by(subjid) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  dplyr::select(delta) %>%
  summary()

# only 15.41 % had the outcome observed during the study

# total number of visits
nrow(selected_df) #9866
# table of numbesr of visits per patient
table(selected_df%>%group_by(subjid)%>%summarise(max(visit))%>%dplyr::select('max(visit)'))
# avg number of visits
mean(selected_df$visit) # 2.33
# median number of visits
median(selected_df$visit) # 2
# median and mean time between visits
df <- selected_df %>%
  group_by(subjid) %>%
  mutate(Diff = visdy - lag(visdy))%>%filter(visdy!=0)
mean(df$Diff)
median(df$Diff)

summary(selected_df$TMS)

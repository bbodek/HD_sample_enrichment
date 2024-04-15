library(tidyverse)
library(ggplot2)
library(stringr)

#### Reading & Processing Simulation Results ####
#################################################

# create overall list of files
files<-list.files('data/',pattern = "results")
# grab files which are for the alt hypothesis not null hypothesis
files<-files[-grep('null',files)]
files

sim_results<-data.frame()
for(file in files){
  data<-read.csv(paste0("data/",file))
  enrichment<-str_extract(file, "^[^_]+")
  data$enrichment<-enrichment
  sim_results<-rbind(sim_results,data)
}
# calculate monte carlo SE as sqrt(power*(1-power)/n.sim)
sim_results$se<-(sqrt(sim_results$power*(1-sim_results$power)/3600))
sim_results$lower.ci<-sim_results$power-sim_results$se*1.96
sim_results$upper.ci<-sim_results$power+sim_results$se*1.96
sim_results$enrichment<-as.character(sim_results$enrichment)
sim_results$enrichment.type<-ifelse(sim_results$enrichment=='cap','CAP Threshold',
                                     ifelse(sim_results$enrichment=='pin','PIN Threshold','Unenriched'))
# plot resulting power curves
# color blind friendly pallet
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# filter only estimates for sample size smaller than 7521
sim_results<-sim_results%>%filter(n<7520)
# plot power curve
power_plt<-ggplot(sim_results, aes(x = n, y = power,color=enrichment.type)) +
  geom_point(shape=1) +
  labs(title = "Power Curves by Enrichment Strategy", x = "Sample Size", y = "Power") +
  scale_color_manual(values=cbbPalette,
                     name = "Enrichment Critera",
                     labels = c('cap'='CAP Threshold',
                                'pin'='PIN Threshold',
                                'unenriched'='Unenriched'))+
  scale_y_continuous(breaks=seq(0,1.1,0.1))+
  scale_x_continuous(breaks=seq(0,7600,400),minor_breaks = T)+
  # add smoothed curve fit to points
  geom_smooth(method = 'gam',size=0.5,se=F)+
  # add a line at 90% power
  geom_hline(yintercept = 0.9,color='red',linetype='dashed')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.3, hjust=0.2),
        plot.title=element_text(hjust=0.5))

power_plt
ggsave("./figures/powercurve.png",plot=power_plt, width = 10, height = 5.5)

sim_results_2<-sim_results%>%filter(lower.ci>0.85 & upper.ci<0.95)

power_plt_2<-ggplot(sim_results_2, aes(x = n, y = power,color=enrichment.type)) +
  geom_point() +
  geom_errorbar(aes(ymin = power-se,ymax=power+se),width=20)+
  #geom_crossbar(aes(y=power, ymin = lower.ci, ymax = upper.ci), width = 1)+
  labs(title = "Power Estimates with Standard Error by Enrichment Type", x = "Sample Size", y = "Power") +
  scale_color_manual(values=cbbPalette,guide='none')+
  theme(plot.title=element_text(hjust=0.5))+
  ylim(0.86,0.94)+
  facet_wrap(~enrichment.type, scales = 'free_x',ncol = 3)+
  # add a line at 90% power
  geom_hline(yintercept = 0.9,color='red',linetype='dashed')
power_plt_2
ggsave("./figures/power.est.png",plot=power_plt_2,width = 10, height = 5.5)

# calculate approximate sample size range for 90% power
# as the smallest n such that estimated power + se < 90% to the largest n such that 
# estimated power - se > 90%
samplesize.range<-function(enrichment.strat,results_df){
  df<-results_df%>%filter(enrichment == enrichment.strat)
  lower.df<-df%>%filter(power+se>0.9)
  upper.df<-df%>%filter(power-se<0.9)
  lower<-lower.df%>%select(n)%>%min()
  upper<-upper.df%>%select(n)%>%max()
  return(list(lower = lower,upper = upper))
}
# estimated sample size range for cap
samplesize.range("cap",sim_results)
# estimated sample size range for pin
samplesize.range("pin",sim_results)
# estimated sample size range for unenriched
samplesize.range("unenriched",sim_results)


### Analyze Type I error
# create overall list of files
files_null<-list.files('data/',pattern = "results")
# grab files which are for the alt hypothesis not null hypothesis
files_null<-files_null[grep('null',files_null)]
files_null
null_results<-data.frame()
for(file in files_null){
  data<-read.csv(paste0("data/",file))
  enrichment<-str_extract(file, "^[^_]+")
  data$enrichment<-enrichment
  null_results<-rbind(null_results,data)
}
null_results<-null_results%>%filter(power>0)
summary(null_results$power)
sd(null_results$power)
                 
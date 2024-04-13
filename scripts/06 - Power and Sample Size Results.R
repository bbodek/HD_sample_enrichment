library(tidyverse)
library(ggplot2)
library(stringr)

#### Reading & Processing Simulation Results ####
#################################################

# create overall list of files
files<-list.files('data/',pattern = "results")

sim_results<-data.frame()
for(file in files){
  data<-read.csv(paste0("data/",file))
  enrichment<-str_extract(file, "^[^_]+")
  data$enrichement<-enrichment
  sim_results<-rbind(sim_results,data)
}
# calculate monte carlo SE as sqrt(power*(1-power)/n.sim)
sim_results$se<-(sqrt(sim_results$power*(1-sim_results$power)/3600))
sim_results$lower.ci<-sim_results$power-sim_results$se*1.96
sim_results$upper.ci<-sim_results$power+sim_results$se*1.96
sim_results$enrichement.type<-as.character(sim_results$enrichement)
# plot resulting power curves
# color blind friendly pallet
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# plot power curve
power_plt<-ggplot(sim_results, aes(x = n, y = power,color=enrichement.type)) +
  geom_point(shape=1) +
  labs(title = "Power Curves by Enrichment Strategy", x = "Sample Size", y = "Power") +
  scale_color_manual(values=cbbPalette,
                     name = "Enrichment Critera",
                     labels = c('cap'='CAP Threshold',
                                'pin'='PIN Threshold',
                                'unenriched'='Unenriched'))+
  scale_y_continuous(breaks=seq(0,1.1,0.1))+
  scale_x_continuous(breaks=seq(0,3000,200),minor_breaks = NULL)+
  # add smoothed curve fit to points
  geom_smooth(method = 'loess',size=0.5,se=F)+
  # add a line at 90% power
  geom_hline(yintercept = 0.9,color='red',linetype='dashed')+
  theme_minimal()
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

power_plt




                 
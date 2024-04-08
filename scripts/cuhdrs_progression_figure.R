library(ggplot2)

# create dataframe of times at which we want predicted cUHDRS values
predict.df<-data.frame(subjid = '1',vis.yr = c(0,0.5,1,1.5,2,2.5,3,3.5,4))

# generate predicted average cUHDRS values for each subcohort
predict.pin.df<-data.frame(subjid = '1',vis.yr = predict.df$vis.yr,
                           cuhdrs = predict(pin.model,newdata = predict.df,
                                            re.form = NA))
predict.cap.df<-data.frame(subjid = '1',vis.yr = predict.df$vis.yr,
                           cuhdrs = predict(cap.model,newdata = predict.df,
                                            re.form = NA))
predict.unenriched.df<-data.frame(subjid='1',vis.yr = predict.df$vis.yr,
                           cuhdrs = predict(unenriched.model,newdata = predict.df,
                                            re.form = NA))

# function to plot individual progression for a sample of participants along with 
# overall projected population project
plot.cuhdrs.progression<-function(df,sample_size,pred_df){
  sample_subj <- sample(unique(df$subjid), sample_size)
  plt<-ggplot(df[df$subjid %in% sample_subj, ], aes(x = vis.yr, y = cuhdrs, group = subjid)) +
    geom_line(color='grey') +  # Plot individual lines
    geom_point(color='grey',size=1) +  # Plot individual points
    geom_line(aes(y = cuhdrs, x = vis.yr), color = "black", linewidth = 1,pred_df) +  # Plot average line
    labs(y = "cUHDRS",x='Year in Study') +
    scale_x_continuous(name = "Year in Study",breaks = c(0,0.5,1,1.5,2,2.5,3,3.5,4),limits=c(0,4),expand=c(0,0))+
    #xlim(0,4)+
    scale_y_continuous(limits=c(5,20),expand=c(0,0))+
    ggtitle("cUHDRS Progression Over Time In Study") +
    theme_bw()
  plt}

plt<-plot.cuhdrs.progression(pinenriched_2viz,sample_size = 200,predict.pin.df)
ggsave("./figures/cuhdrs_time_pin.png",plot=plt)

plt<-plot.cuhdrs.progression(capenriched_2viz,sample_size = 200,predict.cap.df)
ggsave("./figures/cuhdrs_time_cap.png",plot=plt)

plt<-plot.cuhdrs.progression(unenriched_2viz,sample_size = 200,predict.unenriched.df)
ggsave("./figures/cuhdrs_time_unenriched.png",plot=plt)

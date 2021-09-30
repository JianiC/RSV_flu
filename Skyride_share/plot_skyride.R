# Script to visualize csv output file(export from Tracer )from BEAST Bayesian Skyride coalescent analysis
# Input: csv file (remove first line with	GMRF Skyride: xxxxxxxx.log)
# Output: ggplot of skyride wiht 95% HPD 


library(ggplot2)
library(lubridate)


dat <- read.table('cluster_result/rsvb_12_14/rsvb_12_14/rsvb_12_14_skyride.txt', skip = 1, header =TRUE, sep ='\t')

dat$date <- as.Date(paste(date_decimal(dat$Time)),"%Y-%m-%d")

p1 <- ggplot(dat, aes(y=Mean, x=date)) + geom_line() + 
  geom_ribbon(aes(ymin=Lower, ymax=Upper, alpha = 0.2)) +
  theme(legend.position = "bottom") + labs(alpha="95% hpd") + 
  labs(title = "Skyride coalescent for RSV-B in 2fi12-2014 season") + 
  scale_y_continuous(name ="Log mean pop (N)", trans = 'log10') 


p1

## estimate R0


dat_12 <- dat %>%
  filter(as.numeric(Time) >= 2012 & as.numeric(Time) <=2012.7 )

ggplot(dat_12, aes(y=Mean, x=date)) + geom_line()
model_12 <- lm(Median ~ Time, data = dat_12)
print(model_12)
k = 24.82

## function
plot_skyride<-function(skyride_df,title){
  dat <- read.table(skyride_df, skip = 1, header =TRUE, sep ='\t')
  dat$date <- as.Date(paste(date_decimal(dat$Time)),"%Y-%m-%d")
  p1 <- ggplot(dat, aes(y=Mean, x=date)) + geom_line() + 
    geom_ribbon(aes(ymin=Lower, ymax=Upper, alpha = 0.2)) +
    theme(legend.position = "bottom") + labs(alpha="95% hpd") + 
    labs(title = title) + 
    scale_y_continuous(name ="Log mean pop (N)", trans = 'log10')+
    theme_bw()
  
  return(p1)
}
H1<-plot_skyride("cluster_result/Skyride_0916/gisaid_IAV_H1_12_19/gisaid_IAV_H1_skyride.txt","IAV-H1 GMRF Skyride")
H3<-plot_skyride("cluster_result/Skyride_0916/gisaid_IAV_H3_12_19/gisaid_IAV_H3_skyride.txt","IAV-H3 GMRF Skyride")
rsva<-plot_skyride("cluster_result/Skyride_0916/rsva_12_19season/rsva_12_19_season_skyride.txt","RSV-A GMRF Skyride")
IBV_Y<-plot_skyride("cluster_result/Skyride_0916/gisaid_IBV_Y_12_19/gisaid_IBV_Y_skyride.txt","IBV-Yamagata GMRF Skyride")
IBV_V<-plot_skyride("cluster_result/Skyride_0916/gisaid_IBV_V_12_19/gisaid_IBV_V_skyride.txt","IBV-Victoria GMRF Skyride")
rsvb<-plot_skyride("cluster_result/Skyride_0916/rsvb_12_19season/rsvb_12_19_season_skyride.txt","RSV-B GMRF Skyride")

ggarrange(H1,H3,IBV_V,IBV_Y,ncol=1,common.legend = TRUE)
ggarrange(rsva,rsvb,ncol=1,common.legend = TRUE)
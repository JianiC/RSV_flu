## Skyride visulization
## RSV skyride
ON_skyride<-read.table('Genetic_data/flu_RSV_15-19_north_hemisphere/rsv_10_19_global/rsv_beast/rsv_beast2/gisaid_ON1_north_100901_190430_north_G/gisaid_ON1_north_100901_190430_skyride', skip = 1, header =TRUE, sep ='\t')
BA_skyride<-read.table('Genetic_data/flu_RSV_15-19_north_hemisphere/rsv_10_19_global/rsv_beast/rsv_beast2/gisaid_BA_north_100901_190401_G/gisaid_BA_north_100901_190430_skyride.txt', skip = 1, header =TRUE, sep ='\t')

ON_skyride%>%
  dplyr::rename(ON_mean=Mean,ON_median=Median,ON_upper=Upper,ON_lower=Lower,ON_time=Time)->ON_skyride2

BA_skyride%>%
  dplyr::rename(BA_mean=Mean,BA_median=Median,BA_upper=Upper,BA_lower=Lower,BA_time=Time)->BA_skyride2
col2<-c("ON"="#66A61E","BA"="#1B9E77")
ON_skyride2%>%
  cbind(BA_skyride2)%>%
  ggplot()+
  geom_ribbon(aes(x = ON_time, ymin = ON_lower, ymax = ON_upper), fill = "#66A61E", alpha = 0.4) +
  geom_line(aes(x = ON_time, y = ON_mean, color = "ON")) +
  geom_ribbon(aes(x = BA_time, ymin = BA_lower, ymax = BA_upper),  fill = "#1B9E77", alpha = 0.4) +
  geom_line(aes(x = BA_time, y = BA_mean, color = "BA"))+
  theme_bw()+
  ylab("Effective number of RSV infections")+
  xlab("Time")+
  scale_x_continuous(breaks = seq(2010,2019,1),position = 'bottom',limits =c(2009,2019) )+
  scale_colour_manual(name="",values = col2)+
  theme(legend.position=c(0.1,0.8))->skyride_rsv

## fluA skyride
H3_skyride<-read.table('Genetic_data/flu_RSV_15-19_north_hemisphere/rsv_10_19_global/flu_beast2/gisaid_H3_100901_190430/gisaid_H3_100901_190401_skyride.txt', skip = 1, header =TRUE, sep ='\t')
H1_skyride<-read.table('Genetic_data/flu_RSV_15-19_north_hemisphere/rsv_10_19_global/flu_beast2/gisaid_H1_100901_190430/gisaid_H1_100901_190430_skyride.txt', skip = 1, header =TRUE, sep ='\t')

H3_skyride%>%
  dplyr::rename(H3_mean=Mean,H3_median=Median,H3_upper=Upper,H3_lower=Lower,H3_time=Time)->H3_skyride2

H1_skyride%>%
  dplyr::rename(H1_mean=Mean,H1_median=Median,H1_upper=Upper,H1_lower=Lower,H1_time=Time)->H1_skyride2
col2<-c("H1"="#E6AB02","H3"="#D95F02")
H3_skyride2%>%
  cbind(H1_skyride2)%>%
  ggplot()+
  geom_ribbon(aes(x = H3_time, ymin = H3_lower, ymax = H3_upper), fill = "#D95F02", alpha = 0.25) +
  geom_line(aes(x = H3_time, y = H3_mean,, color = "H3")) +
  geom_ribbon(aes(x = H1_time, ymin = H1_lower, ymax = H1_upper), 
              fill = "#E6AB02", alpha = 0.25) +
  geom_line(aes(x = H1_time, y = H1_mean, color = "H1"))+
  ylim(0,900)+
  theme_bw()+
  scale_x_continuous(breaks = seq(2010,2019,1),position = 'bottom',limits =c(2009,2019) )+
  ylab("Effective number of FluA infections")+
  xlab("Time")+
  scale_colour_manual(name="",values = col2)+
  theme(legend.position=c(0.1,0.8))->skyride_fluA


## fluB skyride
Vic_skyride<-read.table('Genetic_data/flu_RSV_15-19_north_hemisphere/rsv_10_19_global/flu_beast2/gisaid_Vic_100901_190430/gisaid_Vic_100901_190430_skyride.txt', skip = 1, header =TRUE, sep ='\t')
Yam_skyride<-read.table('Genetic_data/flu_RSV_15-19_north_hemisphere/rsv_10_19_global/flu_beast2/gisaid_Yam_100901_190430/gisaid_Yam_100901_190430_skyride.txt', skip = 1, header =TRUE, sep ='\t')

Vic_skyride%>%
  dplyr::rename(Vic_mean=Mean,Vic_median=Median,Vic_upper=Upper,Vic_lower=Lower,Vic_time=Time)->Vic_skyride2

Yam_skyride%>%
  dplyr::rename(Yam_mean=Mean,Yam_median=Median,Yam_upper=Upper,Yam_lower=Lower,Yam_time=Time)->Yam_skyride2

col2<-c("Victoria"="#E7298A","Yamagata"="#7570B3")
Vic_skyride2%>%
  cbind(Yam_skyride2)%>%
  ggplot()+
  geom_ribbon(aes(x = Vic_time, ymin = Vic_lower, ymax = Vic_upper), fill = "#E7298A", alpha = 0.25) +
  geom_line(aes(x = Vic_time, y = Vic_mean,color = "Victoria") ) +
  geom_ribbon(aes(x = Yam_time, ymin = Yam_lower, ymax = Yam_upper), 
              fill = "#7570B3", alpha = 0.25) +
  geom_line(aes(x = Yam_time, y = Yam_mean, color = "Yamagata"))+
  ylim(0,800)+
  theme_bw()+
  scale_x_continuous(breaks = seq(2010,2019,1),position = 'bottom',limits =c(2009,2019))+
  ylab("Effective number of FluB infections")+
  xlab("Time")+
  scale_colour_manual(name="",values = col2)+
  theme(legend.position=c(0.1,0.8))->skyride_fluB

RSV_national_antigen<-RSV_national%>%
  filter(TestType=="1")%>%
  ungroup()%>%
  mutate(type="RSV")

## plot surveliance
flu_national%>%
  ungroup()%>%
  dplyr::select(date,fluApos_national,fluBpos_national)%>%
  melt(id.vars = "date",
       variable.name = "type",
       value.name = "case")%>%
  ggplot(aes(x=date, y=case, fill= type)) +
  geom_area( )+ 
  theme_bw()+
  geom_area(data=RSV_national_antigen,aes(x=date,y=RSVpos_national),alpha=0.6)+
  scale_fill_manual(name = "",values=c("#E6AB02","#E7298A","#66A61E"), labels = c("FluA", "FluB", "RSV"))+
  scale_x_continuous(breaks = seq(2010,2019,1),position = 'bottom',limits =c(2009,2019))+
  theme(legend.position=c(0.1,0.8))+
  ylab("Flu and RSV in US")+
  xlab("Time")->survel_RSV_fluAB

Plot<-cowplot::plot_grid(survel_RSV_fluAB,skyride_rsv,skyride_fluA,skyride_fluB, ncol=1, align='v',rel_heights = c(1,1,1,1))
  
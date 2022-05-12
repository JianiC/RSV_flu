library(ggplot2)
library(ggpubr)
library(ggtree)
library(treeio)
library(scales)
##

setwd("/Users/jianichen/Dropbox/RSV_flu/RSV_flu_git")
#tree file
BA<-read.beast("Genetic_analysis/RSV_flu_Beast/RSV/gisaid_BA_north_100901_190430_anot.tree")
ON<-read.beast("Genetic_analysis/RSV_flu_Beast/RSV/gisaid_ON1_north_100901_190430_anot.tree")
H1<-read.beast("Genetic_analysis/RSV_flu_Beast/Flu/gisaid_H1_100901_190430_anot.tree")
H3<-read.beast("Genetic_analysis/RSV_flu_Beast/Flu/gisaid_H3_100901_190430_anot.tree")
Vic<-read.beast("Genetic_analysis/RSV_flu_Beast/Flu/gisaid_Vic_100901_190430_anot.tree")
Yam<-read.beast("Genetic_analysis/RSV_flu_Beast/Flu/gisaid_Yam_100901_190430_anot.tree")

BA_sta<-read.csv("Genetic_analysis/RSV_flu_Beast/RSV/gisaid_BA_out.skylines",sep="\t")
ON_sta<-read.csv("Genetic_analysis/RSV_flu_Beast/RSV/gisaid_ON1_out.skylines",sep="\t")
H1_sta<-read.csv("Genetic_analysis/RSV_flu_Beast/Flu/gisaid_H1_out.skylines",sep="\t")
H3_sta<-read.csv("Genetic_analysis/RSV_flu_Beast/Flu/gisaid_H3_out.skylines",sep="\t")
Vic_sta<-read.csv("Genetic_analysis/RSV_flu_Beast/Flu/gisaid_Vic_out.skylines",sep="\t")
Yam_sta<-read.csv("Genetic_analysis/RSV_flu_Beast/Flu/gisaid_Yam_out.skylines",sep="\t")


p1<-ggtree(BA, mrsd="2019-4-3",size=0.2) + 
  #geom_range(range='height_0.95_HPD', color='blue', alpha=.3, size=2)+
  #theme_tree2()+
  theme(panel.grid.major.x=element_line(color="grey20", linetype="dotted", size=.3),
        text = element_text(size = 8),
        aspect.ratio = 0.5)+
  scale_x_continuous(breaks = seq(2010,2019,1),position = 'bottom',limits =c(2008,2019))+
  annotate("text",x=0.1,y=-0.8,label="BA",color="red")


  geom_text(aes(label = "BA", x = 0.1, y =0))

  
  
  
  
  
sta_names <- c(
  `tmrca` = "TMRCA",
  `div` = "Diversity")
BA_sta%>%
  filter(statistic=="tmrca")%>%
  ggplot(aes(y=mean, x=time))+
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper, alpha = 0.2))+
  #facet_wrap(~statistic,ncol=1,strip.position="left",labeller = as_labeller(sta_names))+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(breaks = seq(2010,2019,1),position = 'bottom',limits =c(2008,2019))+
  theme(panel.grid.major.x=element_line(color="grey80",size=.3),
        aspect.ratio = 0.2)+
  ylab("TMRCA")->p2

BA_sta%>%
  filter(statistic=="div")%>%
  ggplot(aes(y=mean, x=time))+
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper, alpha = 0.2))+
  #facet_wrap(~statistic,ncol=1,strip.position="left",labeller = as_labeller(sta_names))+
  theme_bw()+
  theme(legend.position = "none")+
  scale_x_continuous(breaks = seq(2010,2019,1),position = 'bottom',limits =c(2008,2019))+
  theme(panel.grid.major.x=element_line(color="grey80", linetype="dotted", size=.3))+
  ylab("Diversity")+
  gg.theme +
  theme(text = element_text(size = 8),
        legend.spacing.y = unit(0.1, "lines"),
        legend.key = element_blank(), 
        legend.background = element_blank(),
        legend.position = "None")->p3

Plot<-cowplot::plot_grid(p1,p2,p3, ncol=1, align='v',rel_heights = c(2,1,1))
############
## plot everything together
## plot everything together
tree_sta<-function(tree,recent_d,sta){
  
  p1<-ggtree(tree, mrsd=recent_d,size=0.2) + 
    #geom_range(range='height_0.95_HPD', color='blue', alpha=.3, size=2)+
    theme(panel.grid.major.x=element_line(color="grey80",  size=.3))+
    scale_x_continuous(breaks = seq(2010,2019,2),position = 'bottom',limits =c(2008,2019))+
    theme(text = element_text(size = 8),
          
          plot.margin = margin(t=1,b = 1,,l=1,r=1) )
  
  sta%>%
    filter(statistic=="tmrca")%>%
    ggplot(aes(y=mean, x=time))+
    geom_line() + 
    geom_ribbon(aes(ymin=lower, ymax=upper, alpha = 0.2))+
    #facet_wrap(~statistic,ncol=1,strip.position="left",labeller = as_labeller(sta_names))+
    theme_bw()+
    theme(legend.position = "none")+
    scale_x_continuous(breaks = seq(2010,2019,2),position = 'bottom',limits =c(2008,2019))+
    theme(panel.grid.major.x=element_line(color="grey80",  size=.3))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "none")+
    ylab("TMRCA")+
    theme(text = element_text(size = 8),
          axis.title=element_text(size=8),
          plot.margin = margin(t=1,b = 1,,l=1,r=1) )+ 
    scale_y_continuous(
      labels = label_number(accuracy = 1),
      limits = c(0,8)
    )->p2
  
  sta%>%
    filter(statistic=="div")%>%
    ggplot(aes(y=mean, x=time))+
    geom_line() + 
    geom_ribbon(aes(ymin=lower, ymax=upper, alpha = 0.2))+
    #facet_wrap(~statistic,ncol=1,strip.position="left",labeller = as_labeller(sta_names))+
    theme_bw()+
    theme(legend.position = "none")+
    scale_x_continuous(breaks = seq(2010,2019,2),position = 'bottom',limits =c(2008,2019))+
    theme(panel.grid.major.x=element_line(color="grey80", size=.3))+
    ylab("Diversity")+
    xlab("Time")+
    scale_y_continuous(
      labels = label_number(accuracy = 1),
      limits = c(0,15)) +
    theme(text = element_text(size = 8),
          axis.title=element_text(size=8),
          plot.margin = margin(t=1,b = 1,,l=1,r=1) )->p3
  
  Plot<-cowplot::plot_grid(p1,p2,p3, ncol=1, align='v',
                           rel_heights = c(1.45,1,1.35))
  return(Plot)
  
}


p_BA<-tree_sta(BA,"2019-04-03",BA_sta)
p_ON<-tree_sta(ON,"2019-03-26",ON_sta)
p_H3<-tree_sta(H3,"2019-04-30",H3_sta)
p_H1<-tree_sta(H1,"2019-04-25",H1_sta)
p_Vic<-tree_sta(Vic,"2019-04-30",Vic_sta)
p_Yam<-tree_sta(Yam,"2019-04-30",Yam_sta)

plot_grid(p_ON, p_BA, p_H1,p_H3,p_Vic,p_Yam,
          nrow = 3,
          labels = c("RSVA:ON","RSVB:BA","FluA:H1","FluA:H3","FluB:Victoria","FluB:Yamagata"), label_size = 8)->tree_pact

##########################################################################################################################
## Skyride plot
###########################################################################################################################
## Skyride visulization
## RSV skyride
ON_skyride<-read.table('Genetic_analysis/RSV_flu_Beast/RSV/gisaid_ON1_north_100901_190430_skyride.txt', skip = 1, header =TRUE, sep ='\t')
BA_skyride<-read.table('Genetic_analysis/RSV_flu_Beast/RSV/gisaid_BA_north_100901_190430_skyride.txt', skip = 1, header =TRUE, sep ='\t')

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
  ylab("Effective number of\nRSV infections")+
  xlab("Time")+
  scale_x_continuous(breaks = seq(2010,2019,2),position = 'bottom',limits =c(2009,2019) )+
  scale_colour_manual(name="",values = col2)+
  theme(legend.position=c(0.1,0.8),
        text = element_text(size = 8),
        plot.margin = margin(t=1,b = 1,,l=1,r=1)
        )->skyride_rsv

## fluA skyride
H3_skyride<-read.table('Genetic_analysis/RSV_flu_Beast/Flu/gisaid_H3_100901_190401_skyride.txt', skip = 1, header =TRUE, sep ='\t')
H1_skyride<-read.table('Genetic_analysis/RSV_flu_Beast/Flu/gisaid_H1_100901_190430_skyride.txt', skip = 1, header =TRUE, sep ='\t')

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
  scale_x_continuous(breaks = seq(2010,2019,2),position = 'bottom',limits =c(2009,2019) )+
  ylab("Effective number of\nFluA infections")+
  xlab("Time")+
  scale_colour_manual(name="",values = col2)+
  theme(legend.position=c(0.1,0.8),
        text = element_text(size = 8),
        plot.margin = margin(t=1,b = 1,,l=1,r=1),
        legend.key = element_rect(fill = "transparent"))->skyride_fluA


## fluB skyride
Vic_skyride<-read.table('Genetic_analysis/RSV_flu_Beast/Flu/gisaid_Vic_100901_190430_skyride.txt', skip = 1, header =TRUE, sep ='\t')
Yam_skyride<-read.table('Genetic_analysis/RSV_flu_Beast/Flu/gisaid_Yam_100901_190430_skyride.txt', skip = 1, header =TRUE, sep ='\t')

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
  scale_x_continuous(breaks = seq(2010,2019,2),position = 'bottom',limits =c(2009,2019))+
  ylab("Effective number of\nFluB infections")+
  xlab("Time")+
  scale_colour_manual(name="",values = col2)+
  theme(legend.position=c(0.1,0.8),
        text = element_text(size = 8),
        plot.margin = margin(t=1,b = 1,,l=1,r=1))->skyride_fluB


## plot surveliance

RSV_flu_national<-read_csv("Description/RSV_flu_combine_national.csv")

RSV_flu_national%>%
  drop_na()%>%
  dplyr::select(date,RSVpos,fluApos_national,fluBpos_national)%>%
  melt(id.vars = "date",
       variable.name = "type",
       value.name = "case")%>%
  ggplot(aes(x=date, y=case, fill= type)) +
  geom_area() +
  scale_fill_manual(name = "",values=c("#E6AB02","#E7298A","#66A61E"), labels = c("FluA", "FluB", "RSV"))+
  scale_x_continuous(breaks = seq(2010,2019,2),position = 'bottom',limits =c(2009,2019))+
  ylab("Flu and RSV in US")+
  xlab("Time")+
  gg.theme+
  theme(legend.position=c(0.1,0.8),
        text = element_text(size = 8),
        plot.margin = margin(t=1,b = 1,,l=1,r=1))->survel_RSV_fluAB
  


p_skyride<-cowplot::plot_grid(survel_RSV_fluAB,skyride_rsv,skyride_fluA,skyride_fluB, ncol=1, align='v',
                         rel_heights = c(1,1,1,1),
                         labels = c("A","B","",""), label_size = 12)

p_skyrdie_pact<-cowplot::plot_grid(p_skyride,tree_pact,
                                   rel_widths = c(1,1.5),
                                   labels=c("","C"),label_size = 12)

















#########################################
## PACT statistic
BA_sti<-read.csv("Genetic_data/flu_RSV_15-19_north_hemisphere/rsv_10_19_global/rsv_beast/rsv_beast2/gisaid_BA_north_100901_190401_G/PACT/out.stats",sep="\t")
ON_sti<-read.csv("Genetic_data/flu_RSV_15-19_north_hemisphere/rsv_10_19_global/rsv_beast/rsv_beast2/gisaid_ON1_north_100901_190430_north_G/PACT/out.stats",sep="\t")
H3_sti<-read.csv("Genetic_data/flu_RSV_15-19_north_hemisphere/rsv_10_19_global/flu_beast2/gisaid_H3_100901_190430/PACT/out.stats",sep="\t")
H1_sti<-read.csv("Genetic_data/flu_RSV_15-19_north_hemisphere/rsv_10_19_global/flu_beast2/gisaid_H1_100901_190430/PACT/out.stats",sep="\t")
Vic_sti<-read.csv("Genetic_data/flu_RSV_15-19_north_hemisphere/rsv_10_19_global/flu_beast2/gisaid_Vic_100901_190430/PACT/out.stats",sep="\t")
Yam_sti<-read.csv("Genetic_data/flu_RSV_15-19_north_hemisphere/rsv_10_19_global/flu_beast2/gisaid_Yam_100901_190430/PACT/out.stats",sep="\t")


BA_sti%>%
  mutate(pathogen="BA")-> BA_sti
ON_sti%>%
  mutate(pathogen="ON")->ON_sti
H3_sti%>%
  mutate(pathogen="H3")->H3_sti
H1_sti%>%
  mutate(pathogen="H1")->H1_sti
Vic_sti%>%
  mutate(pathogen="Victoria")->Vic_sti
Yam_sti%>%
  mutate(pathogen="Yamagata")->Yam_sti
sta_names<-c("tmrca"="TMRCA","div"="Diversity","coal"="Coalescent rate","tajimad"="Tajima's D")
BA_sti%>%
  rbind(ON_sti,H3_sti,H1_sti,Vic_sti,Yam_sti)%>%
  #filter(statistic=="coal"| statistic=="div")%>%
  filter(statistic!="tmrca")%>%
  mutate(pathogen=factor(pathogen,level=c("ON", "BA","H1","H3","Victoria","Yamagata")))%>%
  ggplot(aes(x=pathogen,y=mean))+
  geom_point(shape=3)+
  geom_errorbar(aes(ymin=lower,ymax=upper))+
  facet_wrap(~statistic,scales = "free_y",ncol=1,labeller = as_labeller(sta_names))+
  theme_bw()->pact_mean

ggarrange(tree_pact,pact_mean,ncol=2,widths=c(2.5,1))
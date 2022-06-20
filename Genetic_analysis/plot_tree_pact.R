library(ggplot2)
library(ggpubr)
library(ggtree)
library(treeio)
setwd("/Users/jianichen/Dropbox/RSV_flu/RSV_flu_git")
##
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
        text = element_text(size = 8))+
  scale_x_continuous(breaks = seq(2010,2019,1),position = 'bottom',limits =c(2008,2019))

  
  
  
  
  
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
  theme(panel.grid.major.x=element_line(color="grey20", linetype="dotted", size=.3))+
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
    theme(panel.grid.major.x=element_line(color="grey20", linetype="dotted", size=.3))+
    scale_x_continuous(breaks = seq(2010,2019,1),position = 'bottom',limits =c(2008,2019))+
    theme(text = element_text(size = 8))
  
  sta%>%
    filter(statistic=="tmrca")%>%
    ggplot(aes(y=mean, x=time))+
    geom_line() + 
    geom_ribbon(aes(ymin=lower, ymax=upper, alpha = 0.2))+
    #facet_wrap(~statistic,ncol=1,strip.position="left",labeller = as_labeller(sta_names))+
    theme_bw()+
    theme(legend.position = "none")+
    scale_x_continuous(breaks = seq(2010,2019,1),position = 'bottom',limits =c(2008,2019))+
    theme(panel.grid.major.x=element_line(color="grey20", linetype="dotted", size=.3))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "none")+
    ylab("TMRCA")+
    ylim(0,10)+
    theme(text = element_text(size = 8))->p2
  
  sta%>%
    filter(statistic=="div")%>%
    ggplot(aes(y=mean, x=time))+
    geom_line() + 
    geom_ribbon(aes(ymin=lower, ymax=upper, alpha = 0.2))+
    #facet_wrap(~statistic,ncol=1,strip.position="left",labeller = as_labeller(sta_names))+
    theme_bw()+
    theme(legend.position = "none")+
    scale_x_continuous(breaks = seq(2010,2019,1),position = 'bottom',limits =c(2008,2019))+
    theme(panel.grid.major.x=element_line(color="grey20", linetype="dotted", size=.3))+
    ylab("Diversity")+
    ylim(0,15)+
    theme(text = element_text(size = 8))->p3
  
  
  Plot<-cowplot::plot_grid(p1,p2,p3, ncol=1, align='v',rel_heights = c(2,1,1.4))
  return(Plot)
  
}


p_BA<-tree_sta(BA,"2019-04-03",BA_sta)
p_ON<-tree_sta(ON,"2019-03-26",ON_sta)
p_H3<-tree_sta(H3,"2019-04-30",H3_sta)
p_H1<-tree_sta(H1,"2019-04-25",H1_sta)
p_Vic<-tree_sta(Vic,"2019-04-30",Vic_sta)
p_Yam<-tree_sta(Yam,"2019-04-30",Yam_sta)


ggarrange( p_ON, p_BA, p_H1,p_H3,p_Vic,p_Yam, ncol = 2, nrow = 3, labels=c("ON","BA","H1","H3","Victoria","Yamagata"))->tree_pact

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
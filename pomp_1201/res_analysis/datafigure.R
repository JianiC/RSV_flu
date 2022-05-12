## surveliance data correlation
source("./fit_functions.R", chdir = TRUE) 
library(ggpubr)
library("lubridate")
RSV_flu_surveliance <-read.csv("RSV_flu_combine.csv")

pop<-read.csv("../../demographic-data/nst-est2019-alldata.csv")
region_state<-read.csv("../../demographic-data/state_region.csv")

RSV_flu_national<-read_csv("RSV_flu_combine_national.csv")

pop%>%
  left_join(region_state,by=c("US.STATE"="State_name"))%>%
  group_by(HHS_region)%>%
  summarise(across(6:27, sum))%>%
  filter(HHS_region !=0)%>%
  select(HHS_region,starts_with("POPEST"))%>%
  gather(year,pop,POPESTIMATE2010:POPESTIMATE2019,factor_key=TRUE)%>%
  mutate(year=as.numeric(gsub("^POPESTIMATE", "",year)))->pop_hhsregion

RSV_flu_surveliance%>% 
  mutate(year=year(RepWeekDate.y))%>%
  right_join(pop_hhsregion, by=c("HHS_REGION"="HHS_region","year"="year"))%>%
  mutate(RSVpos_norm=(RSVpos/pop)*100000,
         fluApos_norm=(fluApos/pop)*100000,
         fluBpos_norm=(fluBpos/pop)*100000)%>%
  select(RSVpos_norm,fluApos_norm,fluBpos_norm,HHS_REGION,date)->RSV_flu_surveliance_norm

## check normal distribution
for (i in 1:10){
  HHSdata<-RSV_flu_surveliance_norm%>%filter(as.numeric(HHS_REGION)==i)
  print(shapiro.test(HHSdata$RSVpos_norm))
  print(shapiro.test(HHSdata$fluApos_norm))
  print(shapiro.test(HHSdata$fluBpos_norm))
}


###################################################################

HHS_corre<-function(data=RSV_flu_surveliance_norm,virus,HHSreg){
  HHSdata<-data%>%filter(as.numeric(HHS_REGION)==HHSreg)
  if(virus=="RSV_fluA"){
    res<-cor.test(HHSdata$RSVpos,HHSdata$fluApos_norm,method="pearson")
  }else if(virus=="RSV_fluB"){
    res<-cor.test(HHSdata$RSVpos,HHSdata$fluBpos_norm,method="pearson")
  }else{
    res<-cor.test(HHSdata$fluApos,HHSdata$fluBpos_norm,method="pearson")
  }
  
  result=data_frame(HHS_region=HHSreg,
                    corre=res$estimate,
                    CI_low=res$conf.int[1],
                    CI_high=res$conf.int[2],
                    virus=virus
                    )
  return(result)
  
}
RSV_fluA_cor<-data.frame()
for( i in 1:10){
  out<-HHS_corre(virus="RSV_fluA",HHSreg=i)
  RSV_fluA_cor<-rbind(RSV_fluA_cor,out)
}

RSV_fluB_cor<-data.frame()
for( i in 1:10){
  out<-HHS_corre(virus="RSV_fluB",HHSreg=i)
  RSV_fluB_cor<-rbind(RSV_fluB_cor,out)
}

fluA_fluB_cor<-data.frame()
for( i in 1:10){
  out<-HHS_corre(virus="fluA_fluB",HHSreg=i)
  fluA_fluB_cor<-rbind(fluA_fluB_cor,out)
}

####################################################################################33

library(ggplot2)
library(dplyr)
library(cowplot)
library(RColorBrewer)
#install.packages("mapproj")
library(mapproj)

# read in state map information
all_states <- map_data("state")
names(all_states)[5] <- "State"
all_states$State = toupper(all_states$State)

#read in region breakdown  (gsub reformts column to match the state column of map data)
HHSreg <- read.csv("HHS-regions.csv", header = TRUE )
HHSreg$State <- toupper(HHSreg$State)
HHSreg$State <- gsub("_"," ", as.character(HHSreg$State))

regiondata <- left_join(all_states,HHSreg,by="State")

regiondata%>%
  dplyr::group_by(HHS_REGION)%>%
  summarise(centroid_long=mean(long),
            centroid_lat=mean(lat))%>%
  full_join(RSV_fluA_cor,by=c("HHS_REGION"="HHS_region"))%>%
  mutate(label=paste("(",round(CI_low,2),",",round(CI_high,2),")",sep=""))-> fluA_label
  
regiondata%>%
  left_join(RSV_fluA_cor,by=c("HHS_REGION"="HHS_region"))%>%
  ggplot()+
  geom_polygon( aes(x=long, y=lat, group = group, fill=corre)) + 
  coord_map() +
  theme(axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+
  theme(legend.position="top")+
  theme(text = element_text(size = 8))+
  theme(aspect.ratio = 0.5, 
        legend.position = "top") +
  guides(linetype = guide_legend(ncol = 1, order = 1), 
         colour = guide_legend(ncol = 3)) +
  theme(legend.spacing.y = unit(0.1, "lines"),
        legend.key = element_blank(), 
        legend.background = element_blank(),
        strip.text.y = element_blank())+
  geom_text(data=fluA_label,aes(centroid_long,centroid_lat,label=label),color="black",size=8)+
  #labs(fill="correlation coefficient of RSV and FluA weekly case")+
  scale_fill_gradient("Correation coefficient",limits=c(0,1))->p_RSVfluA_cor




###########RSV-fluB
regiondata%>%
  dplyr::group_by(HHS_REGION)%>%
  summarise(centroid_long=mean(long),
            centroid_lat=mean(lat))%>%
  full_join(RSV_fluB_cor,by=c("HHS_REGION"="HHS_region"))%>%
  mutate(label=paste("(",round(CI_low,2),",",round(CI_high,2),")",sep=""))-> fluB_label

regiondata%>%
  left_join(RSV_fluB_cor,by=c("HHS_REGION"="HHS_region"))%>%
  ggplot()+
  geom_polygon( aes(x=long, y=lat, group = group, fill=corre)) + 
  coord_map() +
  theme(axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+
  theme(legend.position="top")+
  geom_text(data=fluB_label,aes(centroid_long,centroid_lat,label=label),color="black",size=3)+
  #labs(fill="correlation coefficient of RSV and fluB weekly case")+
  scale_fill_gradient("Correation coefficient",limits=c(0,1))->p_RSVfluB_cor


regiondata%>%
  dplyr::group_by(HHS_REGION)%>%
  summarise(centroid_long=mean(long),
            centroid_lat=mean(lat))%>%
  full_join(fluA_fluB_cor,by=c("HHS_REGION"="HHS_region"))%>%
  mutate(label=paste("(",round(CI_low,2),",",round(CI_high,2),")",sep=""))-> fluAB_label


regiondata%>%
  left_join(fluA_fluB_cor,by=c("HHS_REGION"="HHS_region"))%>%
  ggplot()+
  geom_polygon( aes(x=long, y=lat, group = group, fill=corre)) + 
  coord_map() +
  theme(axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+
  theme(legend.position="top")+
  geom_text(data=fluAB_label,aes(centroid_long,centroid_lat,label=label),color="black",size=3)+
  #labs(fill="correlation coefficient of fluA and fluB weekly case")+
  scale_fill_gradient("Correation coefficient",limits=c(0,1))->p_fluAfluB_cor


ggarrange(p_fluAfluB_cor,p_RSVfluA_cor,p_RSVfluB_cor,
          #labels=c("correlation coefficient of FluA and FluB weekly case",
           #        "correlation coefficient of RSV and FluA weekly case",
            #       "correlation coefficient of RSV and FluB weekly case"),
          labels=c("A","B","C"),
          ncol=1,
          common.legend = TRUE
  
)->corre_map

###################################################################
## facet to plot map together
cor_labs<-c("FluA & FluB","RSV & FluA","RSV & FluB")
names(cor_labs)<-c("fluA_fluB","RSV_fluA","RSV_fluB")


regiondata%>%
  left_join(rbind(RSV_fluA_cor,RSV_fluB_cor,fluA_fluB_cor),by=c("HHS_REGION"="HHS_region"))%>%
  ggplot()+
  geom_polygon( aes(x=long, y=lat, group = group, fill=corre)) + 
  coord_map() +
  theme(axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+
  theme(legend.position="top")+
  facet_wrap(~virus,labeller = labeller(virus=cor_labs))+
  theme(text = element_text(size = 8))+
  theme(aspect.ratio = 0.5, 
        legend.position = "bottom") +
  guides(linetype = guide_legend(ncol = 1, order = 1), 
         colour = guide_legend(ncol = 3)) +
  theme(legend.spacing.y = unit(0.1, "lines"),
        legend.key = element_blank(), 
        legend.background = element_blank(),
        strip.text.y = element_blank())+
  scale_fill_gradient("Correation coefficient",limits=c(0,1))->cor_map
  
  



















## propotion heatmap

RSV_flu_surveliance %>% 
  dplyr::group_by(date,TestType) %>%
  dplyr::summarise(RSVpos=sum(RSVpos),
                   fluApos=sum(fluApos),
                   fluBpos=sum(fluBpos))%>%
  mutate(HHS_REGION="U.S.")%>%
  mutate(fluA_p=fluApos/(RSVpos+fluApos+fluBpos),
                                        RSV_p=RSVpos/(RSVpos+fluApos+fluBpos),
                                        fluB_p=fluBpos/(RSVpos+fluApos+fluBpos))%>%
  select(HHS_REGION,date,RSV_p,fluA_p,fluB_p)%>%
  as.data.frame()->RSV_flu_pro_national

order<-c("U.S.","1","2","3","4","5","6","7","8","9","10")
  RSV_flu_surveliance %>% 
  drop_na()%>%
  mutate(fluA_p=fluApos/(RSVpos+fluApos+fluBpos),
         RSV_p=RSVpos/(RSVpos+fluApos+fluBpos),
         fluB_p=fluBpos/(RSVpos+fluApos+fluBpos))%>%
  select(HHS_REGION,date,RSV_p,fluA_p,fluB_p)%>%
  mutate(HHS_REGION=as.character(HHS_REGION))%>%
    rbind(RSV_flu_pro_national)%>%
  mutate(HHS_REGION=factor(HHS_REGION,levels=order))%>%
  gather(pathogen,percentage,RSV_p:fluB_p,factor_key=TRUE)%>%
  ggplot(aes(x=date,y=percentage,fill=pathogen))+
  geom_area()+
  facet_grid(HHS_REGION~.)+
    theme(text = element_text(size = 8))+
    theme( 
          legend.position = "top") +
    guides(linetype = guide_legend(ncol = 1, order = 1), 
           colour = guide_legend(ncol = 3)) +
    theme(legend.spacing.y = unit(0.1, "lines"),
          legend.key = element_blank(), 
          legend.background = element_blank())+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0,1,0.4)
                     )+
  scale_fill_manual(name = "",values=c("#66A61E","#E6AB02","#E7298A"), labels = c("RSV","FluA", "FluB" ))+
  scale_x_continuous(breaks = seq(2011,2019,1),position = 'bottom',limits =c(2011,2019))+
  xlab("Date")+
  ylab("Percentage")->surveliance_propotion
  
  
  
plot_grid( surveliance_propotion, cor_map,
          
            nrow = 2,
            labels = "AUTO", label_size = 12,
          rel_heights = c(2.5,1))->data_figure 
  
  

ggarrange(surveliance_propotion,ncol=2,labels="A",
          ggarrange(p_fluAfluB_cor,p_RSVfluA_cor,p_RSVfluB_cor,
                    #labels=c("correlation coefficient of FluA and FluB weekly case",
                    #        "correlation coefficient of RSV and FluA weekly case",
                    #       "correlation coefficient of RSV and FluB weekly case"),
                    labels=c("B","C","D"),
                    ncol=1,
                    common.legend = TRUE
                    
)
          
)
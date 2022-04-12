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


## load surveliance data

load("./inc_data_add.rds")
inc_data_add%>%
  drop_na()%>%
  filter(virus=="fluA")%>%
  mutate(season = case_when(
  between(date, 2010.5, 2011.5) ~ "10-11",
  between(date, 2011.5, 2012.5) ~ "11-12",
  between(date, 2012.5, 2013.5) ~ "12-13",
  between(date, 2013.5, 2014.5) ~ "13-14",
  between(date, 2014.5, 2015.5) ~ "14-15",
  between(date, 2015.5, 2016.5) ~ "15-16",
  between(date, 2016.5, 2017.5) ~ "16-17",
  between(date, 2017.5, 2018.5) ~ "17-18",
  between(date, 2018.5, 2019.5) ~ "18-19",
  between(date, 2019.5, 2020.5) ~ "19-20",
  TRUE ~ NA_character_
))%>%
  dplyr::group_by(season,HHS_REGION)%>%
  dplyr::summarise(maxcase=max(cases))%>%
  spread(season,maxcase)%>%
  mutate(peakdiff=`14-15`-((`13-14`+`15-16`+`16-17`)/3))%>%
  dplyr::group_by(HHS_REGION) %>%
  mutate(SD = sd(unlist(select(cur_data(), `13-14`:`16-17`))))->fluA_peakstat

inc_data_add%>%
  drop_na()%>%
  filter(virus=="fluB")%>%
  mutate(season = case_when(
    between(date, 2010.5, 2011.5) ~ "10-11",
    between(date, 2011.5, 2012.5) ~ "11-12",
    between(date, 2012.5, 2013.5) ~ "12-13",
    between(date, 2013.5, 2014.5) ~ "13-14",
    between(date, 2014.5, 2015.5) ~ "14-15",
    between(date, 2015.5, 2016.5) ~ "15-16",
    between(date, 2016.5, 2017.5) ~ "16-17",
    between(date, 2017.5, 2018.5) ~ "17-18",
    between(date, 2018.5, 2019.5) ~ "18-19",
    between(date, 2019.5, 2020.5) ~ "19-20",
    TRUE ~ NA_character_
  ))%>%
  dplyr::group_by(season,HHS_REGION)%>%
  summarise(maxcase=max(cases))%>%
  spread(season,maxcase)%>%
  mutate(peakdiff=`14-15`-((`13-14`+`15-16`+`16-17`)/3))%>%
  group_by(HHS_REGION) %>%
  mutate(SD = sd(unlist(select(cur_data(), `13-14`:`16-17`))))->fluB_peakstat

inc_data_add%>%
  drop_na()%>%
  filter(virus=="RSV")%>%
  mutate(season = case_when(
    between(date, 2010.5, 2011.5) ~ "10-11",
    between(date, 2011.5, 2012.5) ~ "11-12",
    between(date, 2012.5, 2013.5) ~ "12-13",
    between(date, 2013.5, 2014.5) ~ "13-14",
    between(date, 2014.5, 2015.5) ~ "14-15",
    between(date, 2015.5, 2016.5) ~ "15-16",
    between(date, 2016.5, 2017.5) ~ "16-17",
    between(date, 2017.5, 2018.5) ~ "17-18",
    between(date, 2018.5, 2019.5) ~ "18-19",
    between(date, 2019.5, 2020.5) ~ "19-20",
    TRUE ~ NA_character_
  ))%>%
  dplyr::group_by(season,HHS_REGION)%>%
  summarise(maxcase=max(cases))%>%
  spread(season,maxcase)%>%
  mutate(peakdiff=`14-15`-((`13-14`+`15-16`+`16-17`)/3))%>%
  group_by(HHS_REGION) %>%
  mutate(SD = sd(unlist(select(cur_data(), `13-14`:`16-17`))))->RSV_peakstat

label<-read.csv("map_model_label.csv")
regiondata%>%
  dplyr::group_by(HHS_REGION)%>%
  summarise(centroid_long=mean(long),
            centroid_lat=mean(lat))%>%
  full_join(label,by=c("HHS_REGION"="HHS_REGION"))->map_label




#plots
regiondata%>%
  left_join(fluA_peakstat,by=c("HHS_REGION"="HHS_REGION"))%>%
  ggplot()+
  geom_polygon( aes(x=long, y=lat, group = group, fill=peakdiff)
               ,color="grey50") + coord_map() +
  theme(axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(legend.position="top")+
  geom_text(data=map_label,aes(centroid_long,centroid_lat,label=fluA_model),color="white")+
  labs(fill="Flu A peak differences (case)")->p_fluA_peakstat

  
regiondata%>%
  left_join(RSV_peakstat,by=c("HHS_REGION"="HHS_REGION"))%>%
  ggplot()+
  geom_polygon( aes(x=long, y=lat, group = group, fill=peakdiff)
                ,color="grey50") + coord_map() +
  theme(axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(legend.position="top")+
  geom_text(data=map_label,aes(centroid_long,centroid_lat,label=fluA_model),color="white")+
  labs(fill="RSV peak differences (case)")->p_RSVa_peakstat 

regiondata%>%
  left_join(fluA_peakstat,by=c("HHS_REGION"="HHS_REGION"))%>%
  ggplot()+
  geom_polygon( aes(x=long, y=lat, group = group, fill=peakdiff)
                ,color="grey50") + coord_map() +
  theme(axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(legend.position="top")+
  geom_text(data=map_label,aes(centroid_long,centroid_lat,label=fluA_model),color="white")+
  labs(fill="Flu A peak differences (case)")->p_fluA_peakstat


regiondata%>%
  left_join(RSV_peakstat,by=c("HHS_REGION"="HHS_REGION"))%>%
  ggplot()+
  geom_polygon( aes(x=long, y=lat, group = group, fill=peakdiff)
                ,color="grey50") + coord_map() +
  theme(axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(legend.position="top")+
  geom_text(data=map_label,aes(centroid_long,centroid_lat,label=fluB_model),color="white")+
  labs(fill="RSV peak differences (case)")->p_RSVb_peakstat 


regiondata%>%
  left_join(fluB_peakstat,by=c("HHS_REGION"="HHS_REGION"))%>%
  ggplot()+
  geom_polygon( aes(x=long, y=lat, group = group, fill=peakdiff)
                ,color="grey50") + coord_map() +
  theme(axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  theme(legend.position="top")+
  geom_text(data=map_label,aes(centroid_long,centroid_lat,label=fluA_model),color="white")+
  labs(fill="Flu B peak differences (case)")->p_fluB_peakstat

library(ggpubr)
ggarrange(p_RSVa_peakstat,p_fluA_peakstat,
          p_RSVb_peakstat,p_fluB_peakstat)

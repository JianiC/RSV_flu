
source("./fit_functions.R", chdir = TRUE) 
################
## plot with the best-fit model
 #########################################################################################
res_arsv_bestmodel<-list.files(path="pomp_longphi_result/best_model_rsva/",pattern=".rds")
for (i in 1:length(res_arsv_bestmodel)){
  load(paste("pomp_longphi_result/best_model_rsva/",res_arsv_bestmodel[i],sep=""))
  print(paste0("load ",res_arsv_bestmodel[i]))
}

res_arsv_bestmodel_list<-list(res_hhs1_arsv_coinfect,res_hhs2_arsv_chi,res_hhs3_arsv_coinfect,res_hhs4_arsv_coinfect,res_hhs5_arsv_coinfect,
  res_hhs6_arsv_coinfect,res_hhs7_arsv_coinfect,res_hhs8_arsv_neutral,res_hhs9_arsv_chi,res_hhs10_arsv_chi)
est_res_arsv_coinfect<-get_est_all(res_arsv_bestmodel_list)

traj_fit_arsv_bestmodel<-get_traj_fitall(res_list=res_arsv_bestmodel_list)

## clear large list to release memory  
rm(res_hhs1_arsv_coinfect,res_hhs2_arsv_chi,res_hhs3_arsv_coinfect,res_hhs4_arsv_coinfect,res_hhs5_arsv_coinfect,
   res_hhs6_arsv_coinfect,res_hhs7_arsv_coinfect,res_hhs8_arsv_neutral,res_hhs9_arsv_chi,res_hhs10_arsv_chi)

patho_order<-c("RSV","FluA")
inc_data_fit<-inc_data_perdic%>%
  mutate(HHSregion=HHS_REGION)%>%
  #pivot_wider(names_from = virus, values_from = cases)%>%
  drop_na()%>%
  mutate(method=case_when(
    time<=6.5 ~ "Within sample",
    time>6.5 ~ "Out of sample"
  ))%>%
  mutate(virus=case_when(
    virus=="RSV" ~ "RSV",
    virus=="fluA" ~ "FluA",
    virus=="fluB" ~ "FluB"
  ))%>%
  mutate(pathogen=factor(virus,levels=patho_order))%>%
  filter(pathogen!="FluB")



traj_fit_arsv_bestmodel%>%
  mutate(pathogen=factor(type,levels=patho_order))%>%
  mutate(method=case_when(
    time<=6.5 ~ "Within sample",
    time>6.5 ~ "Out of sample"
  ))%>%
  ggplot(aes(x=time+2011,y=cases))+
  geom_line(data=inc_data_fit,aes(x=date,y=cases))+
  geom_line(aes(color=method),size=0.5)+
  geom_ribbon(aes(ymin = cases_lb, ymax = cases_ub, fill = method),alpha=0.5)+
  facet_grid(HHSregion~pathogen,scales="free")+
  theme_classic()+
  theme(legend.position="bottom")+
  xlab("Time")+
  ylab("Cases")+
  scale_color_manual(name="",values = c("Within sample"="#1B9E77","Out of sample"="#D95F02","Observed"="black"))+
  guides(fill = "none")+
  gg.theme +
  theme(aspect.ratio = 0.35, 
        legend.position = "top") +
  guides(linetype = guide_legend(ncol = 1, order = 1), 
         colour = guide_legend(ncol = 3)) +
  theme(legend.spacing.y = unit(0.1, "lines"),
        legend.key = element_blank(), 
        legend.background = element_blank())->fluA_besfit



legend_pos = c(0.1, 0.68)
inc_data_fit_repre<-inc_data_fit%>%filter(HHSregion==1)
traj_fit_arsv_bestmodel%>%
  filter(HHSregion==1)%>%
  mutate(pathogen=factor(type,levels=patho_order))%>%
  mutate(method=case_when(
    time<=6.5 ~ "Within sample",
    time>6.5 ~ "Out of sample"
  ))%>%
  ggplot(aes(x=time+2011,y=cases))+
  geom_line(data=inc_data_fit_repre,aes(x=date,y=cases))+
  geom_line(aes(color=method),size=0.5)+
  geom_ribbon(aes(ymin = cases_lb, ymax = cases_ub, fill = method),alpha=0.5)+
  facet_grid(.~pathogen,scales="free")+
  theme_classic()+
  theme(legend.position="bottom")+
  xlab("Time")+
  ylab("Cases")+
  scale_color_manual(name="",values = c("Within sample"="#1B9E77","Out of sample"="#D95F02","Observed"="black"))+
  guides(fill = "none")+
  gg.theme +
  theme(aspect.ratio = 0.5, 
        legend.position = legend_pos) +
  guides(linetype = guide_legend(ncol = 1, order = 1), 
         colour = guide_legend(ncol = 1)) +
  theme(legend.spacing.y = unit(0.1, "lines"),
        legend.key = element_blank(), 
        legend.background = element_blank())->fluA_besfit_repre
  
#######################################################

## RSV and fluB modleing


res_brsv_bestmodel<-list.files(path="pomp_longphi_result/best_model_rsvb/",pattern=".rds")
for (i in 1:length(res_brsv_bestmodel)){
  load(paste("pomp_longphi_result/best_model_rsvb/",res_brsv_bestmodel[i],sep=""))
  print(paste0("load ",res_brsv_bestmodel[i]))
}

res_brsv_bestmodel_list<-list(
  res_hhs1_brsv_chi,res_hhs2_brsv_chi,res_hhs3_brsv_coinfect,res_hhs4_brsv_coinfect,res_hhs5_brsv_chi,
  res_hhs6_brsv_chi, res_hhs7_brsv_coinfect,
  res_hhs8_brsv_neutral,res_hhs9_brsv_coinfect,res_hhs10_brsv_chi)




traj_fit_brsv_bestmodel<-get_traj_fitall(res_list=res_brsv_bestmodel_list)

## clear large list to release memory  
rm(res_hhs1_brsv_chi,res_hhs2_brsv_chi,res_hhs3_brsv_coinfect,res_hhs4_brsv_coinfect,res_hhs5_brsv_chi,
   res_hhs6_brsv_chi, res_hhs7_brsv_coinfect,
   res_hhs8_brsv_neutral,res_hhs9_brsv_coinfect,res_hhs10_brsv_chi)

patho_order<-c("RSV","FluB")
inc_data_fit<-inc_data_perdic%>%
  mutate(HHSregion=HHS_REGION)%>%
  #pivot_wider(names_from = virus, values_from = cases)%>%
  drop_na()%>%
  mutate(method=case_when(
    time<=6.5 ~ "Within sample",
    time>6.5 ~ "Out of sample"
  ))%>%
  mutate(virus=case_when(
    virus=="RSV" ~ "RSV",
    virus=="fluA" ~ "FluA",
    virus=="fluB" ~ "FluB"
  ))%>%
  mutate(pathogen=factor(virus,levels=patho_order))%>%
  filter(pathogen!="FluA")



traj_fit_brsv_bestmodel%>%
  mutate(pathogen=factor(type,levels=patho_order))%>%
  mutate(method=case_when(
    time<=6.5 ~ "Within sample",
    time>6.5 ~ "Out of sample"
  ))%>%
  ggplot(aes(x=time+2011,y=cases))+
  geom_line(data=inc_data_fit,aes(x=date,y=cases))+
  geom_line(aes(color=method),size=0.5)+
  geom_ribbon(aes(ymin = cases_lb, ymax = cases_ub, fill = method),alpha=0.5)+
  facet_grid(HHSregion~pathogen,scales="free")+
  theme_classic()+
  theme(legend.position="bottom")+
  xlab("Time")+
  ylab("Cases")+
  scale_color_manual(name="",values = c("Within sample"="#1B9E77","Out of sample"="#D95F02","Observed"="black"))+
  guides(fill = "none")+
  gg.theme +
  theme(aspect.ratio = 0.35, 
        legend.position = "top") +
  guides(linetype = guide_legend(ncol = 1, order = 1), 
         colour = guide_legend(ncol = 3)) +
  theme(legend.spacing.y = unit(0.1, "lines"),
        legend.key = element_blank(), 
        legend.background = element_blank())->fluB_besfit

inc_data_fit_repre<-inc_data_fit%>%filter(HHSregion==1)
traj_fit_brsv_bestmodel%>%
  filter(HHSregion==1)%>%
  mutate(pathogen=factor(type,levels=patho_order))%>%
  mutate(method=case_when(
    time<=6.5 ~ "Within sample",
    time>6.5 ~ "Out of sample"
  ))%>%
  ggplot(aes(x=time+2011,y=cases))+
  #geom_text(aes(label = "B", x = -16.5, y = 590))+
  geom_line(data=inc_data_fit_repre,aes(x=date,y=cases))+
  geom_line(aes(color=method),size=0.5)+
  geom_ribbon(aes(ymin = cases_lb, ymax = cases_ub, fill = method),alpha=0.5)+
  facet_grid(.~pathogen,scales="free")+
  theme_classic()+
  theme(legend.position="bottom")+
  xlab("Time")+
  ylab("Cases")+
  scale_color_manual(name="",values = c("Within sample"="#1B9E77","Out of sample"="#D95F02","Observed"="black"))+
  guides(fill = "none",color="none")+
  gg.theme +
  theme(aspect.ratio = 0.5, 
        legend.position = "none") +
  guides(linetype = guide_legend(ncol = 1, order = 1), 
         colour = guide_legend(ncol = 1)) +
  theme(legend.spacing.y = unit(0.1, "lines"),
        legend.key = element_blank(), 
        legend.background = element_blank())->fluB_besfit_repre



library(cowplot)
plot_grid(fluA_besfit_repre, fluB_besfit_repre, nrow = 2,
          labels = "AUTO", label_size = 12) 


legend <- get_legend(
  fluB_besfit + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

plot_grid(fluA_besfit+theme(legend.position="none"), 
          fluB_besfit+ theme(legend.position="none"), 
          nrow = 1,
          labels = "AUTO", label_size = 12)->best_fit

plot_grid(best_fit,legend,
          ncol=1,
          rel_heights = c(1, .1))
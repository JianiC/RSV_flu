res_brsv_coinfect<-list.files(path="pomp_result_1213/res_brsv_coinfect/",pattern=".rds")
for (i in 1:length(res_brsv_coinfect)){
   load(paste("pomp_result_1213/res_brsv_coinfect/",res_brsv_coinfect[i],sep=""))
   print(paste0("load ",res_brsv_coinfect[i]))
}

res_brsv_coinfect_list<-list(
   res_hhs1_brsv_coinfect,res_hhs2_brsv_coinfect,res_hhs3_brsv_coinfect,
   res_hhs4_brsv_coinfect,res_hhs5_brsv_coinfect,res_hhs6_brsv_coinfect,
   res_hhs7_brsv_coinfect,res_hhs8_brsv_coinfect,res_hhs9_brsv_coinfect,res_hhs10_brsv_coinfect)
est_res_brsv_coinfect<-get_est_all(res_brsv_coinfect_list)

traj_fit_brsv_coinfect<-get_traj_fitall(res_brsv_coinfect_list)

## clear large list to release memory  
rm(res_brsv_coinfect_list,res_hhs1_brsv_coinfect,res_hhs2_brsv_coinfect,res_hhs3_brsv_coinfect,
   res_hhs4_brsv_coinfect,res_hhs5_brsv_coinfect,res_hhs6_brsv_coinfect,
   res_hhs7_brsv_coinfect,res_hhs8_brsv_coinfect,res_hhs9_brsv_coinfect,res_hhs10_brsv_coinfect)


#########################################################################################  


res_brsv_neutral<-list.files(path="pomp_result_1213/res_brsv_neutral/",pattern=".rds")
for (i in 1:length(res_brsv_neutral)){
   load(paste("pomp_result_1213/res_brsv_neutral/",res_brsv_neutral[i],sep=""))
   print(paste0("load ",res_brsv_neutral[i]))
}

res_brsv_neutral_list<-list(res_hhs1_brsv_neutral,res_hhs2_brsv_neutral,res_hhs3_brsv_neutral,
                            res_hhs4_brsv_neutral,res_hhs5_brsv_neutral,res_hhs6_brsv_neutral,
                            res_hhs7_brsv_neutral,res_hhs8_brsv_neutral,res_hhs9_brsv_neutral,res_hhs10_brsv_neutral)
est_res_brsv_neutral<-get_est_all(res_brsv_neutral_list)
traj_fit_brsv_neutral<-get_traj_fitall(res_brsv_neutral_list)

rm(res_hhs1_brsv_neutral,res_hhs2_brsv_neutral,res_hhs3_brsv_neutral,
   res_hhs4_brsv_neutral,res_hhs5_brsv_neutral,res_hhs6_brsv_neutral,
   res_hhs7_brsv_neutral,res_hhs8_brsv_neutral,res_hhs9_brsv_neutral,res_hhs10_brsv_neutral,res_brsv_neutral_list)


#########################################################################################   

res_brsv_psi<-list.files(path="pomp_result_1213/res_brsv_psi/",pattern=".rds")
for (i in 1:length(res_brsv_psi)){
   load(paste("pomp_result_1213/res_brsv_psi/",res_brsv_psi[i],sep=""))
   print(paste0("load ",res_brsv_psi[i]))
}

res_brsv_psi_list<-list(res_hhs1_brsv_psi,res_hhs2_brsv_psi,res_hhs3_brsv_psi,
                        res_hhs4_brsv_psi,res_hhs5_brsv_psi,res_hhs6_brsv_psi,
                        res_hhs7_brsv_psi,res_hhs8_brsv_psi,res_hhs9_brsv_psi,res_hhs10_brsv_psi)
est_res_brsv_psi<-get_est_all(res_brsv_psi_list)
traj_fit_brsv_psi<-get_traj_fitall(res_brsv_psi_list)

rm(res_hhs1_brsv_psi,res_hhs2_brsv_psi,res_hhs3_brsv_psi,
   res_hhs4_brsv_psi,res_hhs5_brsv_psi,res_hhs6_brsv_psi,
   res_hhs7_brsv_psi,res_hhs8_brsv_psi,res_hhs9_brsv_psi,res_hhs10_brsv_psi,res_brsv_psi_list)

#########################################################################################   


res_brsv_chi<-list.files(path="pomp_result_1213/res_brsv_chi/",pattern=".rds")
for (i in 1:length(res_brsv_chi)){
   load(paste("pomp_result_1213/res_brsv_chi/",res_brsv_chi[i],sep=""))
   print(paste0("load ",res_brsv_chi[i]))
   
}


res_brsv_chi_list<-list(res_hhs1_brsv_chi,res_hhs2_brsv_chi,res_hhs3_brsv_chi,
                        res_hhs4_brsv_chi,res_hhs5_brsv_chi,res_hhs6_brsv_chi,
                        res_hhs7_brsv_chi,res_hhs8_brsv_chi,res_hhs9_brsv_chi,res_hhs10_brsv_chi)
est_res_brsv_chi<-get_est_all(res_brsv_chi_list)
traj_fit_brsv_chi<-get_traj_fitall(res_brsv_chi_list)

rm(res_hhs1_brsv_chi,res_hhs2_brsv_chi,res_hhs3_brsv_chi,
   res_hhs4_brsv_chi,res_hhs5_brsv_chi,res_hhs6_brsv_chi,
   res_hhs7_brsv_chi,res_hhs8_brsv_chi,res_hhs9_brsv_chi,res_hhs10_brsv_chi,res_brsv_chi_list)
#########################################################################################   
est_res_brsv<-rbind(est_res_brsv_neutral,est_res_brsv_psi,est_res_brsv_chi,est_res_brsv_coinfect)
hypo_levels=c("neutral","psi","chi","co-infect")

est_res_brsv%>%
   mutate(hyphothesis=factor(hyphothesis,levels=hypo_levels))%>%
   ggplot(aes(x=hyphothesis,y=AIC,color=hyphothesis))+
   geom_point(size=3)+
   facet_wrap(~HHSregion,scales="free_y",nrow = 1)+
   theme_bw()+
   theme(legend.position="bottom")+
   scale_colour_brewer(palette = "Dark2")+
   scale_linetype_manual(values=c("solid", "dotted"))+
   theme(legend.position="bottom")+
   #scale_color_discrete(name = "Hypothesis", labels = c("neutral", "inhibition on co-infection", "cross-protection","co-ifection + cross-protection"))+
   scale_colour_brewer(palette = "Dark2",name = "Hypothesis", 
                       labels = c("neutral", "inhibition on co-infection", "cross-protection","co-ifection + cross-protection"))+
   theme(axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank())->brsv_AIC



## plot trajectory fit
inc_data_fit<-inc_data_perdic%>%
   mutate(HHSregion=HHS_REGION)%>%
   pivot_wider(names_from = virus, values_from = cases)%>%
   drop_na()%>%
   mutate(method=case_when(
      time<=6.5 ~ "Fit",
      time>6.5 ~ "predict"
   ))

traj_fit_brsv<-rbind(traj_fit_brsv_neutral,traj_fit_brsv_coinfect,traj_fit_brsv_psi,traj_fit_brsv_chi)
traj_fit_brsv_coinfect%>%
   mutate(hypothesis="co_infect")
hypo_levels=c("neutral","psi","chi","co_infect")
patho_order<-c("RSV","fluB")

inc_data_fit<-inc_data_perdic%>%
   mutate(HHSregion=HHS_REGION)%>%
   pivot_wider(names_from = virus, values_from = cases)%>%
   drop_na()%>%
   mutate(method=case_when(
      time<=6.5 ~ "Model fit",
      time>6.5 ~ "Project"
   ))


   traj_fit_brsv%>%
      mutate(hypothesis=factor(hypothesis,levels=hypo_levels))%>%
      mutate(pathogen=factor(type,levels=patho_order))%>%
      mutate(method=case_when(
         time<=6.5 ~ "Model fit",
         time>6.5 ~ "Project"
      ))%>%
      #filter(HHSregion!=1 | hypothesis!="chi")%>%
      #filter(HHSregion!=6 | hypothesis!="psi")%>%
      #filter(HHSregion!=6 | hypothesis!="chi")%>%
      #filter(HHSregion!=7 | hypothesis!="chi")%>%
      #filter(HHSregion!=7 | hypothesis!="neutral")%>%
      #filter(HHSregion!=6)%>%
      #filter(HHSregion!=7)%>%
      ggplot(aes(x=time+2011,y=cases))+
      geom_line(aes(color=hypothesis,linetype=pathogen ),alpha=0.5)+
      
      geom_area(data=inc_data_fit,aes(x=date,y=fluB),fill="gray",alpha=0.6)+
      geom_area(data=inc_data_fit,aes(x=date,y=RSV),fill="brown",alpha=0.4)+
      facet_grid(HHSregion~method,scales="free_x",space="free")+
      scale_colour_brewer(palette = "Dark2")+
      theme_bw()+
      xlab("time(weeks)")+
      scale_colour_brewer(palette = "Dark2",name = "Hypothesis", 
                          labels = c("neutral", "inhibition on co-infection", "cross-protection","co-ifection + cross-protection"))->brsv_trajfit
#      theme(legend.position="bottom")
      
   
 library(ggpubr)  
   ggarrange(
      arsv_AIC, brsv_AIC, labels = c("A", "B"),
      nrow=2)
   
   
   ggarrange(
      arsv_trajfit, brsv_trajfit, labels = c("A", "B"),
      nrow=2)
   
   
   
   
   
    
   
   
   
load("pomp_result_1213/res_brsv_neutral/res_hhs7_brsv_neutral.rds")
test_traj(res_hhs=res_hhs7_brsv_neutral)%>%
   slice( 2:n()) %>% 
   select( -`.id`) %>% 
   gather( "comp", "count", -time) %>% 
   ggplot( aes(x = time, y = count)) +
   geom_line()+
   facet_wrap(.~comp, scales = "free")


res_hhs7_brsv_neutral$DEobj$optim$bestmem

source("./fit_functions.R", chdir = TRUE) 


res_brsv_coinfect<-list.files(path="pomp_longphi_result/res_brsv_coinfect/",pattern=".rds")
for (i in 1:length(res_brsv_coinfect)){
   load(paste("pomp_longphi_result/res_brsv_coinfect/",res_brsv_coinfect[i],sep=""))
   print(paste0("load ",res_brsv_coinfect[i]))
}

res_brsv_coinfect_list<-list(
   res_hhs1_brsv_coinfect,res_hhs2_brsv_coinfect,res_hhs3_brsv_coinfect,res_hhs4_brsv_coinfect,
   res_hhs5_brsv_coinfect,res_hhs6_brsv_coinfect,
   res_hhs7_brsv_coinfect,res_hhs8_brsv_coinfect,res_hhs9_brsv_coinfect,res_hhs10_brsv_coinfect)
est_res_brsv_coinfect<-get_est_all(res_brsv_coinfect_list)
get_fitmeasure_all(res_brsv_coinfect_list)%>%
   full_join(est_res_brsv_coinfect,
             by = c("HHSregion" =  "HHSregion", "Pathogen2"="pathogen2","Hypothesis" = "hyphothesis"))->result_res_brsv_coinfect


traj_fit_brsv_coinfect<-get_traj_fitall(data=inc_data_add,res_brsv_coinfect_list)

## clear large list to release memory  
rm(res_brsv_coinfect_list,res_hhs1_brsv_coinfect,res_hhs2_brsv_coinfect,res_hhs3_brsv_coinfect,
   res_hhs4_brsv_coinfect,res_hhs5_brsv_coinfect,res_hhs6_brsv_coinfect,
   res_hhs7_brsv_coinfect,res_hhs8_brsv_coinfect,res_hhs9_brsv_coinfect,res_hhs10_brsv_coinfect)


#########################################################################################  


res_brsv_neutral<-list.files(path="pomp_longphi_result/res_brsv_neutral/",pattern=".rds")
for (i in 1:length(res_brsv_neutral)){
   load(paste("pomp_longphi_result/res_brsv_neutral/",res_brsv_neutral[i],sep=""))
   print(paste0("load ",res_brsv_neutral[i]))
}

res_brsv_neutral_list<-list(res_hhs1_brsv_neutral,res_hhs2_brsv_neutral,res_hhs3_brsv_neutral,
                            res_hhs4_brsv_neutral,res_hhs5_brsv_neutral,res_hhs6_brsv_neutral,
                            res_hhs7_brsv_neutral,res_hhs8_brsv_neutral,res_hhs9_brsv_neutral,res_hhs10_brsv_neutral)
est_res_brsv_neutral<-get_est_all(res_brsv_neutral_list)

get_fitmeasure_all(res_brsv_neutral_list)%>%
   full_join(est_res_brsv_neutral,
             by = c("HHSregion" =  "HHSregion", "Pathogen2"="pathogen2","Hypothesis" = "hyphothesis"))->result_res_brsv_neutral
traj_fit_brsv_neutral<-get_traj_fitall(data=inc_data_add,res_brsv_neutral_list)

rm(res_hhs1_brsv_neutral,res_hhs2_brsv_neutral,res_hhs3_brsv_neutral,
   res_hhs4_brsv_neutral,res_hhs5_brsv_neutral,res_hhs6_brsv_neutral,
   res_hhs7_brsv_neutral,res_hhs8_brsv_neutral,res_hhs9_brsv_neutral,res_hhs10_brsv_neutral,res_brsv_neutral_list)


#########################################################################################   

res_brsv_psi<-list.files(path="pomp_longphi_result/res_brsv_psi/",pattern=".rds")
for (i in 1:length(res_brsv_psi)){
   load(paste("pomp_longphi_result/res_brsv_psi/",res_brsv_psi[i],sep=""))
   print(paste0("load ",res_brsv_psi[i]))
}

res_brsv_psi_list<-list(res_hhs1_brsv_psi,res_hhs2_brsv_psi,res_hhs3_brsv_psi,
                        res_hhs4_brsv_psi,res_hhs5_brsv_psi,res_hhs6_brsv_psi,
                        res_hhs7_brsv_psi,res_hhs8_brsv_psi,res_hhs9_brsv_psi,res_hhs10_brsv_psi)
est_res_brsv_psi<-get_est_all(res_brsv_psi_list)
get_fitmeasure_all(res_brsv_psi_list)%>%
   full_join(est_res_brsv_psi,
             by = c("HHSregion" =  "HHSregion", "Pathogen2"="pathogen2","Hypothesis" = "hyphothesis"))->result_res_brsv_psi

traj_fit_brsv_psi<-get_traj_fitall(data=inc_data_add,res_brsv_psi_list)

rm(res_hhs1_brsv_psi,res_hhs2_brsv_psi,res_hhs3_brsv_psi,
   res_hhs4_brsv_psi,res_hhs5_brsv_psi,res_hhs6_brsv_psi,
   res_hhs7_brsv_psi,res_hhs8_brsv_psi,res_hhs9_brsv_psi,res_hhs10_brsv_psi,res_brsv_psi_list)

#########################################################################################   


res_brsv_chi<-list.files(path="pomp_longphi_result/res_brsv_chi/",pattern=".rds")
for (i in 1:length(res_brsv_chi)){
   load(paste("pomp_longphi_result/res_brsv_chi/",res_brsv_chi[i],sep=""))
   print(paste0("load ",res_brsv_chi[i]))
   
}


res_brsv_chi_list<-list(res_hhs1_brsv_chi,res_hhs2_brsv_chi,res_hhs3_brsv_chi,
                        res_hhs4_brsv_chi,res_hhs5_brsv_chi,res_hhs6_brsv_chi,
                        res_hhs7_brsv_chi,res_hhs8_brsv_chi,res_hhs9_brsv_chi,res_hhs10_brsv_chi)
est_res_brsv_chi<-get_est_all(res_brsv_chi_list)
get_fitmeasure_all(res_brsv_chi_list)%>%
   full_join(est_res_brsv_chi,
             by = c("HHSregion" =  "HHSregion", "Pathogen2"="pathogen2","Hypothesis" = "hyphothesis"))->result_res_brsv_chi
traj_fit_brsv_chi<-get_traj_fitall(data=inc_data_add,res_brsv_chi_list)

rm(res_hhs1_brsv_chi,res_hhs2_brsv_chi,res_hhs3_brsv_chi,
   res_hhs4_brsv_chi,res_hhs5_brsv_chi,res_hhs6_brsv_chi,
   res_hhs7_brsv_chi,res_hhs8_brsv_chi,res_hhs9_brsv_chi,res_hhs10_brsv_chi,res_brsv_chi_list)
#########################################################################################   
#est_res_brsv<-rbind(est_res_brsv_neutral,est_res_brsv_psi,est_res_brsv_chi,est_res_brsv_coinfect)

rbind(result_res_brsv_neutral,result_res_brsv_psi,result_res_brsv_chi,result_res_brsv_coinfect)%>%
   select("Pathogen2","HHSregion","Hypothesis",
          "R01","R02","amplitude1","amplitude2","tpeak1","tpeak2","rho1","rho2", "psi","chi",
          "loglik","AIC","R2","RMSE" )%>%
   mutate(psi=1-psi,
          chi=1-chi)%>%
   melt(id.vars=c("Pathogen2","HHSregion","Hypothesis"))%>%
   dcast(Pathogen2 + HHSregion+variable~ Hypothesis)%>%
   select(Pathogen2, HHSregion, variable, neutral,psi,chi,co_infect) -> result_res_brsv

write.csv(result_res_brsv,"./figures/result_res_brsv_longphhi_check.csv",row.names = FALSE)


result_res_brsv<-read.csv("./figures/result_res_brsv_longphhi.csv")

result_res_brsv%>%
   filter(variable=="AIC")%>%
   gather(hypothesis, AIC, neutral:co_infect, factor_key=TRUE)%>%
   group_by(HHSregion)%>%
   summarise(min_AIC=min(AIC))->brsv_minAIC

result_res_brsv%>%
   filter(variable=="AIC")%>%
   gather(hypothesis, AIC, neutral:co_infect, factor_key=TRUE)%>%
   left_join(brsv_minAIC,by=c("HHSregion"="HHSregion"))%>%
   mutate(delta_AIC=AIC-min_AIC)%>%
   select(Pathogen2,HHSregion,hypothesis,delta_AIC)%>%
   spread(hypothesis,delta_AIC)%>%
   mutate(variable="deltaAIC")->brsv_deltaAIC

result_res_brsv%>%
   filter(variable=="R2"| variable=="RMSE")%>%
   rbind(brsv_deltaAIC)->brsv_modelcompare

write.csv(brsv_modelcompare,"./figures/brsv_modelcompare.csv",row.names = FALSE)



hypo_levels=c("neutral","psi","chi","co_infect")



 




#      theme(legend.position="bottom")
      
   
 library(ggpubr)  
   ggarrange(
      arsv_AIC, brsv_AIC, labels = c("A", "B"),
      nrow=2)
   
   
   ggarrange(
      arsv_trajfit, brsv_trajfit, labels = c("A", "B"),
      nrow=2)
   
   
   
   
   
    
   
   
   




   traj_fit_brsv%>%
      filter(HHSregion!=1 | hypothesis!="chi")%>%
      filter(HHSregion!=6 | hypothesis!="psi")%>%
      filter(HHSregion!=6 | hypothesis!="chi")%>%
      filter(HHSregion!=7 | hypothesis!="chi")%>%
      filter(HHSregion!=7 | hypothesis!="neutral")%>%
      #filter(HHSregion!=6)%>%
      #filter(HHSregion!=7)%>%
      ggplot(aes(x=time+2011,y=cases))+
      geom_line(aes(color=hypothesis,linetype=type))+
      
      geom_area(data=inc_data_fit,aes(x=date,y=fluB),fill="gray",alpha=0.4)+
      geom_area(data=inc_data_fit,aes(x=date,y=RSV),fill="brown",alpha=0.4)+
      facet_grid(HHSregion~.)+
      scale_colour_brewer(palette = "Dark2")+
      theme_bw()

load("pomp_result_1213/res_brsv_neutral/res_hhs7_brsv_neutral.rds")
test_traj(res_hhs=res_hhs7_brsv_neutral)%>%
   slice( 2:n()) %>% 
   select( -`.id`) %>% 
   gather( "comp", "count", -time) %>% 
   ggplot( aes(x = time, y = count)) +
   geom_line()+
   facet_wrap(.~comp, scales = "free")


res_hhs7_brsv_neutral$DEobj$optim$bestmem
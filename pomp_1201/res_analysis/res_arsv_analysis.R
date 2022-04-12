
#Sys.setenv('R_MAX_VSIZE'=64000000000)
#Sys.getenv('R_MAX_VSIZE')
Sys.setenv('R_MAX_NUM_DLL'=1000)
Sys.getenv('R_MAX_NUM_DLL')
system("ulimit -n 1000")
# first load prerequisites
source("./fit_functions.R", chdir = TRUE) 
#########################################################################################
  res_arsv_coinfect<-list.files(path="pomp_longphi_result/res_arsv_coinfect/",pattern=".rds")
  for (i in 1:length(res_arsv_coinfect)){
    load(paste("pomp_longphi_result/res_arsv_coinfect/",res_arsv_coinfect[i],sep=""))
    print(paste0("load ",res_arsv_coinfect[i]))
  }
  
  res_arsv_coinfect_list<-list(
    res_hhs1_arsv_coinfect,res_hhs2_arsv_coinfect,res_hhs3_arsv_coinfect,
    res_hhs4_arsv_coinfect,res_hhs5_arsv_coinfect,res_hhs6_arsv_coinfect,
    res_hhs7_arsv_coinfect,res_hhs8_arsv_coinfect,res_hhs9_arsv_coinfect,res_hhs10_arsv_coinfect)
  est_res_arsv_coinfect<-get_est_all(res_arsv_coinfect_list)
  get_fitmeasure_all(res_arsv_coinfect_list)%>%
    full_join(est_res_arsv_coinfect,
              by = c("HHSregion" =  "HHSregion", "Pathogen2"="pathogen2","Hypothesis" = "hyphothesis"))->result_res_arsv_coinfect

  
  traj_fit_arsv_coinfect<-get_traj_fitall(data=inc_data_add,res_arsv_coinfect_list)
  
## clear large list to release memory  
  rm(res_arsv_coinfect_list,res_hhs1_arsv_coinfect,res_hhs2_arsv_coinfect,res_hhs3_arsv_coinfect,
     res_hhs4_arsv_coinfect,res_hhs5_arsv_coinfect,res_hhs6_arsv_coinfect,
     res_hhs7_arsv_coinfect,res_hhs8_arsv_coinfect,res_hhs9_arsv_coinfect,res_hhs10_arsv_coinfect)

  
  #########################################################################################  
  
  
  res_arsv_neutral<-list.files(path="pomp_longphi_result/res_arsv_neutral/",pattern=".rds")
  for (i in 1:length(res_arsv_neutral)){
    load(paste("pomp_longphi_result/res_arsv_neutral/",res_arsv_neutral[i],sep=""))
    print(paste0("load ",res_arsv_neutral[i]))
  }
  
  res_arsv_neutral_list<-list(res_hhs1_arsv_neutral,res_hhs2_arsv_neutral,res_hhs3_arsv_neutral,
                              res_hhs4_arsv_neutral,res_hhs5_arsv_neutral,res_hhs6_arsv_neutral,
                              res_hhs7_arsv_neutral,res_hhs8_arsv_neutral,res_hhs9_arsv_neutral,res_hhs10_arsv_neutral)
  est_res_arsv_neutral<-get_est_all(res_arsv_neutral_list)
  get_fitmeasure_all(res_arsv_neutral_list)%>%
    full_join(est_res_arsv_neutral,
              by = c("HHSregion" =  "HHSregion", "Pathogen2"="pathogen2","Hypothesis" = "hyphothesis"))->result_res_arsv_neutral
  
  
  traj_fit_arsv_neutral<-get_traj_fitall(data=inc_data_add,res_arsv_neutral_list)
  
  rm(res_hhs1_arsv_neutral,res_hhs2_arsv_neutral,res_hhs3_arsv_neutral,
     res_hhs4_arsv_neutral,res_hhs5_arsv_neutral,res_hhs6_arsv_neutral,
     res_hhs7_arsv_neutral,res_hhs8_arsv_neutral,res_hhs9_arsv_neutral,res_hhs10_arsv_neutral,res_arsv_neutral_list)
  
  
  #########################################################################################   
  
  res_arsv_psi<-list.files(path="pomp_longphi_result/res_arsv_psi/",pattern=".rds")
  for (i in 1:length(res_arsv_psi)){
    load(paste("pomp_longphi_result/res_arsv_psi/",res_arsv_psi[i],sep=""))
    print(paste0("load ",res_arsv_psi[i]))
  }
  
  res_arsv_psi_list<-list(res_hhs1_arsv_psi,res_hhs2_arsv_psi,res_hhs3_arsv_psi,
                          res_hhs4_arsv_psi,res_hhs5_arsv_psi,res_hhs6_arsv_psi,
                          res_hhs7_arsv_psi,res_hhs8_arsv_psi,res_hhs9_arsv_psi,res_hhs10_arsv_psi)
  est_res_arsv_psi<-get_est_all(res_arsv_psi_list)
  get_fitmeasure_all(res_arsv_psi_list)%>%
    full_join(est_res_arsv_psi,
              by = c("HHSregion" =  "HHSregion", "Pathogen2"="pathogen2","Hypothesis" = "hyphothesis"))->result_res_arsv_psi
  
  traj_fit_arsv_psi<-get_traj_fitall(data=inc_data_add,res_arsv_psi_list)
  
  rm(res_hhs1_arsv_psi,res_hhs2_arsv_psi,res_hhs3_arsv_psi,
     res_hhs4_arsv_psi,res_hhs5_arsv_psi,res_hhs6_arsv_psi,
     res_hhs7_arsv_psi,res_hhs8_arsv_psi,res_hhs9_arsv_psi,res_hhs10_arsv_psi,res_arsv_psi_list)
  
  #########################################################################################   
  
  
  res_arsv_chi<-list.files(path="pomp_longphi_result/res_arsv_chi/",pattern=".rds")
  for (i in 1:length(res_arsv_chi)){
    load(paste("pomp_longphi_result/res_arsv_chi/",res_arsv_chi[i],sep=""))
    print(paste0("load ",res_arsv_chi[i]))
    
  }
  
  res_hhs2_arsv_chi$Hypothesis="chi"
  res_arsv_chi_list<-list(res_hhs1_arsv_chi,res_hhs2_arsv_chi,res_hhs3_arsv_chi,
                          res_hhs4_arsv_chi,res_hhs5_arsv_chi,res_hhs6_arsv_chi,
                          res_hhs7_arsv_chi,res_hhs8_arsv_chi,res_hhs9_arsv_chi,res_hhs10_arsv_chi)
  est_res_arsv_chi<-get_est_all(res_arsv_chi_list)
  get_fitmeasure_all(res_arsv_chi_list)%>%
    full_join(est_res_arsv_chi,
              by = c("HHSregion" =  "HHSregion", "Pathogen2"="pathogen2","Hypothesis" = "hyphothesis"))->result_res_arsv_chi
  traj_fit_arsv_chi<-get_traj_fitall(data=inc_data_add,res_arsv_chi_list)
  
  rm(res_hhs1_arsv_chi,res_hhs2_arsv_chi,res_hhs3_arsv_chi,
     res_hhs4_arsv_chi,res_hhs5_arsv_chi,res_hhs6_arsv_chi,
     res_hhs7_arsv_chi,res_hhs8_arsv_chi,res_hhs9_arsv_chi,res_hhs10_arsv_chi,res_arsv_chi_list)
  #########################################################################################   
  #est_res_arsv<-rbind(est_res_arsv_neutral,est_res_arsv_psi,est_res_arsv_chi,est_res_arsv_coinfect)
  
  rbind(result_res_arsv_neutral,result_res_arsv_psi,result_res_arsv_chi,result_res_arsv_coinfect)%>%
    select("Pathogen2","HHSregion","Hypothesis",
           "R01","R02","amplitude1","amplitude2","tpeak1","tpeak2","rho1","rho2", "psi","chi",
           "loglik","AIC","R2","RMSE" )%>%
    mutate(psi=1-psi,
           chi=1-chi)%>%
    melt(id.vars=c("Pathogen2","HHSregion","Hypothesis"))%>%
    dcast(Pathogen2 + HHSregion+variable~ Hypothesis)%>%
    select(Pathogen2, HHSregion, variable, neutral,psi,chi,co_infect) -> result_res_arsv
  
  write.csv(result_res_arsv,"./figures/result_res_arsv_longphhi_check.csv",row.names = FALSE)
  
 result_res_arsv<-read.csv("./figures/result_res_arsv_longphhi.csv")
 
 result_res_arsv%>%
   filter(variable=="AIC")%>%
   gather(hypothesis, AIC, neutral:co_infect, factor_key=TRUE)%>%
   group_by(HHSregion)%>%
   summarise(min_AIC=min(AIC))->arsv_minAIC
   
 result_res_arsv%>%
   filter(variable=="AIC")%>%
   gather(hypothesis, AIC, neutral:co_infect, factor_key=TRUE)%>%
   left_join(arsv_minAIC,by=c("HHSregion"="HHSregion"))%>%
   mutate(delta_AIC=AIC-min_AIC)%>%
   select(Pathogen2,HHSregion,hypothesis,delta_AIC)%>%
   spread(hypothesis,delta_AIC)%>%
   mutate(variable="deltaAIC")->arsv_deltaAIC
 
 result_res_arsv%>%
   filter(variable=="R2"| variable=="RMSE")%>%
   rbind(arsv_deltaAIC)->arsv_modelcompare

 write.csv(arsv_modelcompare,"./figures/arsv_modelcompare.csv",row.names = FALSE)
 
 
 
 
 
 
  
## plot compartment for sepcific HHS region
## HHS5: both psi and chi <1
test_traj(res_hhs = res_hhs3_arsv_neutral)%>%
  slice(., 2:n()) %>% 
  select(., -`.id`) %.>% 
  gather(., "comp", "count", -time) %>% 
  ggplot(., aes(x = time + 2011, y = count)) +
  geom_line()+
  facet_wrap(.~comp, scales = "free")+
  theme_classic()+
  xlab("Time")+
  ylab("Simulated cases under neutral model")+
  scale_x_continuous(breaks = seq(2014, 2018, by = 2))->HHS3_neutral
  
  
test_traj(res_hhs = res_hhs3_arsv_coinfect)%>%
  slice(., 2:n()) %>% 
  select(., -`.id`) %.>% 
  gather(., "comp", "count", -time) %>% 
  ggplot(., aes(x = time + 2011 , y = count)) +
  geom_line()+
  facet_wrap(.~comp, scales = "free")  +
  theme_classic()+
  xlab("Time")+
  ylab("Simulated cases under inhibition of co-infection and cross-immunity model")+
  scale_x_continuous(breaks = seq(2014, 2018, by = 2))->HHS3_coinfect

ggarrange(
  HHS3_neutral, HHS3_coinfect, labels = c("A", "B"),
  nrow=2)

  
  
  
  
## HHS 10: only Psi <1
  
  
  
  
  
  
    
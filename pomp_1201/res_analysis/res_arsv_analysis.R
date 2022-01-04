
#Sys.setenv('R_MAX_VSIZE'=64000000000)
#Sys.getenv('R_MAX_VSIZE')
# first load prerequisites
source("./fit_functions.R", chdir = TRUE) 
#########################################################################################
  res_arsv_coinfect<-list.files(path="pomp_result_1213/res_arsv_coinfect/",pattern=".rds")
  for (i in 1:length(res_arsv_coinfect)){
    load(paste("pomp_result_1213/res_arsv_coinfect/",res_arsv_coinfect[i],sep=""))
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

  
  traj_fit_arsv_coinfect<-get_traj_fitall(res_arsv_coinfect_list)
  
## clear large list to release memory  
  rm(res_arsv_coinfect_list,res_hhs1_arsv_coinfect,res_hhs2_arsv_coinfect,res_hhs3_arsv_coinfect,
     res_hhs4_arsv_coinfect,res_hhs5_arsv_coinfect,res_hhs6_arsv_coinfect,
     res_hhs7_arsv_coinfect,res_hhs8_arsv_coinfect,res_hhs9_arsv_coinfect,res_hhs10_arsv_coinfect)

  
  #########################################################################################  
  
  
  res_arsv_neutral<-list.files(path="pomp_result_1213/res_arsv_neutral/",pattern=".rds")
  for (i in 1:length(res_arsv_neutral)){
    load(paste("pomp_result_1213/res_arsv_neutral/",res_arsv_neutral[i],sep=""))
    print(paste0("load ",res_arsv_neutral[i]))
  }
  
  res_arsv_neutral_list<-list(res_hhs1_arsv_neutral,res_hhs2_arsv_neutral,res_hhs3_arsv_neutral,
                              res_hhs4_arsv_neutral,res_hhs5_arsv_neutral,res_hhs6_arsv_neutral,
                              res_hhs7_arsv_neutral,res_hhs8_arsv_neutral,res_hhs9_arsv_neutral,res_hhs10_arsv_neutral)
  est_res_arsv_neutral<-get_est_all(res_arsv_neutral_list)
  get_fitmeasure_all(res_arsv_neutral_list)%>%
    full_join(est_res_arsv_neutral,
              by = c("HHSregion" =  "HHSregion", "Pathogen2"="pathogen2","Hypothesis" = "hyphothesis"))->result_res_arsv_neutral
  
  
  traj_fit_arsv_neutral<-get_traj_fitall(res_arsv_neutral_list)
  
  rm(res_hhs1_arsv_neutral,res_hhs2_arsv_neutral,res_hhs3_arsv_neutral,
     res_hhs4_arsv_neutral,res_hhs5_arsv_neutral,res_hhs6_arsv_neutral,
     res_hhs7_arsv_neutral,res_hhs8_arsv_neutral,res_hhs9_arsv_neutral,res_hhs10_arsv_neutral,res_arsv_neutral_list)
  
  
  #########################################################################################   
  
  res_arsv_psi<-list.files(path="pomp_result_1213/res_arsv_psi/",pattern=".rds")
  for (i in 1:length(res_arsv_psi)){
    load(paste("pomp_result_1213/res_arsv_psi/",res_arsv_psi[i],sep=""))
    print(paste0("load ",res_arsv_psi[i]))
  }
  
  res_arsv_psi_list<-list(res_hhs1_arsv_psi,res_hhs2_arsv_psi,res_hhs3_arsv_psi,
                          res_hhs4_arsv_psi,res_hhs5_arsv_psi,res_hhs6_arsv_psi,
                          res_hhs7_arsv_psi,res_hhs8_arsv_psi,res_hhs9_arsv_psi,res_hhs10_arsv_psi)
  est_res_arsv_psi<-get_est_all(res_arsv_psi_list)
  get_fitmeasure_all(res_arsv_psi_list)%>%
    full_join(est_res_arsv_psi,
              by = c("HHSregion" =  "HHSregion", "Pathogen2"="pathogen2","Hypothesis" = "hyphothesis"))->result_res_arsv_psi
  
  traj_fit_arsv_psi<-get_traj_fitall(res_arsv_psi_list)
  
  rm(res_hhs1_arsv_psi,res_hhs2_arsv_psi,res_hhs3_arsv_psi,
     res_hhs4_arsv_psi,res_hhs5_arsv_psi,res_hhs6_arsv_psi,
     res_hhs7_arsv_psi,res_hhs8_arsv_psi,res_hhs9_arsv_psi,res_hhs10_arsv_psi,res_arsv_psi_list)
  
  #########################################################################################   
  
  
  res_arsv_chi<-list.files(path="pomp_result_1213/res_arsv_chi/",pattern=".rds")
  for (i in 1:length(res_arsv_chi)){
    load(paste("pomp_result_1213/res_arsv_chi/",res_arsv_chi[i],sep=""))
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
  traj_fit_arsv_chi<-get_traj_fitall(res_arsv_chi_list)
  
  rm(res_hhs1_arsv_chi,res_hhs2_arsv_chi,res_hhs3_arsv_chi,
     res_hhs4_arsv_chi,res_hhs5_arsv_chi,res_hhs6_arsv_chi,
     res_hhs7_arsv_chi,res_hhs8_arsv_chi,res_hhs9_arsv_chi,res_hhs10_arsv_chi,res_arsv_chi_list)
  #########################################################################################   
  est_res_arsv<-rbind(est_res_arsv_neutral,est_res_arsv_psi,est_res_arsv_chi,est_res_arsv_coinfect)
  
  rbind(result_res_arsv_neutral,result_res_arsv_psi,result_res_arsv_chi,result_res_arsv_coinfect)%>%
    select("Pathogen2","HHSregion","Hypothesis",
           "R01","R02","amplitude1","amplitude2","tpeak1","tpeak2","rho1","rho2", "psi","chi",
           "loglik","AIC","R2","RMSE" )-> result_res_arsv
  
  write.csv(result_res_arsv,"./figures/result_res_arsv.csv",row.names = FALSE)
  
  hypo_levels=c("neutral","psi","chi","co-infect")
 
  est_res_arsv%>%
    mutate(hyphothesis=factor(hyphothesis,levels=hypo_levels))%>%
    ggplot(aes(x=hyphothesis,y=AIC,color=hyphothesis))+
    geom_point(size=3)+
    facet_wrap(~HHSregion,scales="free_y",nrow = 1)+
    theme_bw()+
    theme(legend.position="bottom")+
    #scale_color_discrete(name = "Hypothesis", labels = c("neutral", "inhibition on co-infection", "cross-protection","co-ifection + cross-protection"))+
    scale_colour_brewer(palette = "Dark2",name = "Hypothesis", 
                        labels = c("neutral", "inhibition on co-infection", "cross-protection","co-ifection + cross-protection"))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())->arsv_AIC
  
  
  
## plot trajectory fit
  inc_data_fit<-inc_data_perdic%>%
    mutate(HHSregion=HHS_REGION)%>%
    pivot_wider(names_from = virus, values_from = cases)%>%
    drop_na()%>%
    mutate(method=case_when(
      time<=6.5 ~ "Model fit",
      time>6.5 ~ "Project"
    ))
  
  traj_fit_arsv<-rbind(traj_fit_arsv_neutral,traj_fit_arsv_psi,traj_fit_arsv_chi,traj_fit_arsv_coinfect)
  traj_fit_arsv_coinfect%>%
    mutate(hypothesis="co_infect")
  hypo_levels=c("neutral","psi","chi","co_infect")
  patho_order<-c("RSV","fluA")
  traj_fit_arsv%>%
    mutate(hypothesis=factor(hypothesis,levels=hypo_levels))%>%
    mutate(pathogen=factor(type,levels=patho_order))%>%
    mutate(method=case_when(
      time<=6.5 ~ "Model fit",
      time>6.5 ~ "Project"
    ))%>%
    ggplot(aes(x=time+2011,y=cases))+
    geom_area(data=inc_data_fit,aes(x=date,y=RSV),fill="brown",alpha=0.4)+

    geom_area(data=inc_data_fit,aes(x=date,y=fluA),fill="gray",alpha=0.6)+
    geom_line(aes(color=hypothesis,linetype=pathogen),size=0.5)+
    facet_grid(HHSregion~method,scales="free_x",space="free")+


    geom_area(data=inc_data_fit,aes(x=date,y=fluA),fill="gray",alpha=0.4)+
    geom_line(aes(color=hypothesis,linetype=type))+
    facet_grid(HHSregion~method,scales= "free_x", space="free_x")+
    scale_colour_brewer(palette = "Dark2")+
    theme_bw()+
    xlab("time(weeks)")+
    scale_colour_brewer(palette = "Dark2",name = "Hypothesis", 
                        labels = c("neutral", "inhibition on co-infection", "cross-protection","co-infection + cross-protection")) ->arsv_trajfit
#+
 #   theme(legend.position="bottom"
    

  
## plot compartment for sepcific HHS region
## HHS5: both psi and chi <1
test_traj(res_hhs = res_hhs5_arsv_coinfect)%>%
  slice(., 2:n()) %>% 
  select(., -`.id`) %.>% 
  gather(., "comp", "count", -time) %>% 
  ggplot(., aes(x = time + 2011, y = count)) +
  geom_line()+
  facet_wrap(.~comp, scales = "free")+
  xlab("Time")+
  ylab("Simulated cases under neutral model")
  
  
test_traj(res_hhs = res_hhs5_arsv_neutral)%>%
  slice(., 2:n()) %>% 
  select(., -`.id`) %.>% 
  gather(., "comp", "count", -time) %>% 
  ggplot(., aes(x = time + 2011 , y = count)) +
  geom_line()+
  facet_wrap(.~comp, scales = "free")  +
  xlab("Time")+
  ylab("Simulated cases under co-infection and cross-protection model")

  
  
  
  
## HHS 10: only Psi <1
  
  
  
  
  
  
    
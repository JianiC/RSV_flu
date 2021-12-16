
#Sys.setenv('R_MAX_VSIZE'=64000000000)
#Sys.getenv('R_MAX_VSIZE')
source("./setup.R", chdir = TRUE)
source("./util_funs.R", chdir = TRUE)
source("./pomp_1204.R", chdir = TRUE)
source("./fit_functions.R")
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
  
  
  rm(res_arsv_coinfect_list,res_hhs1_arsv_coinfect,res_hhs2_arsv_coinfect,res_hhs3_arsv_coinfect,
     res_hhs4_arsv_coinfect,res_hhs5_arsv_coinfect,res_hhs6_arsv_coinfect,
     res_hhs7_arsv_coinfect,res_hhs8_arsv_coinfect,res_hhs9_arsv_coinfect,res_hhs10_arsv_coinfect)
  ######
  
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
  
  rm(res_hhs1_arsv_psi,res_hhs2_arsv_psi,res_hhs3_arsv_psi,
     res_hhs4_arsv_psi,res_hhs5_arsv_psi,res_hhs6_arsv_psi,
     res_hhs7_arsv_psi,res_hhs8_arsv_psi,res_hhs9_arsv_psi,res_hhs10_arsv_psi,res_arsv_psi_list)
  
  #########################################################################################   
  
  
  res_arsv_chi<-list.files(path="pomp_result_1213/res_arsv_chi/",pattern=".rds")
  for (i in 1:length(res_arsv_chi)){
    load(paste("pomp_result_1213/res_arsv_chi/",res_arsv_chi[i],sep=""))
    print(paste0("load ",res_arsv_chi[i]))
    
  }
  
  
  res_arsv_chi_list<-list(res_hhs1_arsv_chi,res_hhs2_arsv_chi,res_hhs3_arsv_chi,
                          res_hhs4_arsv_chi,res_hhs5_arsv_chi,res_hhs6_arsv_chi,
                          res_hhs7_arsv_chi,res_hhs8_arsv_chi,res_hhs9_arsv_chi,res_hhs10_arsv_chi)
  est_res_arsv_chi<-get_est_all(res_arsv_chi_list)
  rm(res_hhs1_arsv_chi,res_hhs2_arsv_chi,res_hhs3_arsv_chi,
     res_hhs4_arsv_chi,res_hhs5_arsv_chi,res_hhs6_arsv_chi,
     res_hhs7_arsv_chi,res_hhs8_arsv_chi,res_hhs9_arsv_chi,res_hhs10_arsv_chi,res_arsv_chi_list)
  #########################################################################################   
  est_res_arsv<-rbind(est_res_arsv_neutral,est_res_arsv_psi,est_res_arsv_chi,est_res_arsv_coinfect)
  hypo_levels=c("neutral","psi","chi","co-infect")
 
  est_res_arsv%>%
    mutate(hyphothesis=factor(hyphothesis,levels=hypo_levels))%>%
    ggplot(aes(x=hyphothesis,y=loglik,color=hyphothesis))+
    geom_point(size=3)+
    facet_wrap(~HHSregion,scales="free_y",nrow = 1)+
    theme_bw()+
    theme(legend.position="bottom")+
    scale_colour_brewer(palette = "Dark2")
    
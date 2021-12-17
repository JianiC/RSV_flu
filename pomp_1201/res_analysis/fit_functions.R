## function to get estimated param values
source("./setup.R", chdir = TRUE)
source("./util_funs.R", chdir = TRUE)



get_rp_vals<-function(data=incdata,res_hhs){
  
  HHS_region=res_hhs$HHS
  if(res_hhs$total2=="fluA"){
    virus_combo=c("RSV","fluA")
  }else{
    virus_combo=c("RSV","fluB")
  }
  ## pompdata
  pomp_data_hhs <- (
    inc_data %.>% 
      make_data_pomp_ready(., virus_combo =virus_combo, HHS_region = HHS_region)
  )
  
  ##
  if (res_hhs$Hypothesis == "neutral"){
    param_vector<-as.data.frame(res_hhs$DEobj$optim$bestmem)%>%
      cbind(row.names(.))%>%
      pivot_wider(names_from = `row.names(.)`, values_from = `res_hhs$DEobj$optim$bestmem`)%>%
      mutate(phi1=365/30, phi2=365/30,
             psi=1, chi =1)
  }else if (res_hhs$Hypothesis == "psi"){
    param_vector<-as.data.frame(res_hhs$DEobj$optim$bestmem)%>%
    cbind(row.names(.))%>%
    pivot_wider(names_from = `row.names(.)`, values_from = `res_hhs$DEobj$optim$bestmem`)%>%
    mutate(phi1=365/30, phi2=365/30,
           chi =1)
    }else if (res_hhs$Hypothesis == "chi") {
      param_vector<-as.data.frame(res_hhs$DEobj$optim$bestmem)%>%
        cbind(row.names(.))%>%
        pivot_wider(names_from = `row.names(.)`, values_from = `res_hhs$DEobj$optim$bestmem`)%>%
        mutate(phi1=365/30, phi2=365/30,
               psi =1)
    }else{
    param_vector<-as.data.frame(res_hhs$DEobj$optim$bestmem)%>%
      cbind(row.names(.))%>%
      pivot_wider(names_from = `row.names(.)`, values_from = `res_hhs$DEobj$optim$bestmem`)%>%
      mutate(phi1=365/30, phi2=365/30)
    
  }
  rp_vals <- c(R01=param_vector$R01, gamma1=365./9, w1=1/1,
               R02=param_vector$R02, gamma2=365./3, w2=1/1,
               amplitude1 = param_vector$amplitude1,
               amplitude2 = param_vector$amplitude2,
               tpeak1 = param_vector$tpeak1,
               tpeak2 = param_vector$tpeak2,
               phi1=param_vector$phi1, phi2=param_vector$phi2,
               psi=param_vector$psi, chi=param_vector$chi,
               eta1=365., eta2=365.,
               rho1=param_vector$rho1, rho2=param_vector$rho2,
               sigmaSE=0.0000,
               pop= pomp_data_hhs$N[1],
               mu=1/80)
  
  return(rp_vals)
  
  
}



##get testtraject using param function

test_traj<-function(data=incdata, virus_combo,HHS_region,res_hhs){
  
  params=get_rp_vals(data = incdata,res_hhs)
  
  ##pomp model
  pomp_data_hhs <- (
    inc_data %.>% 
      make_data_pomp_ready(., virus_combo =virus_combo, HHS_region = HHS_region)
  )
  pomp_data_hhs<-pomp_data_hhs%>%drop_na()
  hhs_po <- make_pomp(df = pomp_data_hhs, time_start_sim = -100)
  
  ## trajectory
  test_traj <- trajectory(object = hhs_po, params=params, format = "d", method = "ode23")
  
  return(test_traj)
  
}





test_traj_pseudo<-function(data=incdata, virus_combo,HHS_region,res_hhs, n=10){
  
  params=get_rp_vals(data = incdata,virus_combo, HHS_region ,res_hhs)
  
  ##pomp model
  pomp_data_hhs <- (
    inc_data %.>% 
      make_data_pomp_ready(., virus_combo =virus_combo, HHS_region = HHS_region)
  )
  
   pseudo_data <- tibble(time = seq(2.5, n, by = 1/52), 
                         total1 = NA, 
                         total2 = NA, N = pomp_data_hhs$N[1])
  
  # pseudo_data<-tibble(time=pomp_data_hhs$time,
  #                    total1=NA,
  #                    total2 = NA, N = pomp_data_hhs$N[1])
  # 
  hhs_po <- make_pomp(df = pseudo_data, time_start_sim = -100)
  
  ## trajectory
  test_traj <- trajectory(object = hhs_po, params=params, format = "d", method = "ode23")
  
  return(test_traj)
  
}


##get testtraject using param function

traj_fit<-function(data=inc_data_perdic,res_hhs){
  HHS_region=res_hhs$HHS
  if(res_hhs$total2=="fluA"){
    virus_combo=c("RSV","fluA")
  }else{
    virus_combo=c("RSV","fluB")
  }
  
  params=get_rp_vals(data=inc_data_perdic,virus_combo, HHS_region ,res_hhs)
  
  ##pomp model
  pomp_data_hhs <- (
    data %.>% 
      make_data_pomp_ready(., virus_combo =virus_combo, HHS_region = HHS_region)
      
  )
  
  hhs_po <- make_pomp(df = pomp_data_hhs, time_start_sim = -100)
  
  ## trajectory
  test_traj <- trajectory(object = hhs_po, params=params, format = "d", method = "ode23")
  test_traj %>% 
    slice(., 2:n()) %>% 
    #select(., -`.id`) %.>% 
    select(time, K1, K2) %>%
    mutate(scases1 = res_hhs$DEobj$optim$bestmem["rho1"]*K1, 
           scases2 = res_hhs$DEobj$optim$bestmem["rho2"]*K2)%>%
    select(time, scases1, scases2) %>%
    gather(key = "type", value = "cases", -time)%>%
    mutate(hypothesis=res_hhs$Hypothesi)%>%
    mutate(HHSregion=res_hhs$HHS,
           pathogen2=res_hhs$total2)->traj_fit 

   # pomp_data_hhs%>%
   #   mutate(scases1=total1,
   #          scases2=total2)%>%
   #   select(time, scases1, scases2) %>%
   #   gather(key = "type", value = "obscases", -time) %>%
   #   full_join(test_traj_data, by =c ("time" = "time","type"="type"))->fit_withdata
  
  
  return(traj_fit)
  
}

calculate_aic<-function(loglik,npar){
  aic=-2*(loglik)+2*npar
  return(aic)
}

getparam<-function(res_hhs){
  if (res_hhs$Hypothesis == "neutral"){
    param_vector<-as.data.frame(res_hhs$DEobj$optim$bestmem)%>%
      cbind(row.names(.))%>%
      pivot_wider(names_from = `row.names(.)`, values_from = `res_hhs$DEobj$optim$bestmem`)%>%
      mutate(phi1=365/30, phi2=365/30,
             psi=1, chi =1)%>%
      mutate(loglik = -(res_hhs$DEobj$optim$bestval))%>%
      mutate(AIC=calculate_aic(loglik,npar=8),
             hyphothesis="neutral",
             pathogen1=res_hhs$total1,
             pathogen2=res_hhs$total2,
             HHSregion=res_hhs$HHS)
    
  }else if (res_hhs$Hypothesis == "psi"){
    param_vector<-as.data.frame(res_hhs$DEobj$optim$bestmem)%>%
      cbind(row.names(.))%>%
      pivot_wider(names_from = `row.names(.)`, values_from = `res_hhs$DEobj$optim$bestmem`)%>%
      mutate(phi1=365/30, phi2=365/30,
             chi=1)%>%
      mutate(loglik = -(res_hhs$DEobj$optim$bestval))%>%
      mutate(AIC=calculate_aic(loglik,npar=8),
             hyphothesis="psi",
             pathogen1=res_hhs$total1,
             pathogen2=res_hhs$total2,
             HHSregion=res_hhs$HHS)
    
  }else if (res_hhs$Hypothesis == "chi") {
    param_vector<-as.data.frame(res_hhs$DEobj$optim$bestmem)%>%
      cbind(row.names(.))%>%
      pivot_wider(names_from = `row.names(.)`, values_from = `res_hhs$DEobj$optim$bestmem`)%>%
      mutate(phi1=365/30, phi2=365/30,
             psi=1)%>%
      mutate(loglik = -(res_hhs$DEobj$optim$bestval))%>%
      mutate(AIC=calculate_aic(loglik,npar=8),
             hyphothesis="chi",
             pathogen1=res_hhs$total1,
             pathogen2=res_hhs$total2,
             HHSregion=res_hhs$HHS,
             )
    
  }else{
    param_vector<-as.data.frame(res_hhs$DEobj$optim$bestmem)%>%
      cbind(row.names(.))%>%
      pivot_wider(names_from = `row.names(.)`, values_from = `res_hhs$DEobj$optim$bestmem`)%>%
      mutate(phi1=365/30, phi2=365/30)%>%
      mutate(loglik = -(res_hhs$DEobj$optim$bestval))%>%
      mutate(AIC=calculate_aic(loglik,npar=10),
             hyphothesis="co-infect",
             pathogen1=res_hhs$total1,
             pathogen2=res_hhs$total2,
             HHSregion=res_hhs$HHS)
    
  }
  
  return(param_vector)
  
}

get_est_all<-function(res_list){
  df_est<-data.frame()
  for( i in 1:length(res_list)){
    out<-getparam(res_list[[i]])
    df_est=rbind(df_est,out)
    
  }
  return(df_est)
}





get_traj_fitall<-function(res_list){
  
  df_traj_fit<-data.frame()
  
  for( i in 1:length(res_list)){
    out<-traj_fit(res_hhs=res_list[[i]])
    df_traj_fit<-rbind(df_traj_fit,out)
    
  }
  return(df_traj_fit)
}




loglik<-function(fit_par,pomp_data){
  
  model_variables = c("S1_S2","I1_S2", "C1_S2","R1_S2",
                      "S1_I2","I1_I2", "C1_I2","R1_I2",
                      "S1_C2","I1_C2", "C1_C2","R1_C2",
                      "S1_R2","I1_R2", "C1_R2","R1_R2")
  model_params = c("R01","R02", "gamma1", "gamma2",
                   "w1","w2","eta1", "eta2", "phi1", "phi2",
                   "psi","chi", "rho1", "rho2",
                   "amplitude1", "tpeak1", "amplitude2", "tpeak2", "pop","mu")
  
  pomp_objfun <- (
    make_pomp(pomp_data,time_start_sim = -100) %.>%
      # define the objective function 
      traj_objfun(., 
                  est = names(fit_par), 
                  params = fit_par ,
                  dmeasure=dmeas_poisson,fail.value = 1e20,
                  statenames = c(model_variables, c("Kprim1","Ksec1","Kprim2","Ksec2","K1", "K2")),
                  ode_control=list(method="ode23"),
                  paramnames = c(model_params))
  )
  
  return(pomp_objfun(fit_par))
  
}


par_loglik<-function(res_hhs,par){
  ## get pomp_data and estimated parm for loglik
  
  HHS_region=res_hhs$HHS
  if(res_hhs$total2=="fluA"){
    virus_combo=c("RSV","fluA")
  }else{
    virus_combo=c("RSV","fluB")
  }
  
  
  ##pomp model
  pomp_data_hhs <- (
    inc_data_add %.>% 
      make_data_pomp_ready(., virus_combo =virus_combo, HHS_region = HHS_region))
  
  fit_par <- get_rp_vals(data = inc_data_add,res_hhs)
  result<-data.frame()
  ## change parms to get fit par
  
  if (par=="psi"){
    for (i in seq(0, 1, by = 0.01)) {
      
      fit_par["psi"]<-i
      loglik_value <-loglik(fit_par = fit_par,pomp_data = pomp_data_hhs)
      
      out<-list("psi"=i,"loglik"=-loglik_value)
      #print(out)
      result<-rbind(result,out)
    }
    print("loglik estimation for psi is done")
  }else{
    for (i in seq(0, 1, by = 0.01)) {
      
      fit_par["chi"]<-i
      loglik_value <-loglik(fit_par = fit_par,pomp_data = pomp_data_hhs)
      
      out<-list("chi"=i,"loglik"=-loglik_value)
      #print(out)
      result<-rbind(result,out)
      
    }
    print("loglik estimation for chi is done")
  }
  
  
  
  result%<>%
    mutate(HHS_region=res_hhs$HHS,
           Pathens2=res_hhs$total2)
  
  return(result)
  
}

## from parmeter loglike profile dataframe to get confidence interval
CI_loglik_psi<-function(loglik_p){
  
  maxloglik <- max(loglik_p$loglik)
  cutoff <- maxloglik-qchisq(p=0.95,df=1)/2
  lower <-range(subset(loglik_p,loglik>cutoff)$psi)[1]
  upper <-range(subset(loglik_p,loglik>cutoff)$psi)[2]
  
  loglik_p%>% mutate(loglik_cutoff=cutoff,
                     maxloglik=maxloglik,
                     lower=lower,
                     upper=upper)->CI_parm
  
  return(CI_parm)
  
  
}

CI_loglik_chi<-function(loglik_p){
  
  maxloglik <- max(loglik_p$loglik)
  cutoff <- maxloglik-qchisq(p=0.95,df=1)/2
  lower <-range(subset(loglik_p,loglik>cutoff)$chi)[1]
  upper <-range(subset(loglik_p,loglik>cutoff)$chi)[2]
  
  loglik_p%>% mutate(loglik_cutoff=cutoff,
                     maxloglik=maxloglik,
                     lower=lower,
                     upper=upper)->CI_parm
  
  return(CI_parm)
  
  
}


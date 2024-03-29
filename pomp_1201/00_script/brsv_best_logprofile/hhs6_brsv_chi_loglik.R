# this script runs different functions to make sure that they are running alright
## estimate the loglik profile with different psi and chi value
# first load prerequisites
source("./src.R", chdir = TRUE)

## specify the virus and HHS region

# make data ready for pomp
pomp_data_hhs6_brsv <- (
inc_data_add %.>%
  make_data_pomp_ready(., virus_combo = c("RSV", "fluB"), HHS_region = 6)
  )

## from estimation
load("./res_hhs6_brsv_chi.rds")

  best_past_est<-c(R01=1.279869,
                          R02=1.14500,
                          amplitude1=0.258257,
                          amplitude2 =0.110144,
                          tpeak1=0.805541,
                          tpeak2=0.06007,
                          rho1=0.000703,
                          rho2=0.000324)
                          
chi_vals<-seq(0,1,by=0.1)

## loop different chi values
res_hhs6_brsv_chi<-data.frame()

for (i in 1:length(chi_vals)){
  
  rp_vals_def <- c(R01 = 1, gamma1=365./9, w1=2,
                   R02 = 1, gamma2=365./3, w2=2,
                   phi1=365/180, phi2=365/180, psi =1.0, chi=chi_vals[i],
                   eta1=365., eta2=365.,rho1 = 0, rho2 = 0,
                   amplitude1=0.0, tpeak1=0.0, amplitude2=0.0, tpeak2=0.0,
                   pop=pomp_data_hhs6_brsv$N[1],
                   mu=1/80)
  
  res_hhs6_brsv_chi_tmp <- (
    DE_traj_match(df = pomp_data_hhs6_brsv, 
                  param_constraints = neutral_param_constraints, 
                  params = rp_vals_def,
                  ode_control = list(method = "ode23"), 
                  best_past_est = best_past_est,
                  hypo_name = "chi_loglik", 
                  hhs_reg = 6, 
                  tot1_name = "RSV", 
                  tot2_name = "fluB")
  )
  
  #if(res_hhs6_brsv_coinfect$total2 == "fluB") message("Code itegration complete!")
  print(chi_vals[i])
  print(res_hhs6_brsv_chi_tmp$DEobj$optim$bestval)
  
  res_hhs6_brsv_chi_tmp <- data.frame (HHS_region =6,
                                       pathogen2 = "fluB",
                                       chi=chi_vals[i],
                                       loglik=-res_hhs6_brsv_chi_tmp$DEobj$optim$bestva)
  res_hhs6_brsv_chi <-rbind(res_hhs6_brsv_chi,res_hhs6_brsv_chi_tmp)
  
}

write.csv(res_hhs6_brsv_chi,"res_hhs6_brsv_chi_loglik_s2.csv",row.names=FALSE)


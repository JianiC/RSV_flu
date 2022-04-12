# this script runs different functions to make sure that they are running alright
## estimate the loglik profile with different psi and chi value
# first load prerequisites
source("./src.R", chdir = TRUE)

## specify the virus and HHS region

# make data ready for pomp
pomp_data_hhs2_brsv <- (
inc_data_add %.>%
  make_data_pomp_ready(., virus_combo = c("RSV", "fluB"), HHS_region = 2)
  )

## from estimation
load("res_hhs2_brsv_chi.rds")

best_past_est<-c(R01=past_est(res_hhs2_brsv_chi)$R01,
                          R02=past_est(res_hhs2_brsv_chi)$R02,
                          amplitude1=past_est(res_hhs2_brsv_chi)$amplitude1,
                          amplitude2 =past_est(res_hhs2_brsv_chi)$amplitude2,
                          tpeak1=past_est(res_hhs2_brsv_chi)$tpeak1,
                          tpeak2=past_est(res_hhs2_brsv_chi)$tpeak2,
                          rho1=past_est(res_hhs2_brsv_chi)$rho1,
                          rho2=past_est(res_hhs2_brsv_chi)$rho2)




chi_mle_low <- 0.75*chi_mle
chi_mle_up <- 1.25*chi_mle
chi_vals<-c(seq(chi_mle_low,chi_mle_up,by=0.03),chi_mle)

## loop different chi values
res_hhs2_brsv_chi<-data.frame()

for (i in 1:length(chi_vals)){
  
  rp_vals_def <- c(R01 = 1, gamma1=365./9, w1=2,
                   R02 = 1, gamma2=365./3, w2=2,
                   phi1=365/180, phi2=365/180, psi =1.0, chi=chi_vals[i],
                   eta1=365., eta2=365.,rho1 = 0, rho2 = 0,
                   amplitude1=0.0, tpeak1=0.0, amplitude2=0.0, tpeak2=0.0,
                   pop=pomp_data_hhs2_brsv$N[1],
                   mu=1/80)
  
  res_hhs2_brsv_chi_tmp <- (
    DE_traj_match(df = pomp_data_hhs2_brsv, 
                  param_constraints = neutral_param_constraints, 
                  params = rp_vals_def,
                  ode_control = list(method = "ode23"),
                  best_past_est = best_past_est, 
                  hypo_name = "chi_loglik", 
                  hhs_reg = 2, 
                  tot1_name = "RSV", 
                  tot2_name = "fluB")
  )
  
  if(res_hhs2_brsv_chi_tmp$total2 == "fluB") message("Code itegration complete!")
  print(chi_vals[i])
  print(res_hhs2_brsv_chi_tmp$DEobj$optim$bestval)
  
  res_hhs2_brsv_chi_tmp <- data.frame (HHS_region =2,
                                       pathogen2 = "fluB",
                                       chi=chi_vals[i],
                                       loglik=-res_hhs2_brsv_chi_tmp$DEobj$optim$bestva)
  res_hhs2_brsv_chi <-rbind(res_hhs2_brsv_chi,res_hhs2_brsv_chi_tmp)
  
}

write.csv(res_hhs2_brsv_chi,"res_hhs2_brsv_chi_loglik_s2.csv",row.names=FALSE)


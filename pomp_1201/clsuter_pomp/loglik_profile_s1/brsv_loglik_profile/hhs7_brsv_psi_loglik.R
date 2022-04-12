# this script runs different functions to make sure that they are running alright
## estimate the loglik profile with different psi and psi value
# first load prerequisites
source("./src.R", chdir = TRUE)

## specify the virus and HHS region

# make data ready for pomp
pomp_data_hhs7_brsv <- (
inc_data_add %.>%
  make_data_pomp_ready(., virus_combo = c("RSV", "fluB"), HHS_region = 7)
  )

## from estimation
load("./res_hhs7_brsv_coinfect.rds")

best_past_est<-data.frame(past_est(res_hhs7_brsv_coinfect)$R01,
                          past_est(res_hhs7_brsv_coinfect)$R02,
                          past_est(res_hhs7_brsv_coinfect)$amplitude1,
                          past_est(res_hhs7_brsv_coinfect)$amplitude2,
                          past_est(res_hhs7_brsv_coinfect)$tpeak1,
                          past_est(res_hhs7_brsv_coinfect)$tpeak2,
                          past_est(res_hhs7_brsv_coinfect)$rho1,
                          past_est(res_hhs7_brsv_coinfect)$rho2)

psi_vals<-c(seq(0,1,by=0.1),past_est(res_hhs7_brsv_coinfect)$psi)

## loop different psi values
res_hhs7_brsv_psi<-data.frame()

for (i in 1:length(psi_vals)){
  
  rp_vals_def <- c(R01 = 1, gamma1=365./9, w1=2,
                   R02 = 1, gamma2=365./3, w2=2,
                   phi1=365/180, phi2=365/180, psi =psi_vals[i], chi=1,
                   eta1=365., eta2=365.,rho1 = 0, rho2 = 0,
                   amplitude1=0.0, tpeak1=0.0, amplitude2=0.0, tpeak2=0.0,
                   pop=pomp_data_hhs7_brsv$N[1],
                   mu=1/80)
  
  res_hhs7_brsv_psi_tmp <- (
    DE_traj_match(df = pomp_data_hhs7_brsv, 
                  param_constraints = chi_param_constraints, 
                  params = rp_vals_def,
                  ode_control = list(method = "ode23"), 
                  hypo_name = "psi_loglik", 
                  hhs_reg = 7, 
                  tot1_name = "RSV", 
                  tot2_name = "fluB")
  )
  
  if(res_hhs7_brsv_coinfect$total2 == "fluB") message("Code itegration complete!")
  print(psi_vals[i])
  print(res_hhs7_brsv_psi_tmp$DEobj$optim$bestval)
  
  res_hhs7_brsv_psi_tmp <- data.frame (HHS_region =7,
                                       pathogen2 = "fluB",
                                       psi=psi_vals[i],
                                       loglik=-res_hhs7_brsv_psi_tmp$DEobj$optim$bestva)
  res_hhs7_brsv_psi <-rbind(res_hhs7_brsv_psi,res_hhs7_brsv_psi_tmp)
  
}

write.csv(res_hhs7_brsv_psi,"res_hhs7_brsv_psi_loglik.csv",row.names=FALSE)


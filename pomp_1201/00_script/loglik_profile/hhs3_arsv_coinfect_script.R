# this script runs different functions to make sure that they are running alright
## estimate the loglik profile with different psi and chi value
# first load prerequisites
source("./src.R", chdir = TRUE)

## specify the virus and HHS region

# make data ready for pomp
pomp_data_hhs3_arsv <- (
inc_data_add %.>%
  make_data_pomp_ready(., virus_combo = c("RSV", "fluA"), HHS_region = 3)
  )

## from estimation
load("./res_hhs3_arsv_coinfect.rds")

best_past_est<-data.frame(past_est(res_hhs3_arsv_coinfect)$R01,
                          past_est(res_hhs3_arsv_coinfect)$R02,
                          past_est(res_hhs3_arsv_coinfect)$amplitude1,
                          past_est(res_hhs3_arsv_coinfect)$amplitude2,
                          past_est(res_hhs3_arsv_coinfect)$tpeak1,
                          past_est(res_hhs3_arsv_coinfect)$tpeak2,
                          past_est(res_hhs3_arsv_coinfect)$rho1,
                          past_est(res_hhs3_arsv_coinfect)$rho2)

chi_vals<-c(seq(0,1,by=0.1),past_est(res_hhs3_arsv_coinfect)$chi)

## loop different chi values
res_hhs3_arsv_chi<-data.frame()

for (i in 1:length(chi_vals)){
  
  rp_vals_def <- c(R01 = 1, gamma1=365./9, w1=2,
                   R02 = 1, gamma2=365./3, w2=2,
                   phi1=365/180, phi2=365/180, psi =1.0, chi=chi_vals[i],
                   eta1=365., eta2=365.,rho1 = 0, rho2 = 0,
                   amplitude1=0.0, tpeak1=0.0, amplitude2=0.0, tpeak2=0.0,
                   pop=pomp_data_hhs3_arsv$N[1],
                   mu=1/80)
  
  res_hhs3_arsv_chi_tmp <- (
    DE_traj_match(df = pomp_data_hhs3_arsv, 
                  param_constraints = psi_param_constraints, 
                  params = rp_vals_def,
                  ode_control = list(method = "ode23"), 
                  hypo_name = "chi_loglik", 
                  hhs_reg = 3, 
                  tot1_name = "RSV", 
                  tot2_name = "fluA")
  )
  
  if(res_hhs3_arsv_coinfect$total2 == "fluA") message("Code itegration complete!")
  print(chi_vals[i])
  print(res_hhs3_arsv_chi_tmp$DEobj$optim$bestval)
  
  res_hhs3_arsv_chi_tmp <- data.frame (HHS_region =3,
                                       pathogen2 = "fluA",
                                       chi=chi_vals[i],
                                       loglik=-res_hhs3_arsv_chi_tmp$DEobj$optim$bestva)
  res_hhs3_arsv_chi <-rbind(res_hhs3_arsv_chi,res_hhs3_arsv_chi_tmp)
  
}

write_csv(res_hhs3_arsv_chi,"res_hhs3_arsv_chi_loglik.csv",row.names=FALSE)

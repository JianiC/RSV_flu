# this script runs different functions to make sure that they are running alright

# first load prerequisites
source("./src.R", chdir = TRUE)

## specify the virus and HHS region

# make data ready for pomp
pomp_data_hhs8_brsv <- (
inc_data_add %.>%
  make_data_pomp_ready(., virus_combo = c("RSV", "fluB"), HHS_region = 8)
  )
load("res_hhs8_brsv_chi.rds")

best_past_est<-c(R01=past_est(res_hhs8_brsv_chi)$R01,
                          R02=past_est(res_hhs8_brsv_chi)$R02,
                          amplitude1=past_est(res_hhs8_brsv_chi)$amplitude1,
                          amplitude2 =past_est(res_hhs8_brsv_chi)$amplitude2,
                          tpeak1=past_est(res_hhs8_brsv_chi)$tpeak1,
                          tpeak2=past_est(res_hhs8_brsv_chi)$tpeak2,
                          rho1=past_est(res_hhs8_brsv_chi)$rho1,
                          rho2=past_est(res_hhs8_brsv_chi)$rho2,
                          psi=past_est(res_hhs8_brsv_chi)$psi)


# loading parameter constraints 

 # regular parameters for the full model
  # Npop variable will be updated based on the data
  rp_vals_def <- c(R01 = 1, gamma1=365./9, w1=2,
               R02 = 1, gamma2=365./3, w2=2,
               phi1=365/180, phi2=365/180, psi =1.0, chi=1.0,
               eta1=365., eta2=365.,rho1 = 0, rho2 = 0,
               amplitude1=0.0, tpeak1=0.0, amplitude2=0.0, tpeak2=0.0,
               pop=pomp_data_hhs8_brsv$N[1],
               mu=1/80)

res_hhs8_brsv_coinfect <- (
  DE_traj_match(df = pomp_data_hhs8_brsv, 
                param_constraints = co_infect_param_constraints, 
                params = rp_vals_def,
                ode_control = list(method = "ode23"),
                best_past_est = best_past_est,  
				hypo_name = "co_infect",
                hhs_reg = 8, 
                tot1_name = "RSV", 
                tot2_name = "fluB")
)

if(res_hhs8_brsv_coinfect$total2 == "fluB") message("Code itegration complete!")

save(res_hhs8_brsv_coinfect, file = "res_hhs8_brsv_coinfect.Rdata")



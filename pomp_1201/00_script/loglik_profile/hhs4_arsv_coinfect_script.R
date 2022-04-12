# this script runs different functions to make sure that they are running alright

# first load prerequisites
source("./src.R", chdir = TRUE)

## specify the virus and HHS region

# make data ready for pomp
pomp_data_hhs4_arsv <- (
inc_data_add %.>%
  make_data_pomp_ready(., virus_combo = c("RSV", "fluA"), HHS_region = 4)
  )


best_past_est<-data.frame(R01=1.513132,
                          R02=1.113691,
                          amplitude1=0.150284,
                          amplitude2=0.088219,
                          tpeak1=0.768802,
                          tpeak2=0.886405,
                          rho1=0.000150,
                          rho2=0.001069,
                          psi=0,
                          chi=0.8140189)


# loading parameter constraints 

 # regular parameters for the full model
  # Npop variable will be updated based on the data
  rp_vals_def <- c(R01 = 1, gamma1=365./9, w1=2,
               R02 = 1, gamma2=365./3, w2=2,
               phi1=365/180, phi2=365/180, psi =1, chi=1,
               eta1=365., eta2=365.,rho1 = 0, rho2 = 0,
               amplitude1=0.0, tpeak1=0.0, amplitude2=0.0, tpeak2=0.0,
               pop=pomp_data_hhs4_arsv$N[1],
               mu=1/80)

                          
res_hhs4_arsv_coinfect <- (
  DE_traj_match(df = pomp_data_hhs4_arsv, 
                param_constraints = co_infect_param_constraints, 
                params = rp_vals_def,
                ode_control = list(method = "ode23"), 
                hypo_name = "co_infect", 
                hhs_reg = 4, 
                tot1_name = "RSV", 
                tot2_name = "fluA")
)

if(res_hhs4_arsv_coinfect$total2 == "fluA") message("Code itegration complete!")

save(res_hhs4_arsv_coinfect, file = "res_hhs4_arsv_coinfect.Rdata")



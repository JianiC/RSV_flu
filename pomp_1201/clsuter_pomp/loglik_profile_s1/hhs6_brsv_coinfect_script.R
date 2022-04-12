# this script runs different functions to make sure that they are running alright

# first load prerequisites
source("./src.R", chdir = TRUE)

## specify the virus and HHS region

# make data ready for pomp
pomp_data_hhs6_brsv <- (
inc_data_add %.>%
  make_data_pomp_ready(., virus_combo = c("RSV", "fluB"), HHS_region = 6)
  )



  best_past_est<-data.frame( R01 = 1.502630, R02 = 1.133589,
                             amplitude1= 0.196019, amplititude2 = 0.147097,
                             tpeak1= 0.814616 , tpeak2=0.068305,
                             rho1=0.000462, rho2=0.000304,psi=0,chi=0.6) 
# loading parameter constraints 

 # regular parameters for the full model
  # Npop variable will be updated based on the data
  rp_vals_def <- c(R01 = 1, gamma1=365./9, w1=2,
               R02 = 1, gamma2=365./3, w2=2,
               phi1=365/180, phi2=365/180, psi =1.0, chi=1.0,
               eta1=365., eta2=365.,rho1 = 0, rho2 = 0,
               amplitude1=0.0, tpeak1=0.0, amplitude2=0.0, tpeak2=0.0,
               pop=pomp_data_hhs6_brsv$N[1],
               mu=1/80)

res_hhs6_brsv_coinfect <- (
  DE_traj_match(df = pomp_data_hhs6_brsv, 
                param_constraints = co_infect_param_constraints, 
                params = rp_vals_def,
                ode_control = list(method = "ode23"), 
				hypo_name = "co_infect",
                hhs_reg = 6, 
                tot1_name = "RSV", 
                tot2_name = "fluB")
)

if(res_hhs6_brsv_coinfect$total2 == "fluB") message("Code itegration complete!")

save(res_hhs6_brsv_coinfect, file = "res_hhs6_brsv_coinfect.Rdata")



# this script runs different functions to make sure that they are running alright

# first load prerequisites
source("./src.R", chdir = TRUE)

## specify the virus and HHS region

# make data ready for pomp
pomp_data_hhs6_brsv <- (
  inc_data_add %.>% 
  make_data_pomp_ready(., virus_combo = c("RSV", "fluB"), HHS_region = 6)
  )



   best_past_est<-c( R01 = 1.505710, R02 = 1.106169,
                             amplitude1= 0.194734, amplititude2 = 0.152132 ,
                             tpeak1= 0.816242 , tpeak2=0.051223,
                             rho1= 0.000455, rho2=0.000276,chi=0.75)

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

res_hhs6_brsv_chi <- (
  DE_traj_match(df = pomp_data_hhs6_brsv, 
                param_constraints = chi_param_constraints, 
                params = rp_vals_def,
                ode_control = list(method = "ode23"), 
                best_past_est = best_past_est,
                hypo_name = "chi", 
                hhs_reg = 6, 
                tot1_name = "RSV", 
                tot2_name = "fluB")
)

if(res_hhs6_brsv_chi$total2 == "fluB") message("Code itegration complete!")

save(res_hhs6_brsv_chi, file = "res_hhs6_brsv_chi.Rdata")



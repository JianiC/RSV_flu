# this script runs different functions to make sure that they are running alright

# first load prerequisites
source("./src.R", chdir = TRUE)

## specify the virus and HHS region

# make data ready for pomp
pomp_data_hhs6_arsv <- (
inc_data_add %.>%
  make_data_pomp_ready(., virus_combo = c("RSV", "fluA"), HHS_region = 6)
  )



best_past_est<-data.frame(R01=1.758388,
                          R02=1.180184,
                          amplitude1=0.133626,
                          amplitude2=0.141616,
                          tpeak1=0.867631,
                          tpeak2=0.936452,
                          rho1=0.000361,
                          rho2=0.000804,
                          psi=0.000010,
                          chi=0.765)


# loading parameter constraints 

 # regular parameters for the full model
  # Npop variable will be updated based on the data
  rp_vals_def <- c(R01 = 1, gamma1=365./9, w1=2,
               R02 = 1, gamma2=365./3, w2=2,
               phi1=365/180, phi2=365/180, psi =1, chi=1,
               eta1=365., eta2=365.,rho1 = 0, rho2 = 0,
               amplitude1=0.0, tpeak1=0.0, amplitude2=0.0, tpeak2=0.0,
               pop=pomp_data_hhs6_arsv$N[1],
               mu=1/80)

                          
res_hhs6_arsv_coinfect <- (
  DE_traj_match(df = pomp_data_hhs6_arsv, 
                param_constraints = co_infect_param_constraints, 
                params = rp_vals_def,
                ode_control = list(method = "ode23"), 
                hypo_name = "co_infect", 
                hhs_reg = 6, 
                tot1_name = "RSV", 
                tot2_name = "fluA")
)

if(res_hhs6_arsv_coinfect$total2 == "fluA") message("Code itegration complete!")

save(res_hhs6_arsv_coinfect, file = "res_hhs6_arsv_coinfect.Rdata")



# this script runs different functions to make sure that they are running alright

# first load prerequisites
source("./src.R", chdir = TRUE)

## specify the virus and HHS region

# make data ready for pomp
pomp_data_hhs3_brsv <- (
inc_data_add %.>%
  make_data_pomp_ready(., virus_combo = c("RSV", "fluB"), HHS_region = 3)
  )


# loading parameter constraints 

 # regular parameters for the full model
  # Npop variable will be updated based on the data
  rp_vals_def <- c(R01 = 1, gamma1=365./9, w1=2,
               R02 = 1, gamma2=365./3, w2=2,
               phi1=365/180, phi2=365/180, psi =1.0, chi=1.0,
               eta1=365., eta2=365.,rho1 = 0, rho2 = 0,
               amplitude1=0.0, tpeak1=0.0, amplitude2=0.0, tpeak2=0.0,
               pop=pomp_data_hhs3_brsv$N[1],
               mu=1/80)

res_hhs3_brsv_chi <- (
  DE_traj_match(df = pomp_data_hhs3_brsv, 
                param_constraints = chi_param_constraints, 
                params = rp_vals_def,
                ode_control = list(method = "ode23"), 
                hypo_name = "chi", 
                hhs_reg = 3, 
                tot1_name = "RSV", 
                tot2_name = "fluB")
)

if(res_hhs3_brsv_chi$total2 == "fluB") message("Code itegration complete!")

save(res_hhs3_brsv_chi, file = "res_hhs3_brsv_chi.Rdata")



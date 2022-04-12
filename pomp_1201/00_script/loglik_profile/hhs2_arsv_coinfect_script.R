# this script runs different functions to make sure that they are running alright

# first load prerequisites
source("./src.R", chdir = TRUE)

## specify the virus and HHS region

# make data ready for pomp
pomp_data_hhs2_arsv <- (
inc_data_add %.>%
  make_data_pomp_ready(., virus_combo = c("RSV", "fluA"), HHS_region = 2)
  )




# loading parameter constraints 

 # regular parameters for the full model
  # Npop variable will be updated based on the data
  rp_vals_def <- c(R01 = 1, gamma1=365./9, w1=2,
               R02 = 1, gamma2=365./3, w2=2,
               phi1=365/180, phi2=365/180, psi =1, chi=1,
               eta1=365., eta2=365.,rho1 = 0, rho2 = 0,
               amplitude1=0.0, tpeak1=0.0, amplitude2=0.0, tpeak2=0.0,
               pop=pomp_data_hhs2_arsv$N[1],
               mu=1/80)

# loading parameter constraints

  best_past_est<-data.frame( R01 = 1.681374, R02 = 1.222560 ,
                             amplitude1= 0.16976, amplititude2 =  0.121962,
                             tpeak1= 0.839223 , tpeak2=0.030687,
                             rho1=0.000257, rho2=0.000742,
                             psi=0.999999,chi=0.6)                          
                          
res_hhs2_arsv_coinfect <- (
  DE_traj_match(df = pomp_data_hhs2_arsv, 
                param_constraints = co_infect_param_constraints, 
                params = rp_vals_def,
                ode_control = list(method = "ode23"), 
                hypo_name = "co_infect", 
                hhs_reg = 2, 
                tot1_name = "RSV", 
                tot2_name = "fluA")
)

if(res_hhs2_arsv_coinfect$total2 == "fluA") message("Code itegration complete!")

save(res_hhs2_arsv_coinfect, file = "res_hhs2_arsv_coinfect.Rdata")



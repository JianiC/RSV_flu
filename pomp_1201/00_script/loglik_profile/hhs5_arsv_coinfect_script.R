# this script runs different functions to make sure that they are running alright

# first load prerequisites
source("./src.R", chdir = TRUE)

## specify the virus and HHS region

# make data ready for pomp
pomp_data_hhs5_arsv <- (
inc_data_add %.>%
  make_data_pomp_ready(., virus_combo = c("RSV", "fluA"), HHS_region = 5)
  )

load("./res_hhs5_arsv_coinfect.Rdata")
best_past_est<-data.frame(past_est(res_hhs5_arsv_coinfect)$R01,
                          past_est(res_hhs5_arsv_coinfect)$R02,
                          past_est(res_hhs5_arsv_coinfect)$amplitude1,
                          past_est(res_hhs5_arsv_coinfect)$amplitude2,
                          past_est(res_hhs5_arsv_coinfect)$tpeak1,
                          past_est(res_hhs5_arsv_coinfect)$tpeak2,
                          past_est(res_hhs5_arsv_coinfect)$rho1,
                          past_est(res_hhs5_arsv_coinfect)$rho2,
                          past_est(res_hhs5_arsv_coinfect)$psi,
                          past_est(res_hhs5_arsv_coinfect)$chi)


# loading parameter constraints 

 # regular parameters for the full model
  # Npop variable will be updated based on the data
  rp_vals_def <- c(R01 = 1, gamma1=365./9, w1=2,
               R02 = 1, gamma2=365./3, w2=2,
               phi1=365/180, phi2=365/180, psi =1, chi=1,
               eta1=365., eta2=365.,rho1 = 0, rho2 = 0,
               amplitude1=0.0, tpeak1=0.0, amplitude2=0.0, tpeak2=0.0,
               pop=pomp_data_hhs5_arsv$N[1],
               mu=1/80)

# loading parameter constraints

                          
                          
res_hhs5_arsv_coinfect <- (
  DE_traj_match(df = pomp_data_hhs5_arsv, 
                param_constraints = co_infect_param_constraints, 
                params = rp_vals_def,
                ode_control = list(method = "ode23"), 
                hypo_name = "co_infect", 
                hhs_reg = 5, 
                tot1_name = "RSV", 
                tot2_name = "fluA")
)

if(res_hhs5_arsv_coinfect$total2 == "fluA") message("Code itegration complete!")

save(res_hhs5_arsv_coinfect, file = "res_hhs5_arsv_coinfect.Rdata")



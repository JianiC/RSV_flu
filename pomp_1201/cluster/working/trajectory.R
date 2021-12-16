rp_vals <- c(R01=param_vector$R01, gamma1=365./9, w1=1/1,
             R02=param_vector$R0_B, gamma2=365./3, w2=1/1,
             amplitude1 = param_vector$amplitude1,
             amplitude2 = param_vector$amplitude2,
             tpeak1 = param_vector$tpeak1,
             tpeak2 = param_vector$tpeak2,
             phi1=param_vector$phi1, phi2=param_vector$phi2,
             psi=param_vector$psi, chi=param_vector$chi,
             eta_A=365., eta_B=365.,
             rho1=param_vector$rho1, rho2=param_vector$rho2,
             sigmaSE=0.0000,
             pop=Npop,
             mu=1/80)
ic_vals <- SIRS2_independent_endemic_equilibrium(rp_vals)
params_all <- c(rp_vals,ic_vals)


test_traj <- trajectory(object = hhs1_a_rsv_po, params = params_all, format = "d", method = "ode23")
test_traj %.>% 
  slice(., 2:n()) %.>% 
  select(., -`.id`) %.>% 
  gather(., "comp", "count", -time)
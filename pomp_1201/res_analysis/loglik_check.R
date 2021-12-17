## get loglikelihood for a estimate

fit_par<-get_rp_vals(data=inc_data,res_hhs1_arsv_coinfect)
model_variables = c("S1_S2","I1_S2", "C1_S2","R1_S2",
                    "S1_I2","I1_I2", "C1_I2","R1_I2",
                    "S1_C2","I1_C2", "C1_C2","R1_C2",
                    "S1_R2","I1_R2", "C1_R2","R1_R2")
model_params = c("R01","R02", "gamma1", "gamma2",
                 "w1","w2","eta1", "eta2", "phi1", "phi2",
                 "psi","chi", "rho1", "rho2",
                 "amplitude1", "tpeak1", "amplitude2", "tpeak2", "pop","mu")
pomp_data_hhs1_arsv <- (
  inc_data %.>% 
    make_data_pomp_ready(., virus_combo = c("RSV", "fluA"), HHS_region = 1)
)


pomp_objfun <- (
  make_pomp(pomp_data_hhs1_arsv,time_start_sim = -100) %.>%
    # define the objective function 
    traj_objfun(., 
                est = names(co_infect_param_constraints$lower), 
                params = fit_par ,
                dmeasure=dmeas_poisson,fail.value = 1e20,
                statenames = c(model_variables, c("Kprim1","Ksec1","Kprim2","Ksec2","K1", "K2")),
                paramnames = c(model_params))
)

pomp_objfun(fit_par)
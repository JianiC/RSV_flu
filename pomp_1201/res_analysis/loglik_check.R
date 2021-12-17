## get loglikelihood for a estimate

params<-get_rp_vals(data=incdata, virus_combo,HHS_region,res_hhs)

pomp_objfun <- (
  make_pomp(inc_data_add, time_start_sim = -100) %.>%
    # define the objective function 
    traj_objfun(., 
                est = names(param_constraints$lower), 
                params = params, fail.value = 1e20, 
                ode_control = ode_control)
)
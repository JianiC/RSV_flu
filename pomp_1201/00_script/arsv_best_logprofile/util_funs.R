# this script contains util functions for pomp related analysis - estimation, simulation and prodiction

# this 'convenience' function modifies vector of default params to set hypothesis specific param estimates
sim_p_vals <- function(estm_vect, default_p_vals = param_vals_est) {
  
  replace_these <- names(default_p_vals)[names(default_p_vals) %in% names(estm_vect)] 
  
  tmp_p_vals <- default_p_vals
  tmp_p_vals[replace_these] <- estm_vect[replace_these]
  
  tmp_p_vals
  
}


# a wrapper to unlist and unname a list object.
un_list_name <- function(x) {
  x %.>% 
    unlist(.) %.>% 
    unname(.)
}  


# operator for easy analysis
`%nin%` <- Negate(`%in%`)



# wrapper - this function produces pomp style objective function for DEoptim()
DE_traj_objfun <- function(x, objfun, est, 
                           seed = 986747881L) {
  
  # produce an empty vector
  par_v <- rep(NA, times = length(est))
  # name the vector
  names(par_v) <- est
  
  # assign the guess to a new internal vector
  par_v[est] <- un_list_name(x)
  
  # supply to the pomp objective function
  objfun(par = par_v)
}


#  this is a wrapper function that carries out parameter estimation using DEoptim()
DE_traj_match <- function(param_constraints,
                          params = p_vals,
                          ninit = np_val, ode_control = NULL, 
                          hypo_name, hhs_reg, tot1_name, tot2_name, 
                          best_past_est = NULL,
                          seed = 986747881L,
                          other_DE_controls = my_controls,
                          ...) {
  
  message(cat(c("Est: ", names(param_constraints$lower))))
  # browser()
  # generate a pomp objective function - this step also includes defining the pomp object
  # NOTE: names of the parameters estimated are taken from the lower constraint vector   
  pomp_objfun <- (
    make_pomp(...) %.>%
      # define the objective function 
      traj_objfun(., 
                  est = names(param_constraints$lower), 
                  params = params, fail.value = 1e20, 
                  ode_control = ode_control)
  )
  # browser()
  # generate a grid of initial guesses
  if(is.null(best_past_est)) {
    
    init_guess_grid <- sobol_design(lower = param_constraints$lower, 
                                    upper = param_constraints$upper, 
                                    nseq = ninit)   
  } else {
    
    init_guess_grid <- sobol_design(lower = param_constraints$lower, 
                                    upper = param_constraints$upper, 
                                    nseq = ninit)%>% 
    bind_rows(., best_past_est) %.>% 
      replace(., is.na(.), 1) 
    
  }
  
  # browser()
  
  # set seed for reproducible parallel computation
  set.seed(986474881L)
  
  RNGkind("L'Ecuyer-CMRG")
  
  # set multi-core cluster for parallel computation optimal solution
  #no_cores <- detectCores() - 1
  no_cores <- 93
  
  registerDoParallel(cores = no_cores)  
  
  cl <- makeCluster(no_cores, type="FORK")
  
  
  # feed all this info to the evolutionary optimizer
  DEobj <- DEoptim(fn = DE_traj_objfun, 
                   est = names(param_constraints$lower), 
                   objfun = pomp_objfun,
                   seed = seed,
                   lower = param_constraints$lower, 
                   upper = param_constraints$upper, 
                   control = c(other_DE_controls, 
                               list(cluster = cl, 
                                    NP = nrow(init_guess_grid), 
                                    initialpop = init_guess_grid %.>% 
                                      as.matrix(.))
                   )
  ) 
  
  stopCluster(cl)
  
  # collect results here
  result <- list(initial_pop = init_guess_grid, 
                 DEobj = DEobj, 
                 Hypothesis = hypo_name, 
                 HHS = hhs_reg, 
                 total1 = tot1_name, 
                 total2 = tot2_name)
  
  # written the result
  result  
  
}


# this function prepares incidence data for pomp
make_data_pomp_ready <- function(data = inc_data_add, virus_combo = c("RSV", "fluA"), HHS_region = 1) {
  #browser()
  data %.>% 
    filter(., HHS_REGION == HHS_region & virus %in% virus_combo) %.>% 
    mutate(., 
           virus = ifelse(virus == virus_combo[1], "total1", "total2")) %.>% 
    select(., -c(HHS_REGION)) %.>% 
    spread(., key = virus, value = cases)
}

# make_data_pomp_ready()


## get past estimate
past_est<-function(res_hhs){
  
  param_vector<-as.data.frame(res_hhs$DEobj$optim$bestmem)%>%
    cbind(row.names(.))%>%
    pivot_wider(names_from = `row.names(.)`, values_from = `res_hhs$DEobj$optim$bestmem`)
  
  
  return(param_vector)
  
  
}





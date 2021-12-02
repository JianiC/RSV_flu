## model setup


## create a pomp object for RSV and influenza co-circulation without cross protection
library(dplyr)
library(pomp)
library(GA)
library(tictoc)
library(tidyr)
library(foreach)
library(doParallel)
############################################################################
# Code to define  process model -------------------------------------------
############################################################################
# 1-step function of model
rproc_euler_multinomial <- Csnippet("
  
  double betaA, betaB, foiA, foiB, seas_A, seas_B, N; //dw foiB
  double rate[16], trans[16];

  // white noise (extra-demographic stochasticity, 
  // currently asssumed the same for both types(?))
  // dw = rgammawn(sigmaSE,dt);
    
  // sinusoidal seasonality (for both types)
  seas_A = 1 + amplitude_A*cos(2.*3.1415*(t-tpeak_A)/1.); // in units of years
  seas_B = 1 + amplitude_B*cos(2.*3.1415*(t-tpeak_B)/1.); // in units of years

  betaA = R0_A*seas_A*gamma_A;
  betaB = R0_B*seas_B*gamma_B; 


  // forces of infection
  N = S + I_A + I_B + C_A + C_B + R_A + R_B + I_AB + I_BA + R_AB;
  
  foiA = betaA*(I_A + I_AB + eta_A)/N; 
  foiB = betaB*(I_B + I_BA + eta_B)/N;
  
  // Out of S
  rate[0] = foiA; // *dw/dt;        // infection of S with A
  rate[1] = foiB; // *dw/dt;        // infection of S with B
  
  reulermultinom(2,S,&rate[0],dt,&trans[0]);
  
  // Out of I_A
  rate[2] = gamma_A;                // recovery of I_A
  
  reulermultinom(1,I_A,&rate[2],dt,&trans[2]);

  // Out of I_B
  rate[3] = gamma_B;                // recovery of I_B
  
  reulermultinom(1,I_B,&rate[3],dt,&trans[3]);

  
  // Out of R_A
  rate[4] = foiB; //*dw/dt;        // infection of R_A with B
  rate[5] = w_A;                   // loss of immunity of R_A 
  
  reulermultinom(2,R_A,&rate[4],dt,&trans[4]);

  
  // Out of R_B
  rate[6] = foiA; //*dw/dt;       // infection of R_B with A
  rate[7] = w_B;                  // loss of immunity of R_B
  
  reulermultinom(2,R_B,&rate[6],dt,&trans[6]);

  
  // Out of I_BA
  rate[8] = gamma_B;              // recovery of I_BA (recovery from B infection)
  
  reulermultinom(1,I_BA,&rate[8],dt,&trans[8]);

  // Out of I_AB
  rate[9] = gamma_A;              // recovery of I_AB (recovery from A infection)
  
  reulermultinom(1,I_AB, &rate[9],dt,&trans[9]);

  // Out of R_AB
  rate[10] = w_A;                 // loss of immunity of R_AB to A
  rate[11] = w_B;                 // loss of immunity of R_AB to B
  
  reulermultinom(2,R_AB,&rate[10],dt,&trans[10]);
  
  // Out of C_A
  rate[12] = phi_A;                      // loss of crossimmunity of C_A to R_A
  rate[13] = (1-chi_BA)*foiB; //*dw/dt;  // infection of C_A with B
  
  reulermultinom(2,C_A,&rate[12],dt,&trans[12]);


  // Out of C_B
  rate[14] = phi_B;                      // loss of crossimmunity of C_B to R_B
  rate[15] = (1-chi_AB)*foiA; //*dw/dt;  // infection of C_B with A
  
  reulermultinom(2,C_B,&rate[14],dt,&trans[14]);


  S   += trans[5]  + trans[7] - trans[0] - trans[1];
  I_A += trans[0]  - trans[2];
  I_B += trans[1]  - trans[3];
  C_A += trans[2]  - trans[12] - trans[13];
  C_B += trans[3]  - trans[14] - trans[15]; 
  R_A += trans[12] + trans[11] - trans[4] -  trans[5];
  R_B += trans[14] + trans[10] - trans[6] -  trans[7];
  I_BA += trans[4] + trans[13] - trans[8];
  I_AB += trans[6] + trans[15] - trans[9];
  R_AB += trans[8] + trans[9]  - trans[10] - trans[11];

  //W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise

  K_A += trans[2] + trans[9];           // true incidence of A
  K_B += trans[3] + trans[8];           // true incidence
  
")

# C snippet for initial condition specification ---------------------------

rinit_ee <- Csnippet("
  // how best to initialise population to ensure population size is correct
  // what is this baryometric thing
  // variables iS_A etc refer to endemic equilibrium of single type model
  double iS_A = 1/R0_A;
  double iS_B = 1/R0_B;

  
  double wA_prime = w_A + phi_A;
  double wB_prime = w_B + phi_B;

  double k_A = gamma_A/(wA_prime);
  double k_B = gamma_B/wB_prime;

  double r_A = phi_A/w_A;
  double r_B = phi_B/w_B;


  double iI_A = (1-iS_A)/(1 + k_A + k_A*r_A);
  double iI_B = (1-iS_B)/(1 + k_B + k_B*r_B);

  double iC_A = k_A*iI_A;
  double iC_B = k_B*iI_B;
  double iR_A = r_A*iC_A;
  double iR_B = r_B*iC_B;

  
  S = nearbyint(pop*iS_A*iS_B);
  I_A = nearbyint(pop*iI_A*iS_B);
  I_B = nearbyint(pop*iS_A*iI_B);
  C_A = nearbyint(pop*iC_A*iS_B);
  C_B = nearbyint(pop*iS_A*iC_B);
  R_A = nearbyint(pop*iR_A*iS_B);
  R_B = nearbyint(pop*iS_A*iR_B);
  I_AB = nearbyint(pop*iI_A*iR_B);
  I_BA = nearbyint(pop*iR_A*iI_B);
  R_AB = nearbyint(pop*iR_A*iR_B);

  W = 0;
  K_A = 0;
  K_B = 0;

")

det_skel <- Csnippet("

  double betaA, betaB, foiA, foiB, seas_A, seas_B, N; //dw foiB
  double rate[16];

  // white noise (extra-demographic stochasticity, 
  // currently asssumed the same for both types(?))
  //dw = rgammawn(sigmaSE,dt);
    
  // sinusoidal seasonality (for both types)
  seas_A = 1 + amplitude_A*cos(2.*3.1415*(t-tpeak_A)/1.); // in units of years
  seas_B = 1 + amplitude_B*cos(2.*3.1415*(t-tpeak_B)/1.); // in units of years

  betaA = R0_A*seas_A*gamma_A;
  betaB = R0_B*seas_B*gamma_B; 


  // forces of infection
  N = S + I_A + I_B + C_A + C_B + R_A + R_B + I_AB + I_BA + R_AB;
  
  foiA = betaA*(I_A + I_AB + eta_A)/N; 
  foiB = betaB*(I_B + I_BA + eta_B)/N;
  
  // Out of S
  rate[0] = S*foiA; // *dw/dt;               // infection of S with A
  rate[1] = S*foiB; // *dw/dt;               // infection of S with B
  
  // Out 
  rate[2] = I_A*gamma_A;                     // recovery of I_A
  rate[3] = I_B*gamma_B;                     // recovery of I_B
  rate[4] = R_A*foiB; //*dw/dt;              // infection of R_A with B
  rate[5] = R_A*w_A;                         // loss of immunity of R_A 
  rate[6] = R_B*foiA; //*dw/dt;              // infection of R_B with A
  rate[7] = R_B*w_B;                         // loss of immunity of R_B
  rate[8] = I_BA*gamma_B;                    // recovery of I_BA (recovery from B infection)
  rate[9] = I_AB*gamma_A;                    // recovery of I_AB (recovery from A infection)
  rate[10] = R_AB*w_A;                       // loss of immunity of R_AB to A
  rate[11] = R_AB*w_B;                       // loss of immunity of R_AB to B

  rate[12] = C_A*phi_A;                      // loss of crossimmunity of C_A to R_A
  rate[13] = C_A*(1-chi_BA)*foiB; //*dw/dt;  // infection of C_A with B
  
  rate[14] = C_B*phi_B;                      // loss of crossimmunity of C_B to R_B
  rate[15] = C_B*(1-chi_AB)*foiA; //*dw/dt;  // infection of C_B with A
  
  
  // transitions 
  DS    = rate[5]  +  rate[7] - rate[0] - rate[1];
  DI_A  = rate[0]  -  rate[2];
  DI_B  = rate[1]  -  rate[3];
  DC_A  = rate[2]  -  rate[12]  -  rate[13];
  DC_B  = rate[3]  -  rate[14]  -  rate[15];
  DR_A  = rate[12] +  rate[11]  -  rate[4]   - rate[5];
  DR_B  = rate[14] +  rate[10]  -  rate[6]   - rate[7];
  DI_BA = rate[4]  +  rate[13]  -  rate[8];
  DI_AB = rate[6]  +  rate[15]  -  rate[9];
  DR_AB = rate[8]  +  rate[9]   -  rate[10]  - rate[11];

  //W = (dw - dt)/sigmaSE;  // standardized i.i.d. white noise

  DK_A = rate[2] + rate[9];           // true incidence of A
  DK_B = rate[3] + rate[8];           // true incidence of B
  
")
############################################################################
# Code to define estimation components of  model 
############################################################################

# Define likelihood function ----------------------------------------------

# assume reporting probability the same (is this an issue?)
dmeas_poisson <- Csnippet("
  double lik_A, lik_B;
  
  if(ISNA(total_a)) {
    lik_A = (give_log) ? 0:1;
  } else {
      lik_A = dpois(nearbyint(total_a), rho_A*K_A, give_log); 
  }
  
  if(ISNA(total_a)) {
    lik_B = (give_log) ? 0:1;
  } else {
      lik_B = dpois(nearbyint(total_b), rho_B*K_B, give_log); 
  }
  
  lik = lik_A + lik_B;

")



# Define process simulator for observations  ------------------------------
rmeas_poisson <- Csnippet("

  total_a = rpois(rho_A*K_A);
  total_b = rpois(rho_B*K_B);
  
")


############################################################################
# assign values for model parameters and initial conditions -----------------
############################################################################

# Define the pomp model object --------------------------------------------


make_pomp_model <- function(df, time_start_sim, dt=1) {
  
  
  # Model parameter names
  model_params = c("R0_A","R0_B", "gamma_A", "gamma_B",
                   "w_A", "w_B","eta_A", "eta_B", "phi_A", "phi_B",
                   "chi_BA", "chi_AB", "rho_A", "rho_B", "sigmaSE",
                   "amplitude_A", "tpeak_A", "amplitude_B", "tpeak_B", "pop")
  
  # Model variable names
  model_variables = c("S", "I_A", "I_B", "R_A", "R_B", "C_A", "C_B", 
                      "I_BA", "I_AB", "R_AB")
  
  # Model initial conditions parameter names
  model_ic_params = paste(model_variables, "_0", sep="")
  
  process_model <- euler(rproc_euler_multinomial, delta.t = dt/365.25)
  
  dmeas <- dmeas_poisson
  rmeas <- rmeas_poisson
  
  po <- pomp(
    data = df,
    times = "time",
    t0 = time_start_sim,
    obsnames = c("total_a","total_b"),
    rprocess = process_model,
    skeleton = vectorfield(det_skel),
    dmeasure = dmeas,
    rmeasure = rmeas, 
    rinit = rinit_ee,
    accumvars =c ("K_A", "K_B", "W"),
    statenames = c(model_variables, c("K_A", "K_B", "W")),
    #cfile = "2strain_flu_model",
    #cdir = getwd(),
    paramnames = c(model_params, model_ic_params)
  )
  
  return(po)
}

# Load relevant files and packages
set.seed(594709947L)

stopifnot(packageVersion("pomp")>="2.0.9.1")

# specify the kind of seed for parallelizingwhen using mclapply 
RNGkind("L'Ecuyer-CMRG")

t_start <- -100 

TX_data_3s<-read.csv("test_12_13.csv")

TX_data_3s %>% 
  make_pomp_model(time_start_sim=t_start) -> pomp_sirs_3s
###############################################################
print("run pomp likelihood function!")

## loglikelihood function

##############################################################################################################
################################ Objective function definitions ##############################################
##############################################################################################################
# objective function for 3 seasons
# This function calculates the posssion log-liklihood | condtional  log-liklihood 
# for a target-data:model pair 

of_3s <- function(par, pomp_object = pomp_sirs_3s, est, params_all_int = params_all, 
                  give_conditional_log_lik = FALSE, target_data = inc_data, 
                  fail_value = -1e9) {
  
  # replace values for the parameters of interest
  if(class(par) %in% c("numeric")) {
    params_all_int[est] <- par[est]
  } else{
    params_all_int[est] <- unlist(par[,est])  
  }
  # browser()
  # use globally defined pomp object and replace param vector
  pomp_object %>% 
    pomp(params  = params_all_int) %>% 
    # generate trajectories using the replaced parameters
    trajectory(include.data = FALSE, format = "d", method = "ode45") %>%
    # scale the "true incidence" by the reporting probabilities situated in params_all_int
    mutate(K_A = ifelse(time < 2.5, NA, K_A), 
           total_a = params_all_int["rho_A"]*K_A,
           K_B = ifelse(time < 2.5, NA, K_B), 
           total_b = params_all_int["rho_B"]*K_B) %>% 
    # select relevant state veariables: type specific new scaled cases of flu
    select(time, total_a, total_b) %>% 
    # store this in a long-format tibble 
    gather(key = flu_type, value = st_cases, -time) %>% 
    # join this with the target data: to calculate poisson log-likelihood
    right_join(target_data, by = c("time", "flu_type")) %>% 
    # calculate poisson log-density: P(data_t|model_t)
    mutate(poisson_log_density = dpois(x = obs_cases, lambda = st_cases, log = TRUE)) -> log_density_data 
  
  # This feature is for potential post-hoc analysis: if conditional log-liklihoods are called for
  if(give_conditional_log_lik == TRUE) {
    
    log_density_data %>% 
      # NA's are given a conditional log-density of 0
      replace_na(list(poisson_log_density = 0)) %>% 
      # Following two steps are used to generate type-specific conditional log-liklihood
      group_by(flu_type) %>% 
      mutate(conditional_loglikelihood = cumsum(poisson_log_density)) -> result
    
  } else {
    
    log_density_data %>% 
      # calculate the poisson log-likelihood : FOR THE OPTIMIZER
      select(poisson_log_density) %>% 
      summarise_all(sum, na.rm = TRUE) %>% 
      mutate(poisson_log_density = ifelse(is.finite(poisson_log_density) == TRUE & 
                                            poisson_log_density < 0, poisson_log_density, 
                                          fail_value)) %>% 
      unlist() -> result
  }
  
  print(result)
  print(par)
  # browser()
  return(result)  
  
}



# This function is defined to provide pomp-GA integration
of_ga_3s <- function(x, est = est_these) {
  
  x<- unname(unlist(x))
  # browser()
  of_3s(par = c(R0_B = x[1], R0_A = x[2], 
                amplitude_A = x[3], amplitude_B = x[4], 
                tpeak_A = x[5], tpeak_B = x[6],  
                rho_A = x[7], rho_B = x[8]), est = est, 
        params_all_int = neutral_params_all)
}

#######################################################
print("start to set parameter to guess best value")

## set parameters and vairalbes

# Define the date in a long format : To be used in the objective function
TX_data_3s %>%
  gather(key = flu_type, value = obs_cases, -time) -> inc_data

# set a default vector of parameters: requirement of the objective function ::: Can this be improved(?)
# NOTE: Keep the defaults to be at R0s of 1 and neutral model 
# NOTE: phi value shave been fixed to 365/30 for defaults, mostly for convenience
rp_vals <- c(R0_A = 1, gamma_A=365./9, w_A=1./1.,
             R0_B = 1, gamma_B=365./2.5, w_B=1./1.,
             phi_A=365/30., phi_B=365/30.,
             chi_BA=0, chi_AB=0, eta_A=365., eta_B=365.,
             rho_A=0, rho_B=0, 
             sigmaSE=0.0000, 
             amplitude_A=0, tpeak_A=0, amplitude_B=0, tpeak_B=0,
             pop=2.6e7)

SIRS2_independent_endemic_equilibrium <- function(params){
  S_A <- 1/params[["R0_A"]]
  S_B <- 1/params[["R0_B"]]
  
  k_A <- params[["gamma_A"]]/(params[["w_A"]] + params[["phi_A"]])
  k_B <- params[["gamma_B"]]/(params[["w_B"]] + params[["phi_B"]])
  
  r_A <-  params[["phi_A"]]/params[["w_A"]]
  r_B <-  params[["phi_B"]]/params[["w_B"]]
  
  I_A <- (1-S_A)/(1+k_A + k_A*r_A)
  I_B <- (1-S_B)/(1+k_B + k_B*r_B)
  
  ee_A <- c(S_0=S_A, I_0=I_A,  C_0=k_A*I_A, R_0=r_A*k_A*I_A)
  ee_B <- c(S_0=S_B, I_0=I_B,  C_0=k_B*I_B, R_0=r_B*k_B*I_B)
  
  # Note only  classes: model assumes coinfection (I_A * I_B) is neglible
  return(c(S_0=ee_A[["S_0"]]*ee_B[["S_0"]],
           I_A_0 = ee_A[["I_0"]]*ee_B[["S_0"]],
           I_B_0 = ee_A[["S_0"]]*ee_B[["I_0"]],
           C_A_0 = ee_A[["C_0"]]*ee_B[["S_0"]],
           C_B_0 = ee_A[["S_0"]]*ee_B[["C_0"]],
           R_A_0 = ee_A[["R_0"]]*ee_B[["S_0"]],
           R_B_0 = ee_A[["S_0"]]*ee_B[["R_0"]],
           I_BA_0 = ee_A[["R_0"]]*ee_B[["I_0"]],
           I_AB_0 = ee_A[["I_0"]]*ee_B[["R_0"]],
           R_AB_0 = ee_A[["R_0"]]*ee_B[["R_0"]]))
}

ic_vals <- SIRS2_independent_endemic_equilibrium(rp_vals)
params_all <- c(rp_vals,ic_vals)

# setting the duration of crossprotection to zero, implies phi ~ Inf 
# This is set to have an infinitesimal duration of crossprotection
neutral_params_all <- params_all

# parameters of interest
est_these <- c("R0_B", "R0_A", "amplitude_A", "amplitude_B", "tpeak_A", "tpeak_B", 
               "rho_A", "rho_B")

# create a profile design 
pd <- sobol_design(lower = c(R0_B = 0.5, R0_A = 0.5, amplitude_A = 1e-10, 
                            amplitude_B = 1e-10, tpeak_A = 1/365.25, tpeak_B = 1/365.25,  
                            rho_A = 1e-10, rho_B = 1e-10), 
                  upper = c(R0_B = 10, R0_A = 10, 
                            amplitude_A = 1, amplitude_B = 1, 
                            tpeak_A = 1, tpeak_B = 1, 
                            rho_A = 0.01, rho_B = 0.01), nseq = 500)

# combine the dataframe of initial guesses with the mle from previous analysis
pd -> pd_with_fix_w_neutral_mle



# genetic algorithm set up to run in parallel
tic()
no_cores <- detectCores() - 1 

registerDoParallel(cores=no_cores)  

cl <- makeCluster(no_cores, type="FORK")

GA_neutral <- ga(type = "real-valued", 
                 fitness = of_ga_3s,
                 lower = c(R0_B = 0.5, R0_A = 0.5, amplitude_A = 1e-10, 
                           amplitude_B = 1e-10, tpeak_A = 1/365.25, tpeak_B = 1/365.25,  
                           rho_A = 1e-10, rho_B = 1e-10), 
                 upper = c(R0_B = 10, R0_A = 10, 
                           amplitude_A = 1, amplitude_B = 1, 
                           tpeak_A = 1, tpeak_B = 1, 
                           rho_A = 0.01, rho_B = 0.01),
                 elitism = base::max(1, round(100*.1)), 
                 popSize = 500, maxiter = 10000, run = 500, 
                 suggestions = pd_with_fix_w_neutral_mle,
                 optim = TRUE,
                 optimArgs = list(poptim = 0.1, 
                                  control = list(maxit = 10000)),
                 seed = 12345, parallel = cl)

stopCluster(cl)
toc() -> ga_took


test_oldresult_null<- list(it_took  = c(paste((ga_took$toc - ga_took$tic)/3600, "hrs"), 
                                    paste((ga_took$toc - ga_took$tic)/60, "mins"), 
                                    paste((ga_took$toc - ga_took$tic), "secs")), 
                       GAobj = GA_neutral)
toc()

save(test_oldresult_null, file = "test_oldresult_null.Rdata")

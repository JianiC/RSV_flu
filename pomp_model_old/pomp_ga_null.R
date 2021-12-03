## run pomp model with multiple setting
library(dplyr)
library(pomp)
library(GA)
library(tictoc)
library(tidyr)
library(foreach)
library(doParallel)
############################################################################
## load required parm and dataframes

## which region?
HHS_n=1


pop<-read.csv("pop_HHSRegion.csv")
Npop=pop[pop$HHS_region == HHS_n,3] ##population param
case_data<-read.csv("pomp_HHSregion_case.csv") ## case count data
data<-case_data%>%
  filter(HHS_REGION==HHS_n)%>%
  select(time,RSVpos,fluApos)%>%
  dplyr::rename(total_a=RSVpos,total_b=fluApos)%>%
  add_row(time = 2.5, total_a = NA, total_b=NA, .before = 1)

pd<-read.csv("pd.csv") ## sobol design
#############################################################################

## POMP model
## use A to stand for RSV,
## use B to stand for IAV/IBV
#############################################################################



rproc_euler_multinomial <- Csnippet("
  
  double betaA, betaB, foiA, foiB, seas_A, seas_B, N; 
  
  double rate[20], trans[20];

  // white noise (extra-demographic stochasticity, 
  // currently asssumed the same for both types(?))
  // dw = rgammawn(sigmaSE,dt);
    
  // sinusoidal seasonality (for both types)
  seas_A = 1 + amplitude_A*cos(2.*3.1415*(t-tpeak_A)/1.); // in units of years
  seas_B = 1 + amplitude_B*cos(2.*3.1415*(t-tpeak_B)/1.); // in units of years

  betaA = R0_A*seas_A*gamma_A;
  betaB = R0_B*seas_B*gamma_B; 


  // forces of infection
  N = SA_SB + IA_SB + CA_SB + RA_SB + SA_IB + SA_CB + SA_RB + IA_IB + RA_IB  + IA_RB + RA_RB;

  foiA = betaA*(IA_SB + IA_IB + IA_RB + eta_A)/N; 
  foiB = betaB*(SA_IB + IA_IB + RA_IB + eta_B)/N;
  
  // Out of SA_SB
  rate[0] = foiA; // *dw/dt;        // infection of SA_SB to IA_SB
  rate[1] = foiB; // *dw/dt;        // infection of SA_SB ti SA_IB
  
  reulermultinom(2,SA_SB,&rate[0],dt,&trans[0]);
  
  // Out of IA_SB
  rate[2] = gamma_A;                // recovery of IA_SB to CA_SB
  rate[3] = lamda * foiB;         // infection of IA_SB to IA_IB 
  
  reulermultinom(2,IA_SB,&rate[2],dt,&trans[2]);

  // Out of SA_IB
  rate[4] = gamma_B;                // recovery of SA_IB to SA_CB
  rate[5] = lamda * foiA;         // infection of SA_IB to IA_IB
  
  reulermultinom(2,SA_IB,&rate[4],dt,&trans[4]);

  // Out of CA_SB
  
  rate[6] = phi;                    //  loss cross-immunity of CA_SB to RA_SB
  rate[7] = chi * foiB;          // infectious of CA_SB to RA_IB

  reulermultinom(2,CA_SB,&rate[6],dt,&trans[6]);
  
  // Out of SA_CB
  
  rate[8] = phi;                    //  loss immunity of SA_CB to SA_RB
  rate[9] = chi * foiA;          // infectious of SA_CB to IA_RB

  reulermultinom(2,CA_SB,&rate[8],dt,&trans[8]);
  
  // Out of RA_SB
  
  rate[10] = w;                  // loss of immunity to SA_SB
  rate[11] = foiB;                 // infectious of RA_SB to RA_IB

  reulermultinom(2,RA_SB,&rate[10],dt,&trans[10]);
  
  // Out of SA_RB
  
  rate[12] = w;                  // loss of immunity to SA_SB
  rate[13] = foiA;                 // infectious of SA_RB to IA_RB

  reulermultinom(2,RA_SB,&rate[12],dt,&trans[12]);
  
  
  // out of RA_IB
  
  rate[14] = gamma_B;          // recover from RA_IB to RA_RB
  
  reulermultinom(1,RA_IB,&rate[14],dt,&trans[14]);
  
  
  // out of IA_RB

  rate[15] = gamma_A;          // recover from IA_RB to RA_RB
  
  reulermultinom(1,IA_RB,&rate[15],dt,&trans[15]);
  
  // out of RA_RB
  
  rate[16] = w;       // loss of immunity from RA_RB to SA_RB
  rate[17] = w ;      // loss of immunity from RA_RB to RA_SB
  
  reulermultinom(2,IA_RB,&rate[16],dt,&trans[16]);
  
  // out of IA_IB
  
  rate[18]=gamma_A;     // recover from IA_IB to CA_IB
  rate[19]=gamma_B;      // recover from IA_IB to IA_CB
  
  reulermultinom(2,IA_IB,&rate[18],dt,&trans[18]);
  
  SA_SB += trans[10] + trans[12]- trans[0] - trans[1];
  IA_SB += trans[0] - trans[2] - trans[3];
  SA_IB += trans[1] - trans[4] - trans[5];
  CA_SB += trans[2] - trans[6] - trans[7];
  SA_CB += trans[4] - trans[8] - trans[9];
  RA_SB += trans[6] + trans[17] - trans[10] - trans[11];
  SA_RB += trans[8] + trans[16] - trans[12] - trans[13];
  IA_IB += trans[3] + trans[5] - trans[18] - trans[19];
  IA_RB += trans[9] + trans[13] + trans[19] - trans[15];
  RA_IB += trans[7] + trans[11] + trans[18] - trans[14];
  RA_RB += trans[14] + trans[15] - trans[16] - trans[17];
  
  K_A += trans[2] + trans[18] + trans [15];           // true incidence of A
  K_B += trans[4] + trans[19] + trans[14];           // true incidence of B

")

# C snippet for initial condition specification ---------------------------


rinit_ee <- Csnippet("
  // how best to initialise population to ensure population size is correct
  // what is this baryometric thing
  // variables iS_A etc refer to endemic equilibrium of single type model
  double iS_A = 1/R0_A;
  double iS_B = 1/R0_B;

  
  double wA_prime = w + phi;
  double wB_prime = w + phi;

  double k_A = gamma_A/wA_prime;
  double k_B = gamma_B/wB_prime;

  double r_A = phi/w;
  double r_B = phi/w;


  double iI_A = (1-iS_A)/(1 + k_A + k_A*r_A);
  double iI_B = (1-iS_B)/(1 + k_B + k_B*r_B);

  double iC_A = k_A*iI_A;
  double iC_B = k_B*iI_B;
  double iR_A = r_A*iC_A;
  double iR_B = r_B*iC_B;

  
  SA_SB = nearbyint(pop*iS_A*iS_B);
  IA_SB = nearbyint(pop*iI_A*iS_B);
  CA_SB = nearbyint(pop*iC_A*iS_B);
  RA_SB = nearbyint(pop*iR_A*iS_B);
  SA_IB = nearbyint(pop*iS_A*iI_B);
  SA_CB = nearbyint(pop*iS_A*iC_B);
  SA_RB = nearbyint(pop*iS_A*iR_B);
  IA_IB = nearbyint(pop*iI_A*iI_B);
  RA_IB = nearbyint(pop*iR_A*iI_B);
  IA_RB = nearbyint(pop*iI_A*iR_B);
  RA_RB = nearbyint(pop*iR_A*iR_B);
  
  W = 0;
  K_A = 0;
  K_B = 0;

")

det_skel <- Csnippet("

  double betaA, betaB, foiA, foiB, seas_A, seas_B, N; //dw foiB
  double rate[20];

  // white noise (extra-demographic stochasticity, 
  // currently asssumed the same for both types(?))
  //dw = rgammawn(sigmaSE,dt);
    
  // sinusoidal seasonality (for both types)
  seas_A = 1 + amplitude_A*cos(2.*3.1415*(t-tpeak_A)/1.); // in units of years
  seas_B = 1 + amplitude_B*cos(2.*3.1415*(t-tpeak_B)/1.); // in units of years

  betaA = R0_A*seas_A*gamma_A;
  betaB = R0_B*seas_B*gamma_B; 


  // forces of infection
  
  N = SA_SB + IA_SB + CA_SB + RA_SB + SA_IB + SA_CB + SA_RB + IA_IB  + RA_IB + IA_RB + RA_RB;

  foiA = betaA*(IA_SB + IA_IB + IA_RB + eta_A)/N; 
  foiB = betaB*(SA_IB + IA_IB + RA_IB + eta_B)/N;
  
    // Out of SA_SB
  rate[0] = SA_SB*foiA; // *dw/dt;        // infection of SA_SB to IA_SB
  rate[1] = SA_SB*foiB; // *dw/dt;        // infection of SA_SB ti SA_IB

  // Out 
  
  rate[2] = IA_SB*gamma_A;                // recovery of IA_SB to CA_SB
  rate[3] = IA_SB*lamda * foiB;         // infection of IA_SB to IA_IB 
  rate[4] = SA_IB*gamma_B;                // recovery of SA_IB to SA_CB
  rate[5] = SA_IB*lamda * foiA;         // infection of SA_IB to IA_IB
  rate[6] = CA_SB*phi;                    //  loss cross-immunity of CA_SB to RA_SB
  rate[7] = CA_SB*chi * foiB;          // infectious of CA_SB to RA_IB
  rate[8] = SA_CB*phi;                    //  loss immunity of SA_CB to SA_RB
  rate[9] = SA_CB*chi * foiA;          // infectious of SA_CB to IA_RB
  rate[10] = RA_SB*w;                  // loss of immunity of RA_SB to SA_SB
  rate[11] = RA_SB*foiB;                 // infectious of RA_SB to RA_IB
  rate[12] = SA_RB*w;                  // loss of immunity to SA_SB
  rate[13] = SA_RB*foiA;                 // infectious of SA_RB to IA_RB
  
  rate[14] = RA_IB*gamma_B;          // recover from RA_IB to RA_RB
  rate[15] = IA_RB*gamma_A;          // recover from IA_RB to RA_RB
  rate[16] = RA_RB*w;       // loss of immunity from RA_RB to SA_RB
  rate[17] = RA_RB*w;       // loss of immunity from RA_RB to RA_SB
  rate[18]= IA_IB*gamma_A;      // recover from IA_IB to CA_IB
  rate[19]= IA_IB*gamma_B;      // recover from IA_IB to IA_CB


  
  // transitions 
  DSA_SB = rate[10] + rate[12]- rate[0] - rate[1];
  DIA_SB = rate[0] - rate[2] - rate[3];
  DSA_IB = rate[1] - rate[4] - rate[5];
  DCA_SB = rate[2] - rate[6] - rate[7];
  DSA_CB = rate[4] - rate[8] - rate[9];
  DRA_SB = rate[6] + rate[17] - rate[10] - rate[11];
  DSA_RB = rate[8] + rate[16] - rate[12] - rate[13];
  DIA_IB = rate[3] + rate[5] - rate[18] - rate[19];
  DIA_RB = rate[9] + rate[13] + rate[19] - rate[15];
  DRA_IB = rate[7] + rate[11] + rate[18] - rate[14];
  DRA_RB = rate[14] + rate[15] - rate[16] - rate[17] ;
  
  DK_A = rate[2] + rate[18] + rate[15] ;           // true incidence of A
  DK_B = rate[4] + rate[19] + rate[14] ;           // true incidence of B
  
  
")
############################################################################
# Code to define estimation components of  model 
############################################################################

# Define likelihood function ----------------------------------------------

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
                   "w","eta_A", "eta_B", "phi",
                   "lamda",
                   "chi", "rho_A", "rho_B", "sigmaSE",
                   "amplitude_A", "tpeak_A", "amplitude_B", "tpeak_B", "pop")
  
  # Model variable names
  model_variables = c("SA_SB","IA_SB", "CA_SB","SA_IB","SA_CB",
                      "SA_RB","IA_IB","RA_IB","IA_RB","RA_RB","RA_SB")
  
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

set.seed(594709947L)
stopifnot(packageVersion("pomp")>="2.0.9.1")
RNGkind("L'Ecuyer-CMRG")


t_start <- -100 ## discard first 100ys



data %>% 
  make_pomp_model(time_start_sim=t_start) -> pomp_sirs

##############################################################################################################
################################ loglikelihood function ##############################################
##############################################################################################################


of<- function(par, pomp_object = pomp_sirs, est, params_all_int = params_all, 
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
    gather(key = type, value = st_cases, -time) %>% 
    # join this with the target data: to calculate poisson log-likelihood
    right_join(target_data, by = c("time", "type")) %>% 
    # calculate poisson log-density: P(data_t|model_t)
    mutate(poisson_log_density = dpois(x = obs_cases, lambda = st_cases, log = TRUE)) -> log_density_data 
  
  # This feature is for potential post-hoc analysis: if conditional log-liklihoods are called for
  if(give_conditional_log_lik == TRUE) {
    
    log_density_data %>% 
      # NA's are given a conditional log-density of 0
      replace_na(list(poisson_log_density = 0)) %>% 
      # Following two steps are used to generate type-specific conditional log-liklihood
      group_by(type) %>% 
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



## set parameters and vairalbes

# Define the date in a long format : To be used in the objective function
data %>%
  gather(key = type, value = obs_cases, -time) -> inc_data

# set a default vector of parameters: requirement of the objective function 

rp_vals <- c(R0_A = 1, gamma_A=365./9, w=1./1.,
                  R0_B = 1, gamma_B=365./3,
                  phi=365/30,
                  chi=1, eta_A=365., eta_B=365.,
                  lamda=1,
                  rho_A=0, rho_B=0, 
                  sigmaSE=0.0000, 
                  amplitude_A=0, tpeak_A=0, amplitude_B=0, tpeak_B=0,
                  pop=Npop)

SIRS2_independent_endemic_equilibrium <- function(params){
  S_A <- 1/params[["R0_A"]]
  S_B <- 1/params[["R0_B"]]
  
  k_A <- params[["gamma_A"]]/(params[["w"]] + params[["phi"]])
  k_B <- params[["gamma_B"]]/(params[["w"]] + params[["phi"]])
  
  r_A <-  params[["phi"]]/params[["w"]]
  r_B <-  params[["phi"]]/params[["w"]]
  
  I_A <- (1-S_A)/(1+k_A + k_A*r_A)
  I_B <- (1-S_B)/(1+k_B + k_B*r_B)
  
  ee_A <- c(S_0=S_A, I_0=I_A,  C_0=k_A*I_A, R_0=r_A*k_A*I_A)
  ee_B <- c(S_0=S_B, I_0=I_B,  C_0=k_B*I_B, R_0=r_B*k_B*I_B)
  
  # Note only  classes: model assumes coinfection (I_A * I_B) is neglible
  return(c(SA_SB_0 = ee_A[["S_0"]]*ee_B[["S_0"]],
           IA_SB_0 = ee_A[["I_0"]]*ee_B[["S_0"]],
           CA_SB_0 = ee_A[["C_0"]]*ee_B[["S_0"]],
           RA_SB_0 = ee_A[["R_0"]]*ee_B[["S_0"]],
           SA_IB_0 = ee_A[["S_0"]]*ee_B[["I_0"]],
           SA_CB_0 = ee_A[["S_0"]]*ee_B[["C_0"]],
           SA_RB_0 = ee_A[["S_0"]]*ee_B[["R_0"]],
           IA_IB_0 = ee_A[["I_0"]]*ee_B[["I_0"]],
           RA_IB_0 = ee_A[["R_0"]]*ee_B[["I_0"]],
           IA_RB_0 = ee_A[["I_0"]]*ee_B[["R_0"]],
           RA_RB_0 = ee_A[["R_0"]]*ee_B[["R_0"]]))
}

ic_vals <- SIRS2_independent_endemic_equilibrium(rp_vals)
params_all <- c(rp_vals,ic_vals)
##################################

# create a sampled design of initial guesses: Sobol sampling
## this is for full parameter
# pd <- sobol_design(lower = c(R0_B = 1, R0_A = 1, amplitude_A = 1e-10, 
#                              amplitude_B = 1e-10, tpeak_A = 1/365.25, tpeak_B = 1/365.25,  
#                              rho_A = 1e-10, rho_B = 1e-10, chi =0, lamda=0), 
#                    upper = c(R0_B = 5, R0_A = 5, amplitude_A = 1, amplitude_B = 1, 
#                              tpeak_A = 1, tpeak_B = 1, rho_A = 0.01, rho_B = 0.01, chi =1, lamda=1), 
#                    nseq = 499)
# 
#write.csv(pd,"pd.csv",row.names = FALSE )

#####################################################################################
## parameter setting sepecific
#####################################################################################
##


# This function is defined to provide pomp-GA integration
## null model
# parameters of interest
est_these_null <- c("R0_B", "R0_A", "amplitude_A", "amplitude_B", "tpeak_A", "tpeak_B", 
                    "rho_A", "rho_B")

## from loglikelihood function
of_ga_null<- function(x, est = est_these_null) {
  
  x<- unname(unlist(x))
  
  of(par = c(R0_B = x[1], R0_A = x[2], 
             amplitude_A = x[3], amplitude_B = x[4], 
             tpeak_A = x[5], tpeak_B = x[6],  
             rho_A = x[7], rho_B = x[8]), est = est)
}

pd %>% 
  select(-lamda,-chi)-> pd_withfix_null

pd %>% 
  select(-lamda,-chi)%>% 
  full_join(prior)-> pd_withfix_null
################################################################################3
##
## numerical optimizer
##################################################################################


# genetic algorithm set up to run in parallel
tic()
no_cores <- detectCores() - 1  

registerDoParallel(cores=no_cores)  

cl <- makeCluster(no_cores, type="FORK")

GA_null<- ga(type = "real-valued", 
             fitness = of_ga_null,
             lower = c(R0_B = 1, R0_A = 1, amplitude_A = 1e-10, 
                       amplitude_B = 1e-10, tpeak_A = 1/365.25, tpeak_B = 1/365.25,  
                       rho_A = 1e-10, rho_B = 1e-10), 
             upper = c(R0_B = 5, R0_A = 5, 
                       amplitude_A = 1, amplitude_B = 1, 
                       tpeak_A = 1, tpeak_B = 1, 
                       rho_A = 0.01, rho_B = 0.01),
             elitism = base::max(1, round(100*.1)), 
             popSize = 500, maxiter = 10000, run = 499, 
             suggestions = pd_withfix_null,
             optim = TRUE,
             optimArgs = list(poptim = 0.1, 
                              control = list(maxit = 1000)),
             seed = 123456, parallel = cl)

##########################################################################################################
## set with DEoptim
library(DEoptim)
## no need sobol design
outDEoptim <- DEoptim(of_ga_null, ## from pomp_ga integration
                      lower = c(R0_B = 1, R0_A = 1, amplitude_A = 1e-10, 
                                amplitude_B = 1e-10, tpeak_A = 1/365.25, tpeak_B = 1/365.25,  
                                rho_A = 1e-10, rho_B = 1e-10), 
                      upper = c(R0_B = 5, R0_A = 5, 
                                amplitude_A = 1, amplitude_B = 1, 
                                tpeak_A = 1, tpeak_B = 1, 
                                rho_A = 0.01, rho_B = 0.01),
                      control = DEoptim.control(NP = 499, itermax=10000, F=0.8, CR=0.5,p=0.2,
                                                reltol=1e-3,
                                                parallelType=1,cluster=cl ))

#######################################################################################################

stopCluster(cl)
toc() -> ga_took


result_null <- list(it_took  = c(paste((ga_took$toc - ga_took$tic)/3600, "hrs"), 
                                          paste((ga_took$toc - ga_took$tic)/60, "mins"), 
                                          paste((ga_took$toc - ga_took$tic), "secs")), 
                             GAobj = GA_null)


save(result_null, file = "result_null.Rdata")
#######################################################################################3
#   
# #####################################################################################
# ## parameter setting sepecific
# #####################################################################################
# ##
# 
# 
# # This function is defined to provide pomp-GA integration
# ## lamda model
# # parameters of interest
# est_these_lamda <- c("R0_B", "R0_A", "amplitude_A", "amplitude_B", "tpeak_A", "tpeak_B", 
#                      "rho_A", "rho_B","lamda")
# 
# of_ga_lamda<- function(x, est = est_these_lamda) {
#   
#   x<- unname(unlist(x))
#   
#   of(par = c(R0_B = x[1], R0_A = x[2], 
#              amplitude_A = x[3], amplitude_B = x[4], 
#              tpeak_A = x[5], tpeak_B = x[6],  
#              rho_A = x[7], rho_B = x[8], lamda =x[9]), est = est)
# }
# 
# 
# 
# 
# pd %>% 
#   select(-chi)%>% 
#   full_join(prior)%>%
#   replace_na(replace = list(lamda=1))-> pd_withfix_lamda
# 
# # genetic algorithm set up to run in parallel
# tic()
# no_cores <- detectCores() - 1  
# 
# registerDoParallel(cores=no_cores)  
# 
# cl <- makeCluster(no_cores, type="FORK")
# 
# GA_lamda<- ga(type = "real-valued", 
#              fitness = of_ga_lamda,
#              lower = c(R0_B = 1, R0_A = 1, amplitude_A = 1e-10, 
#                        amplitude_B = 1e-10, tpeak_A = 1/365.25, tpeak_B = 1/365.25,  
#                        rho_A = 1e-10, rho_B = 1e-10, lamda=0), 
#              upper = c(R0_B = 5, R0_A = 5, 
#                        amplitude_A = 1, amplitude_B = 1, 
#                        tpeak_A = 1, tpeak_B = 1, 
#                        rho_A = 0.01, rho_B = 0.01, lamda=1),
#              elitism = base::max(1, round(100*.1)), 
#              popSize = 500, maxiter = 10000, run = 500, 
#              suggestions = pd_withfix_lamda,
#              optim = TRUE,
#              optimArgs = list(poptim = 0.1, 
#                               control = list(maxit = 1000)),
#              seed = 123456, parallel = cl)
# 
# stopCluster(cl)
# toc() -> ga_took
# 
# 
# result_lamda <- list(it_took  = c(paste((ga_took$toc - ga_took$tic)/3600, "hrs"), 
#                                           paste((ga_took$toc - ga_took$tic)/60, "mins"), 
#                                           paste((ga_took$toc - ga_took$tic), "secs")), 
#                              GAobj = GA_lamda)
# 
# 
# save(result_lamda, file = "result_lamda.Rdata")
# #######################################################################################3
# 
# 
# #####################################################################################
# ## parameter setting sepecific
# #####################################################################################
# ##
# 
# 
# # This function is defined to provide pomp-GA integration
# ## lamda model
# # parameters of interest
# est_these_chi <- c("R0_B", "R0_A", "amplitude_A", "amplitude_B", "tpeak_A", "tpeak_B", 
#                      "rho_A", "rho_B","chi")
# 
# of_ga_chi<- function(x, est = est_these_chi) {
#   
#   x<- unname(unlist(x))
#   
#   of(par = c(R0_B = x[1], R0_A = x[2], 
#              amplitude_A = x[3], amplitude_B = x[4], 
#              tpeak_A = x[5], tpeak_B = x[6],  
#              rho_A = x[7], rho_B = x[8], chi =x[9]), est = est)
# }
# 
# 
# 
# 
# pd %>% 
#   select(-lamda)%>% 
#   full_join(prior)%>%
#   replace_na(replace = list(chi=1))-> pd_withfix_chi
# 
# # genetic algorithm set up to run in parallel
# tic()
# no_cores <- detectCores() - 1  
# 
# registerDoParallel(cores=no_cores)  
# 
# cl <- makeCluster(no_cores, type="FORK")
# 
# GA_chi<- ga(type = "real-valued", 
#               fitness = of_ga_chi,
#               lower = c(R0_B = 1, R0_A = 1, amplitude_A = 1e-10, 
#                         amplitude_B = 1e-10, tpeak_A = 1/365.25, tpeak_B = 1/365.25,  
#                         rho_A = 1e-10, rho_B = 1e-10, chi=0), 
#               upper = c(R0_B = 5, R0_A = 5, 
#                         amplitude_A = 1, amplitude_B = 1, 
#                         tpeak_A = 1, tpeak_B = 1, 
#                         rho_A = 0.01, rho_B = 0.01, chi=1),
#               elitism = base::max(1, round(100*.1)), 
#               popSize = 500, maxiter = 10000, run = 500, 
#               suggestions = pd_withfix_chi,
#               optim = TRUE,
#               optimArgs = list(poptim = 0.1, 
#                                control = list(maxit = 1000)),
#               seed = 123456, parallel = cl)
# 
# stopCluster(cl)
# toc() -> ga_took
# 
# 
# result_chi<- list(it_took  = c(paste((ga_took$toc - ga_took$tic)/3600, "hrs"), 
#                                   paste((ga_took$toc - ga_took$tic)/60, "mins"), 
#                                   paste((ga_took$toc - ga_took$tic), "secs")), 
#                      GAobj = GA_chi)
# 
# 
# save(result_chi, file = "result_chi.Rdata")
# #######################################################################################3
# 
# 
# #####################################################################################
# ## parameter setting sepecific
# #####################################################################################
# ##
# 
# 
# # This function is defined to provide pomp-GA integration
# ## lamda model
# # parameters of interest
# est_these_fixphi <- c("R0_B", "R0_A", "amplitude_A", "amplitude_B", "tpeak_A", "tpeak_B", 
#                    "rho_A", "rho_B","lamda", "chi")
# 
# of_ga_fixphi<- function(x, est = est_these_fixphi) {
#   
#   x<- unname(unlist(x))
#   
#   of(par = c(R0_B = x[1], R0_A = x[2], 
#              amplitude_A = x[3], amplitude_B = x[4], 
#              tpeak_A = x[5], tpeak_B = x[6],  
#              rho_A = x[7], rho_B = x[8], lamda =x[9], chi = x[10]), est = est)
# }
# 
# 
# 
# 
# pd %>% 
#   full_join(prior)%>%
#   replace_na(replace = list(lamda=1, chi=1))-> pd_withfix_fixphi
# 
# # genetic algorithm set up to run in parallel
# tic()
# no_cores <- detectCores() - 1  
# 
# registerDoParallel(cores=no_cores)  
# 
# cl <- makeCluster(no_cores, type="FORK")
# 
# GA_fixphi<- ga(type = "real-valued", 
#             fitness = of_ga_fixphi,
#             lower = c(R0_B = 1, R0_A = 1, amplitude_A = 1e-10, 
#                       amplitude_B = 1e-10, tpeak_A = 1/365.25, tpeak_B = 1/365.25,  
#                       rho_A = 1e-10, rho_B = 1e-10, lamda=0, chi=0), 
#             upper = c(R0_B = 5, R0_A = 5, 
#                       amplitude_A = 1, amplitude_B = 1, 
#                       tpeak_A = 1, tpeak_B = 1, 
#                       rho_A = 0.01, rho_B = 0.01, lamda=1, chi=1),
#             elitism = base::max(1, round(100*.1)), 
#             popSize = 500, maxiter = 10000, run = 500, 
#             suggestions = pd_withfix_fixphi,
#             optim = TRUE,
#             optimArgs = list(poptim = 0.1, 
#                              control = list(maxit = 1000)),
#             seed = 123456, parallel = cl)
# 
# stopCluster(cl)
# toc() -> ga_took
# 
# 
# result_fixphi<- list(it_took  = c(paste((ga_took$toc - ga_took$tic)/3600, "hrs"), 
#                                paste((ga_took$toc - ga_took$tic)/60, "mins"), 
#                                paste((ga_took$toc - ga_took$tic), "secs")), 
#                   GAobj = GA_fixphi)
# 
# 
# save(result_fixphi, file = "result_fixphi.Rdata")
# #######################################################################################3


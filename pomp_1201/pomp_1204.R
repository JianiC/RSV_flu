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


pop<-read.csv("pomp/pop_HHSRegion.csv")
Npop=pop[pop$HHS_region == HHS_n,3] ##population param
case_data<-read.csv("pomp/pomp_HHSregion_case.csv") ## case count data
data<-case_data%>%
  filter(HHS_REGION==HHS_n)%>%
  select(time,RSVpos,fluApos)%>%
  dplyr::rename(total1=RSVpos,total2=fluApos)%>%
  add_row(time = 2.5, total1 = NA, total2=NA, .before = 1)

pd<-read.csv("pd.csv") ## sobol design
#############################################################################

## POMP model

#############################################################################



rproc_euler_multinomial <- Csnippet("
  
  double beta1, beta2, foi1, foi2, seas1, seas2, N; 
  
  double rate[49], trans[49];

  // sinusoidal seasonality (for both types)
  seas1 = 1 + amplitude1*cos(2.*3.1415*(t-tpeak1)/1.); // in units of years
  seas2 = 1 + amplitude2*cos(2.*3.1415*(t-tpeak2)/1.); // in units of years

  beta1 = R01*seas1*gamma1;
  beta2 = R02*seas2*gamma2; 


  // forces of infection
  N = S1_S2 + I1_S2 + C1_S2 + R1_S2 + S1_I2 + I1_I2 + C1_I2 + R1_I2 + S1_C2 + I1_C2 + C1_C2 + R1_C2 + S1_R2 + I1_R2 + C1_R2 + R1_R2;
  
  //I1_sum = I1_S2 + I1_I2 + I1_C2 + I1_R2;
  //I2_sum = S1_I2 + I1_I2 + C1_I2 + R1_I2;

  foi1 = beta1 * (I1_S2 + I1_I2 + I1_C2 + I1_R2 + eta1)/N; 
  foi2 = beta2 * (S1_I2 + I1_I2 + C1_I2 + R1_I2 + eta2)/N;
  
  //birth (assume birth rate equal to death rate)
 // rate[0] = births;
  rate[0] = mu;
  reulermultinom(1,N,&rate[0],dt,&trans[0]);
  
 // trans[0] = rpois(rate[0] * dt);
  
  // Out of S1_S2
  rate[1] = foi1;         // infection of S1_S1 -> I1_S2
  rate[2] = foi2;         // infection of S1_S2 -> S1_I2
  rate[3] = mu;            // death from S1_S2   
                
  reulermultinom(3,S1_S2,&rate[1],dt,&trans[1]);
  
  // Out of I1_S2
  rate[4] = gamma1;                // recovery of I1_S2 -> C1_S2
  rate[5] = psi * foi2;         // infection of I1_S2 -> I1_I2 
  rate[6] = mu;            // death from I1_S2  
  
  reulermultinom(3,I1_S2,&rate[4],dt,&trans[4]);
  
    // Out of C1_S2
  rate[7] = phi1;              // recovery of C1_S2 -> R1_S2
  rate[8] = chi * foi2;         // infection of C1_S2 -> C1_I2
  rate[9] = mu;            // death from C1_S2  
  
  reulermultinom(3,C1_S2,&rate[7],dt,&trans[7]);
  
  // Out of R1_S2
  rate[10] = w1;                // recovery of R1_S2->S1_S2
  rate[11] = foi2;         // infection of R1_S2->R1_I2
  rate[12] = mu;            // death from R1_S2  
  
  reulermultinom(3,R1_S2,&rate[10],dt,&trans[10]);

  // Out of S1_I2
  rate[13] = psi * foi1;         // infection of S1_I2 -> I1_I2
  rate[14] = gamma2;         // infection of S1_I2 -> S1_C2
  rate[15] = mu;            // death from S1_I2  
  
  reulermultinom(3,S1_I2,&rate[13],dt,&trans[13]);
  
  // Out of I1_I2
  rate[16] = gamma1;                // recovery of I1_I2 -> C1_I2
  rate[17] = gamma2;         // infection of I1_I2 -> I1_C2 
  rate[18] = mu;            // death from I1_I2  
  
  reulermultinom(3,I1_I2,&rate[16],dt,&trans[16]);
  
  
    // Out of C1_I2
  rate[19] = phi1;              // loss cross immunity from C1_I2->R1_I2
  rate[20] = gamma2;            // recoveru of C1_I2-> C1_C2
  rate[21] = mu;            // death from C1_I2  
  
  reulermultinom(3,C1_I2,&rate[19],dt,&trans[19]);


  // Out of R1_I2
  
  rate[22] = w1;                // loss immunity R1_I2->S1_I2
  rate[23] = gamma2;         // recovery of R1_I2->R1_R2
  rate[24] = mu;            // death from R1_I2  
  
  reulermultinom(3,R1_I2,&rate[22],dt,&trans[22]);


  // Out of S1_C2
  rate[25] = chi* foi1;         // infection of S1_C2 -> I1_C2
  rate[26] = phi2;         // loss cross immunity S1_C2->S1_R2
  rate[27] = mu;            // death from S1_C2  
  
  reulermultinom(3,S1_C2,&rate[25],dt,&trans[25]);
  
  // Out of I1_C2
  rate[28] = gamma1;          // recovery of I1_C2 -> C1_C2
  rate[29] = phi2;         // loss cross immunity I1_C2 -> I1_R2 
  rate[30] = mu;            // death from I1_C2  
  
  reulermultinom(3,I1_C2,&rate[28],dt,&trans[28]);
  

    // Out of C1_C2
  rate[31] = phi1;              // loss cross immunity from C1_C2->R1_C2
  rate[32] = phi2;            // loss cross immunity C1_C2-> C1_R2
  rate[33] = mu;            // death from C1_C2  
  
  reulermultinom(3,C1_C2,&rate[31],dt,&trans[31]);


  // Out of R1_C2
  
  rate[34] = w1;                // loss immunity R1_C2->S1_C2
  rate[35] = phi2;         // ross cross immunity R1_C2-> R1_R2
  rate[36] = mu;            // death from R1_C2  
  
  reulermultinom(3,R1_C2,&rate[34],dt,&trans[34]);
  
  
  // Out of S1_R2
  rate[37] = foi1;         // infection of S1_R2 -> I1_R2
  rate[38] =  w2;         // loss  immunity S1_R2->S1_S2
  rate[39] = mu;            // death from S1_R2  
  
  reulermultinom(3,S1_R2,&rate[37],dt,&trans[37]);
  
  // Out of I1_R2
  rate[40] = gamma1;          // recovery of I1_R2 -> C1_R2
  rate[41] = w2;         // loss immunity I1_R2 -> I1_S2 
  rate[42] = mu;            // death from I1_R2  
  
  reulermultinom(3,I1_R2,&rate[40],dt,&trans[40]);
  

   // Out of C1_R2
  rate[43] = phi1;              // loss cross immunity from C1_R2->R1_R2
  rate[44] = w2;            // loss immunity C1_R2-> C1_S2
  rate[45] = mu;            // death from C1_R2  
  
  reulermultinom(3,C1_R2,&rate[43],dt,&trans[43]);


  // Out of R1_R2
  
  rate[46] = w1;                // loss immunity R1_R2->S1_R2
  rate[47] = w2;         // lossimmunity R1_R2-> R1_S2
  rate[48] = mu;            // death from S1_S2  
  
  reulermultinom(3,R1_R2,&rate[46],dt,&trans[46]);
    
  
  S1_S2 += trans[0] + trans[10] + trans[38] - trans[1] - trans[2] - trans[3];
  I1_S2 += trans[1] + trans[41] - trans[4] - trans[5] - trans[6];
  C1_S2 += trans[4] + trans[44] - trans[7] - trans[8] - trans[9];
  R1_S2 += trans[7] + trans[47] - trans[10] - trans[11] - trans[12];
  
  S1_I2 += trans[22] + trans[2] - trans[13] - trans[14] - trans[15];
  I1_I2 += trans[13] + trans[5] - trans[16] - trans[17] - trans[18];
  C1_I2 += trans[16] + trans[8] - trans[19] - trans[20] - trans[21];
  R1_I2 += trans[19] + trans[11] - trans[22] - trans[23] - trans[24];
  
  S1_C2 += trans[34] + trans[14] - trans[25] - trans[26] - trans[27];
  I1_C2 += trans[25] + trans[17] - trans[28] - trans[29] - trans[30];
  C1_C2 += trans[28] + trans[20] - trans[31] - trans[32] - trans[33];
  R1_C2 += trans[31] + trans[23] - trans[34] - trans[35] - trans[36];
  
  S1_R2 += trans[46] + trans[26] - trans[37] - trans[38] - trans[39];
  I1_R2 += trans[37] + trans[29] - trans[40] - trans[41] - trans[42];
  C1_R2 += trans[40] + trans[32] - trans[43] - trans[44] - trans[45];
  R1_R2 += trans[43] + trans[35] - trans[46] - trans[47] - trans[48];
  

  Kprim1 += trans[4];      // from I1_S2->C1_S2
  Kprim2 += trans[14];     // from S1_I2 ->S1_C2
  
  Ksec1 += trans[16] + trans[28] + trans[40]; // from I1_I2, I1_C2, I1_R2
  Ksec2 += trans[17] + trans[20] + trans[23]; // from I1_I2, C1_I2, R1_I2
  
  K1 += Kprim1 + Ksec1;           // true incidence of pathogen1
  K2 += Kprim2 + Ksec2;           // true incidence of pathogen2
  

  

")

# C snippet for initial condition specification ---------------------------


rinit_ee <- Csnippet("
  // how best to initialise population to ensure population size is correct
  // what is this baryometric thing
  // variables iS_A etc refer to endemic equilibrium of single type model
  
  double iS1 = 1/R01;
  double iS2 = 1/R02;

  
  double w1_prime = w1 + phi1;
  double w2_prime = w2 + phi2;

  double k1 = gamma1/w1_prime;
  double k2 = gamma2/w2_prime;

  double r1 = phi1/w1;
  double r2 = phi2/w2;

  double iI1 = (1-iS1)/(1 + k1 + k1*r1);
  double iI2 = (1-iS2)/(1 + k2 + k2*r2);

  double iC1 = k1*iI1;
  double iC2 = k2*iI2;
  
  double iR1 = r1*iC1;
  double iR2 = r2*iC2;

  
  S1_S2 = nearbyint(pop*iS1*iS2);
  I1_S2 = nearbyint(pop*iI1*iS2);
  C1_S2 = nearbyint(pop*iC1*iS2);
  R1_S2 = nearbyint(pop*iR1*iS2);
  
  S1_I2 = nearbyint(pop*iS1*iI2);
  I1_I2 = nearbyint(pop*iI1*iI2);
  C1_I2 = nearbyint(pop*iC1*iI2);
  R1_I2 = nearbyint(pop*iR1*iI2);
  
  S1_C2 = nearbyint(pop*iS1*iC2);
  I1_C2 = nearbyint(pop*iI1*iC2);
  C1_C2 = nearbyint(pop*iC1*iC2);
  R1_C2 = nearbyint(pop*iR1*iC2);
  
  S1_R2 = nearbyint(pop*iS1*iR2);
  I1_R2 = nearbyint(pop*iI1*iR2);
  C1_R2 = nearbyint(pop*iC1*iR2);
  R1_R2 = nearbyint(pop*iR1*iR2);
  

  K1 = 0;
  K2 = 0;
  
  
")

det_skel <- Csnippet("

  double beta1, beta2, foi1, foi2, seas1, seas2, N; 
  double rate[49];

    
 // sinusoidal seasonality (for both types)
  seas1 = 1 + amplitude1*cos(2.*3.1415*(t-tpeak1)/1.); // in units of years
  seas2 = 1 + amplitude2*cos(2.*3.1415*(t-tpeak2)/1.); // in units of years

  beta1 = R01*seas1*gamma1;
  beta2 = R02*seas2*gamma2; 


  // forces of infection
  
  N = S1_S2 + I1_S2 + C1_S2 + R1_S2 + S1_I2 + I1_I2 + C1_I2 + R1_I2 + S1_C2 + I1_C2 + C1_C2 + R1_C2 + S1_R2 + I1_R2 + C1_R2 + R1_R2;

  //I1_sum = I1_S2 + I1_I2 + I1_C2 + I1_R2;
  //I2_sum = S1_I2 + I1_I2 + C1_I2 + R1_I2;

  foi1 = beta1 * (I1_S2 + I1_I2 + I1_C2 + I1_R2 + eta1)/N; 
  foi2 = beta2 * (S1_I2 + I1_I2 + C1_I2 + R1_I2 + eta2)/N;
  
  
  //birth
  rate[0] = N * mu;
  
  // Out of S1_S2
  rate[1] = S1_S2 * foi1;         // infection of S1_S1 -> I1_S2
  rate[2] = S1_S2 * foi2;         // infection of S1_S2 -> S1_I2
  rate[3] = S1_S2 * mu;            // death from S1_S2  

  // Out of I1_S2
  rate[4] = I1_S2 * gamma1;                // recovery of I1_S2 -> C1_S2
  rate[5] = I1_S2 * psi * foi2;         // infection of IA_SB -> IA_IB 
  rate[6] = I1_S2 * mu;            // death from I1_S2  
  
  // Out of C1_S2
  rate[7] = C1_S2 * phi1;              // recovery of C1_S2 -> R1_S2
  rate[8] = C1_S2 * chi * foi2;         // infection of C1_S2 -> C1_I2
  rate[9] = C1_S2 * mu;            // death from C1_S2 
  
  // Out of R1_S2
  rate[10] = R1_S2 * w1;                // recovery of R1_S2->S1_S2
  rate[11] = R1_S2 * foi2;         // infection of R1_S2->R1_I2
  rate[12] = R1_S2 * mu;            // death from R1_S2 
  
   // Out of S1_I2
  rate[13] = S1_I2 * psi * foi1;         // infection of S1_I2 -> I1_I2
  rate[14] = S1_I2 * gamma2;         // infection of S1_I2 -> S1_C2
  rate[15] = S1_I2 * mu;            // death from S1_I2  
  
    // Out of I1_I2
  rate[16] = I1_I2 * gamma1;                // recovery of I1_I2 -> C1_I2
  rate[17] = I1_I2 * gamma2;         // infection of I1_I2 -> I1_C2
  rate[18] = I1_I2 * mu;            // death from I1_I2 
  
  // Out of C1_I2
  rate[19] = C1_I2 * phi1;              // loss cross immunity from C1_I2->R1_I2
  rate[20] = C1_I2 * gamma2;            // recoveru of C1_I2-> C1_C2
  rate[21] = C1_I2 * mu;            // death from C1_I2 
  
    // Out of R1_I2
  rate[22] = R1_I2 * w1;                // loss immunity R1_I2->S1_I2
  rate[23] = R1_I2 * gamma2;         // recovery of R1_I2->R1_R2
  rate[24] = R1_I2 * mu;            // death from R1_I2
  
  // Out of S1_C2
  rate[25] = S1_C2 * chi* foi1;         // infection of S1_C2 -> I1_C2
  rate[26] = S1_C2 * phi2;         // loss cross immunity S1_C2->S1_R2
  rate[27] = S1_C2 * mu;            // death from S1_C2 
  
  // Out of I1_C2
  rate[28] = I1_C2 * gamma1;          // recovery of I1_C2 -> C1_C2
  rate[29] = I1_C2 * phi2;         // loss cross immunity I1_C2 -> I1_R2
  rate[30] = I1_C2 * mu;            // death from I1_C2 
  
  // Out of C1_C2
  rate[31] = C1_C2 * phi1;              // loss cross immunity from C1_C2->R1_C2
  rate[32] = C1_C2 * phi2;            // loss cross immunity C1_C2-> C1_R2
  rate[33] = C1_C2 * mu;            // death from C1_C2  

  // Out of R1_C2
  rate[34] = R1_C2 * w1;                // loss immunity R1_C2->S1_C2
  rate[35] = R1_C2 * phi2;         // loss cross immunity R1_C2-> R1_R2
  rate[36] = R1_C2 * mu;            // death from R1_C2  
  
  // Out of S1_R2
  rate[37] = S1_R2 * foi1;         // infection of S1_R2 -> I1_R2
  rate[38] = S1_R2 * w2;         // loss  immunity S1_R2->S1_S2
  rate[39] = S1_R2 * mu;            // death from S1_R2 
  
  // Out of I1_R2
  rate[40] = I1_R2 * gamma1;          // recovery of I1_R2 -> C1_R2
  rate[41] = I1_R2 * w2;         // loss immunity I1_R2 -> I1_S2
  rate[42] = I1_R2 * mu;            // death from I1_R2 
  
    // Out of C1_R2
  rate[43] = C1_R2 * phi1;              // loss cross immunity from C1_R2->R1_R2
  rate[44] = C1_R2 * w2;            // loss immunity C1_R2-> C1_S2
  rate[45] = C1_R2 * mu;            // death from C1_R2
  
    // Out of R1_R2
  rate[46] = R1_R2 * w1;                // loss immunity R1_R2->S1_R2
  rate[47] = R1_R2 * w2;         // lossimmunity R1_R2-> R1_S2
  rate[48] = R1_R2 * mu;            // death from S1_S2 
  
  
  // transitions 
  
  DS1_S2 = rate[0] + rate[10] + rate[38] - rate[1] - rate[2] - rate[3];
  DI1_S2 = rate[1] + rate[41] - rate[4] - rate[5] - rate[6];
  DC1_S2 = rate[4] + rate[44] - rate[7] - rate[8] - rate[9];
  DR1_S2 = rate[7] + rate[47] - rate[10] - rate[11] - rate[12];
  
  DS1_I2 = rate[22] + rate[2] - rate[13] - rate[14] - rate[15];
  DI1_I2 = rate[13] + rate[5] - rate[16] - rate[17] - rate[18];
  DC1_I2 = rate[16] + rate[8] - rate[19] - rate[20] - rate[21];
  DR1_I2 = rate[19] + rate[11] - rate[22] - rate[23] - rate[24];
  
  DS1_C2 = rate[34] + rate[14] - rate[25] - rate[26] - rate[27];
  DI1_C2 = rate[25] + rate[17] - rate[28] - rate[29] - rate[30];
  DC1_C2 = rate[28] + rate[20] - rate[31] - rate[32] - rate[33];
  DR1_C2 = rate[31] + rate[23] - rate[34] - rate[35] - rate[36];
  
  DS1_R2 = rate[46] + rate[26] - rate[37] - rate[38] - rate[39];
  DI1_R2 = rate[37] + rate[29] - rate[40] - rate[41] - rate[42];
  DC1_R2 = rate[40] + rate[32] - rate[43] - rate[44] - rate[45];
  DR1_R2 = rate[43] + rate[35] - rate[46] - rate[47] - rate[48];
  

  DKprim1 = rate[4];      // from I1_S2->C1_S2
  DKprim2 = rate[14];     // from S1_I2 ->S1_C2
  
  DKsec1 = rate[16] + rate[28] + rate[40]; // from I1_I2, I1_C2, I1_R2
  DKsec2 = rate[17] + rate[20] + rate[23]; // from I1_I2, C1_I2, R1_I2
  
  DK1 = DKprim1 + DKsec1;           // true incidence of pathogen1
  DK2 = DKprim2 + DKsec2;           // true incidence of pathogen2

  
")
############################################################################
# Code to define estimation components of  model 
############################################################################

# Define likelihood function ----------------------------------------------

dmeas_poisson <- Csnippet("
  double lik1, lik2;
  
  if(ISNA(total1)) {
    lik1= (give_log) ? 0:1;
  } else {
      lik1 = dpois(nearbyint(total1), rho1*K1, give_log); 
  }
  
  if(ISNA(total2)) {
    lik2 = (give_log) ? 0:1;
  } else {
      lik2 = dpois(nearbyint(total2), rho2*K2, give_log); 
  }
  
  lik = lik1 + lik2;

")



# Define process simulator for observations  ------------------------------
rmeas_poisson <- Csnippet("

  total1 = rpois(rho1*K1);
  total2 = rpois(rho2*K2);
  
")
############################################################################
# assign values for model parameters and initial conditions -----------------
############################################################################

# Define the pomp model object --------------------------------------------

make_pomp_model <- function(df, time_start_sim, dt=1) {
  
  
  # Model parameter names
  model_params = c("R01","R02", "gamma1", "gamma2",
                   "w1","w2","eta1", "eta2", "phi1", "phi2",
                   "psi","chi", "rho1", "rho2",
                   "amplitude1", "tpeak1", "amplitude2", "tpeak2", "pop","mu")
  
  # Model variable names
  model_variables = c("S1_S2","I1_S2", "C1_S2","R1_S2",
                      "S1_I2","I1_I2", "C1_I2","R1_I2",
                      "S1_C2","I1_C2", "C1_C2","R1_C2",
                      "S1_R2","I1_R2", "C1_R2","R1_R2")
  
  # Model initial conditions parameter names
  model_ic_params = paste(model_variables, "_0", sep="")
  
  process_model <- euler(rproc_euler_multinomial, delta.t = dt/365.25)
  
  dmeas <- dmeas_poisson
  rmeas <- rmeas_poisson
  
  po <- pomp(
    data = df,
    times = "time",
    t0 = time_start_sim,
    obsnames = c("total1","total2"),
    rprocess = process_model,
    skeleton = vectorfield(det_skel),
    dmeasure = dmeas,
    rmeasure = rmeas, 
    rinit = rinit_ee,
    accumvars =c ("Kprim1","Ksec1","Kprim2","Ksec2","K1", "K2" ),
    statenames = c(model_variables, c("Kprim1","Ksec1","Kprim2","Ksec2","K1", "K2")),
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
################################ set paramerters and variables ##############################################
##############################################################################################################
## set parameters and vairalbes

# Define the date in a long format : To be used in the objective function
data %>%
  gather(key = pathogen, value = obs_cases, -time) -> inc_data

# set a default vector of parameters: requirement of the objective function 

rp_vals <- c(R01 = 1, gamma1=365./9, w1=1./1.,
             R02 = 1, gamma2=365./3, w2=1./1.,
             phi1=365/30, phi2=365/30, psi =1, chi=1, 
             eta1=365., eta2=365.,rho1=0, rho2=0, 
             amplitude1=0, tpeak1=0, amplitude2=0, tpeak2=0,pop=Npop, mu=0)

SIRS2_independent_endemic_equilibrium <- function(params){
  S1 <- 1/params[["R01"]]
  S2 <- 1/params[["R02"]]
  
  k1 <- params[["gamma1"]]/(params[["w1"]] + params[["phi1"]])
  k2 <- params[["gamma2"]]/(params[["w2"]] + params[["phi2"]])
  
  r1 <-  params[["phi1"]]/params[["w1"]]
  r2 <-  params[["phi2"]]/params[["w2"]]
  
  I1 <- (1-S1)/(1+k1 + k1*r2)
  I2 <- (1-S2)/(1+k2 + k2*r2)
  
  ee1 <- c(S0=S1, I0=I1,  C0=k1*I1, R0=r1*k1*I1)
  ee2 <- c(S0=S2, I0=I2,  C0=k2*I2, R0=r2*k2*I2)
  
  return(c(S1_S2_0 = ee1[["S0"]]*ee2[["S0"]],
           I1_S2_0 = ee1[["I0"]]*ee2[["S0"]],
           C1_S2_0 = ee1[["C0"]]*ee2[["S0"]],
           R1_S2_0 = ee1[["R0"]]*ee2[["S0"]],
           
           S1_I2_0 = ee1[["S0"]]*ee2[["I0"]],
           I1_I2_0 = ee1[["I0"]]*ee2[["I0"]],
           C1_I2_0 = ee1[["C0"]]*ee2[["I0"]],
           R1_I2_0 = ee1[["R0"]]*ee2[["I0"]],
           
           S1_C2_0 = ee1[["S0"]]*ee2[["C0"]],
           I1_C2_0 = ee1[["I0"]]*ee2[["C0"]],
           C1_C2_0 = ee1[["C0"]]*ee2[["C0"]],
           R1_C2_0 = ee1[["R0"]]*ee2[["C0"]],
           
           S1_R2_0 = ee1[["S0"]]*ee2[["R0"]],
           I1_R2_0 = ee1[["I0"]]*ee2[["R0"]],
           C1_R2_0 = ee1[["C0"]]*ee2[["R0"]],
           R1_R2_0 = ee1[["R0"]]*ee2[["R0"]]))
  
}

ic_vals <- SIRS2_independent_endemic_equilibrium(rp_vals)
params_all <- c(rp_vals,ic_vals)
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
    mutate(K1 = ifelse(time < 2.5, NA, K1), 
           total1 = params_all_int["rho1"]*K1,
           K2 = ifelse(time < 2.5, NA, K2), 
           total2 = params_all_int["rho2"]*K2) %>% 
    # select relevant state veariables: type specific new scaled cases of flu
    select(time, total1, total2) %>% 
    # store this in a long-format tibble 
    gather(key = pathogen, value = st_cases, -time) %>% 
    # join this with the target data: to calculate poisson log-likelihood
    right_join(target_data, by = c("time", "pathogen")) %>% 
    # calculate poisson log-density: P(data_t|model_t)
    mutate(poisson_log_density = dpois(x = obs_cases, lambda = st_cases, log = TRUE)) -> log_density_data 
  
  # This feature is for potential post-hoc analysis: if conditional log-liklihoods are called for
  if(give_conditional_log_lik == TRUE) {
    
    log_density_data %>% 
      # NA's are given a conditional log-density of 0
      replace_na(list(poisson_log_density = 0)) %>% 
      # Following two steps are used to generate type-specific conditional log-liklihood
      group_by(pathogen) %>% 
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





# create a sampled design of initial guesses: Sobol sampling
## this is for full parameter
 pd <- sobol_design(lower = c(R0_B = 1, R0_A = 1, amplitude_A = 1e-10, 
                              amplitude_B = 1e-10, tpeak_A = 1/365.25, tpeak_B = 1/365.25,  
                              rho_A = 1e-10, rho_B = 1e-10, chi =0, psi=0, mu = 1e-5), 
                    upper = c(R0_B = 5, R0_A = 5, amplitude_A = 1, amplitude_B = 1, 
                              tpeak_A = 1, tpeak_B = 1, rho_A = 0.01, rho_B = 0.01, chi =1, psi=1, mu =1e-2), 
                    nseq = 499)
 
write.csv(pd,"pd.csv",row.names = FALSE )

#####################################################################################
## parameter setting sepecific
#####################################################################################
##


# This function is defined to provide pomp-GA integration
## null model
# parameters of interest
est_these_null <- c("R01", "R02", "amplitude1", "amplitude2", "tpeak1", "tpeak2", 
                    "rho1", "rho2","mu")

## from loglikelihood function
of_ga_null<- function(x, est = est_these_null) {
  
  x<- unname(unlist(x))
  
  of(par = c(R01 = x[1], R02 = x[2], 
             amplitude1 = x[3], amplitude2 = x[4], 
             tpeak1 = x[5], tpeak2 = x[6],  
             rho1 = x[7], rho2 = x[8], mu=x[9]), est = est)
}

pd %>% 
  select(-psi,-chi)-> pd_withfix_null

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
             lower = c(R01 = 1, R02 = 1, amplitude1 = 1e-10, 
                       amplitude2 = 1e-10, tpeak1 = 1/365.25, tpeak2 = 1/365.25,  
                       rho1 = 1e-10, rho2 = 1e-10, mu = 1e-5), 
             upper = c(R01 = 5, R02 = 5, 
                       amplitude1 = 1, amplitude2 = 1, 
                       tpeak2 = 1, tpeak2 = 1, 
                       rho1 = 0.01, rho2 = 0.01, mu= 1),
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


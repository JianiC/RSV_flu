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

  
  //population size
  N = S1_S2 + I1_S2 + C1_S2 + R1_S2 + S1_I2 + I1_I2 + C1_I2 + R1_I2 + S1_C2 + 
      I1_C2 + C1_C2 + R1_C2 + S1_R2 + I1_R2 + C1_R2 + R1_R2;

  //I1_sum = I1_S2 + I1_I2 + I1_C2 + I1_R2;
  //I2_sum = S1_I2 + I1_I2 + C1_I2 + R1_I2;

  foi1 = beta1 * (I1_S2 + I1_I2 + I1_C2 + I1_R2 + eta1)/N; 
  foi2 = beta2 * (S1_I2 + I1_I2 + C1_I2 + R1_I2 + eta2)/N;
  
  
  //birth
  rate[0] = N*mu;
  
  // Out of S1_S2
  rate[1] = S1_S2 * foi1;         // infection of S1_S1 -> I1_S2
  rate[2] = S1_S2 * foi2;         // infection of S1_S2 -> S1_I2
  rate[3] = S1_S2 * mu;           // death from S1_S2  

  // Out of I1_S2
  rate[4] = I1_S2 * gamma1;             // recovery of I1_S2 -> C1_S2
  rate[5] = I1_S2 * psi * foi2;         // infection of IA_SB -> IA_IB 
  rate[6] = I1_S2 * mu;                // death from I1_S2  
  
  // Out of C1_S2
  rate[7] = C1_S2 * phi1;              // recovery of C1_S2 -> R1_S2
  rate[8] = C1_S2 * chi * foi2;        // infection of C1_S2 -> C1_I2
  rate[9] = C1_S2 * mu;               // death from C1_S2 
  
  // Out of R1_S2
  rate[10] = R1_S2 * w1;                // recovery of R1_S2->S1_S2
  rate[11] = R1_S2 * foi2;              // infection of R1_S2->R1_I2
  rate[12] = R1_S2 * mu;               // death from R1_S2 
  
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

make_pomp <- function(df, time_start_sim = -100) {
  
  # browser()
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
  
  
  dmeas <- dmeas_poisson
  rmeas <- rmeas_poisson
  
  # use populations sizes from the data
  Npop <- (
    df %.>% 
      select(., N) %.>% 
      slice(., n()) %.>% 
      un_list_name(.)
  )
  
  # regular parameters for the full model
  # Npop variable will be updated based on the data  
  rp_vals <- c(R01 = 1, gamma1=365./9, w1=1,
               R02 = 1, gamma2=365./3, w2=1,
               phi1=365/90, phi2=365/90, psi =1, chi=1, 
               eta1=365., eta2=365.,rho1 = 0, rho2 = 0, 
               amplitude1=0.0, tpeak1=0.0, amplitude2=0.0, tpeak2=0.0, 
               pop=Npop, 
               mu=1/80)
  
  
  
  pomp(
    data = df,
    times = "time",
    t0 = time_start_sim,
    obsnames = c("total1","total2"),
    skeleton = vectorfield(det_skel),
    dmeasure = dmeas,
    rmeasure = rmeas, 
    rinit = rinit_ee,
    accumvars =c("Kprim1","Ksec1","Kprim2","Ksec2","K1", "K2" ),
    statenames = c(model_variables, c("Kprim1","Ksec1","Kprim2","Ksec2","K1", "K2")),
    paramnames = c(model_params),
    params = rp_vals
  )
  
}


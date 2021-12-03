# dynamic simulation
## code to creat pomp model

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


# specify length of burn in for simulations (in years)
t_start <- 0 
# produce the simulation of compartment for three seasons

# add fake data to simulate three seasons 
data.frame(time = seq(0.5, 4.5, by = 1/52), 
           total_a = NA, 
           total_b = NA) %>% 
  make_pomp_model(time_start_sim=t_start) -> pomp_sirs





simulate_tss <- function(params, give_everything = FALSE, show_progress = TRUE,...) {
  # browser()
  if(show_progress == TRUE) {
    pb$tick()$print()  
  } else {
    print("Progress of the job will not be displayed!")
  }
  
  guess1_params <- c(R0_A=unname(params[,"R0"]), gamma_A=365./9,
                     R0_B=unname(params[,"rel"])*unname(params[,"R0"]), gamma_B=365./3.6,
                     w=1,
                     phi=unname(params[,"phi"]), 
                     chi=unname(params[,"chi"]), 
                     lamda=unname(params[,"lamda"]),
                     eta_A=365., eta_B=365., rho_A =0.0014, rho_B =  0.0004, 
                     sigmaSE=0.0001, psi=0.00001, 
                     amplitude_A=0.2238291, amplitude_B=0.183653, 
                     tpeak_A=0.8571905, tpeak_B = 0.1401407,
                     pop=2.6e7)
  
  guess1_params <- unlist(guess1_params) 
  
  guess1_ic <- SIRS2_independent_endemic_equilibrium(guess1_params)
  guess1_all <- c(guess1_params,guess1_ic)
  
  # browser()
  
  pomp_sirs %>%
    trajectory(params=guess1_all, t0=-100, format="d", method = "ode45") %>% 
    slice(2:n()) %>% 
    mutate(mp_A = max(K_A), 
           mp_B = max(K_B), 
           mp_AB_ratio = mp_A/mp_B, 
           pw_A = time[which(K_A == max(K_A))], 
           pw_B = time[which(K_B == max(K_B))], 
           p_AB_diff = (pw_A - pw_B)*52, 
           pop = guess1_all["pop"],
           S_A_tot = SA_SB + SA_IB + SA_CB + SA_RB,
           S_B_tot = SA_SB + IA_SB + CA_SB + RA_SB,
           I = IA_SB + SA_IB + IA_IB + IA_RB + RA_IB, 
           I_A_tot = IA_SB + IA_IB + IA_RB, 
           I_B_tot = SA_IB + IA_IB + RA_IB,   
           C = CA_SB + SA_CB, 
           R = SA_RB + RA_SB + RA_RB, 
           I_A_prop = (IA_SB + IA_IB + IA_RB)/I, 
           I_B_prop = (SA_IB + IA_IB + RA_IB)/I, 
           C_A_prop = CA_SB/C, 
           C_B_prop = SA_CB/C) -> everything 
  
  everything %>% 
    slice(n()) %>% 
    select(mp_AB_ratio, p_AB_diff) -> test_sim
  
  if(give_everything == TRUE) {
    print("All the Compartments are produced")
    return(everything)
  } else {
    return(test_sim)  
  }
  
} 


#function to loop over values 
multi_simulate_tss <- function(counter, params_mat, ...) {
  simulate_tss(params = params_mat[counter,], ...)
}



trajectory_plot<-function(data,test.labs){
  data %>% 
    mutate(
      facet_label = rep(c("a", "b", "c"), each = length(seq(0.5, 4.5-1/52, by = 1/52))),
      Inc_A_tot = K_A, 
      Inc_B_tot = K_B
    ) %>% 
    select(Inc_A_tot, Inc_B_tot, time, facet_label) %>% 
    gather(key = "Compartment", value = "Count", -c(time, facet_label)) %>% 
    ggplot(aes(x = time+2017, y = Count, colour = Compartment, fill = Compartment)) +
    geom_area(position = position_dodge(width = 0), alpha = 0.5) +
    labs(x = "Time ", 
         y = "Cases ") +
    scale_colour_manual(name = "", 
                        values = c("red", "blue")) +
    facet_wrap(.~ facet_label, scales = "fixed", ncol = 1,labeller = labeller(facet_label=test.labs))+
    theme_bw()+
    #theme(legend.position="none")+
    scale_fill_manual(name = "", 
                      values = c("red", "blue"),
                      labels=c("RSV","IAV"))-> c_grid_plot
  
  return(c_grid_plot)
  
  
}


## plot different compartment

S_compart_plot<-function(data){
  c_facet_data %>% 
    mutate(
      facet_label = rep(c("a", "b", "c"), each = length(seq(0.5, 4.5-1/52, by = 1/52))),
      Inc_A_tot = K_A, 
      Inc_B_tot = K_B
    ) %>% 
    select(S_A_tot, S_B_tot, time, facet_label) %>% 
    gather(key = "Compartment", value = "Count", -c(time, facet_label)) %>% 
    ggplot(aes(x = time+2017, y = Count, colour = Compartment, fill = Compartment)) +
    geom_line()+
    #geom_area(position = position_dodge(width = 0), alpha = 0.5) +
    labs(x = "Time ", 
         y = "Cases ") +
    scale_colour_manual(name = "", 
                        values = c("red", "blue")) +
    facet_wrap(.~ facet_label, scales = "fixed", ncol = 1,labeller = labeller(facet_label=test.labs))+
    theme_bw()+
    #theme(legend.position="none")+
    scale_fill_manual(name = "", 
                      values = c("red", "blue"),
                      labels=c("RSV","IAV"))-> c_grid_plot_S
  
  return(c_grid_plot_S)
  
  
}




library(ggplot2)



## 1. impact of R0 of IAV

param_values <- data.frame(R0 = 2., chi =1, lamda =1 ,phi = 365/30,  
                           rel = c(1,1.2,0.8))


c_facet_data <- map_df(1:3, multi_simulate_tss, params = param_values, 
                       give_everything = TRUE, show_progress = FALSE)

test.labs <-c("a"="relative R0 = 1","b"="relative R0 = 1.2","c"="relative R0 = 0.8")
plot1<-trajectory_plot(c_facet_data,test.labs)
plot1


## 3. impact of lamda

param_values <- data.frame(R0 = 2., chi =1, lamda = c(1, 0.5, 0.2), phi = 365/30,  
                           rel = 0.9)

test.labs <-c("a"="lamda = 1","b"="lamda=0.5","c"="lamda = 0.2")

c_facet_data <- map_df(1:3, multi_simulate_tss, params = param_values, 
                       give_everything = TRUE, show_progress = FALSE)


plot1<-trajectory_plot(c_facet_data,test.labs)
plot1


## 3. impact of chi

param_values <- data.frame(R0 = 2., chi = c(1,0.5,0.2), lamda=1, phi =365/30, rel = 0.9)

test.labs <-c("a"="chi = 1","b"="chi=0.5","c"="chi = 0.2")
c_facet_data <- map_df(1:3, multi_simulate_tss, params = param_values, 
                       give_everything = TRUE, show_progress = FALSE)

test.labs <-c("a"="chi = 1","b"="chi =0.9","c"="chi = 0.8")
plot1<-trajectory_plot(c_facet_data,test.labs)
plot1



## 3. impact of chi

param_values <- data.frame(R0 = 2., chi = 0.8, lamda=1, phi =c(365/30,365/90,365/300), rel = 0.9)

test.labs <-c("a"="chi = 1","b"="chi=0.5","c"="chi = 0.2")
c_facet_data <- map_df(1:3, multi_simulate_tss, params = param_values, 
                       give_everything = TRUE, show_progress = FALSE)

test.labs <-c("a"="phi = 365/30","b"="phi =365/90","c"="phi = 365/300")
plot1<-trajectory_plot(c_facet_data,test.labs)
plot1
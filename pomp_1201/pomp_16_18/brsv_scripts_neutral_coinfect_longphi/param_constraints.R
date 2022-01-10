# this script contains param constraints of all the hypothesis 

# hypothesis neutral
neutral_param_constraints <- list(
  lower = c( R01 = 1, R02 = 1, 
            amplitude1 = 1e-10, amplitude2 = 1e-10, 
            tpeak1 = 1/365.25, tpeak2 = 1/365.25,
            rho1 = 1e-10, rho2 = 1e-10),
  upper = c(R01 = 5, R02 =5, 
            amplitude1 = 1, amplitude2 = 1, 
            tpeak1 = 1, tpeak2 = 1,
            rho1 = 0.01, rho2 = 0.01)
)


co_infect_param_constraints <- list(
  lower = c( R01 = 1, R02 = 1, 
             amplitude1 = 1e-10, amplitude2 = 1e-10, 
             tpeak1 = 1/365.25, tpeak2 = 1/365.25,
             rho1 = 1e-10, rho2 = 1e-10,
             psi = 1e-5, chi =1e-5 ),
  upper = c(R01 = 5, R02 = 5, 
            amplitude1 = 1, amplitude2 = 1, 
            tpeak1 = 1, tpeak2 = 1,
            rho1 = 0.01, rho2 = 0.01,
            psi = 1, chi = 1)
)

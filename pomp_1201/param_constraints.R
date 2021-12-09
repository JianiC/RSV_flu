# this script contains param constraints of all the hypothesis 

# hypothesis co-infection
coinf_param_constraints <- list(
  lower = c(psi = 0, R01 = 1, R02 = 1, rho1 = 0, rho2 = 0, 
            amplitude1 = 0, amplitude2 = 0, 
            tpeak1 = 1, tpeak2 = 1),
  upper = c(psi = 10, R01 = 10, R02 = 10, rho1 = 1, rho2 = 1, 
            amplitude1 = 1, amplitude2 = 1, 
            tpeak1 = 1, tpeak2 = 1)
)



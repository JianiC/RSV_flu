# load packages 
packages <- c("tidyverse", "stringr",
              "pomp", "magrittr", "wrapr",
              "tictoc", "parallel", 
              "LaplacesDemon", 
              "rootSolve",
              "DEoptim", "doParallel", "doRNG", 
              "parallel")


lapply(packages, library, character.only = TRUE)

# load the requisite data 
load("./inc_data_add.rds")

# set up controls for deoptim
### controls list for DEoptim - for MLE estimation
np_val = 30
my_controls <- list(itermax = 10,
                    F = 0.6, CR = 0.9, 
                    strategy = 1,
                    steptol = 750, reltol = 1e-8)

est_past_est<-data.frame(R01=2.012,
                         R02=1.243,
                         amplitude1=0.116,
                         amplitude2=0.126,
                         tpeak1=0.938,
                         tpeak2=0.022,
                         rho1=0.000253,
                         rho2=0.000585)
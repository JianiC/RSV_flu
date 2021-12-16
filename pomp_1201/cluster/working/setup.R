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
load("./inc_data.rds")

# set up controls for deoptim
### controls list for DEoptim - for MLE estimation
np_val = 30
my_controls <- list(itermax = 10,
                    F = 0.6, CR = 0.9, 
                    strategy = 1,
                    steptol = 750, reltol = 1e-8)


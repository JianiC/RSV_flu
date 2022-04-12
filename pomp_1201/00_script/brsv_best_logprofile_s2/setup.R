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
# for the first try - range 0, 1 on chi
# steptol = 200
# np_val = 150

# if the parameter is poorly identified i.e. the likelihood surface is flat
# then {
# for the second try - range 0.75*mle, 1.25*mle on chi
# steptol = 200
# np_val = 150
# } else {

#for the second try - range 0.75*mle, 1.25*mle on chi
# increase the strictness of the optimization problem
# steptol = 500
# np_val = 500
# }

# finally check if the increasing further values makes any difference 
# basically 
# steptol = 750
# np_val = 750


np_val = 1700
my_controls <- list(itermax = 1e6,
                    F = 0.6, CR = 0.9, 
                    strategy = 1,
                    steptol = 1500, reltol = 1e-8)

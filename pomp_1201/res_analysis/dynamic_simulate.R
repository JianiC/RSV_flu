## pomp_model simulation plot
## plot the effect of different psi and chi value

## load pomp_module
setwd("pomp_1201/res_analysis")
source("./fit_functions.R", chdir = TRUE) 
source("./util_funs.R", chdir = TRUE) 
## load pomp module

simulate_tss <- function(params, res_hhs=res_hhs3_arsv_coinfect, give_everything = FALSE, show_progress = TRUE,...) {
  # browser()
  if(show_progress == TRUE) {
    pb$tick()$print()  
  } else {
    print("Progress of the job will not be displayed!")
  }
  
  guess_params <- c(R01=unname(params[,"R01"]), gamma1=365./9,  w1=2,
                     R02=unname(params[,"R02"]), gamma2=365./3, w2=2,
                     phi1=365/180, phi2=365/180,
                     psi=unname(params[,"psi"]), chi=unname(params[,"chi"]), 
                     eta1=365., eta2=365., rho1 =unname(params[,"rho1"]), rho2 =unname(params[,"rho2"]), 
                     amplitude1=unname(params[,"amplitude1"]), amplitude2=unname(params[,"amplitude2"]),
                     tpeak1=unname(params[,"tpeak1"]), tpeak2=unname(params[,"tpeak2"]),
                     pop=inc_data_add$N[3],
                     mu=1/80)
  
  #guess_params <- unlist(guess_params) 

  ## make pomp object
  tibble(time = seq(0.5, 1.5, by = 1/52), 
         total1 = NA, 
         total2 = NA, N = inc_data_add$N[3]) %>% 
    make_pomp(time_start_sim = -100) -> pomp_sirs

  
  
  # browser()
  
  
  pomp_sirs %>%
    trajectory(params=guess_params, t0=-100,format = "d", method = "ode23") %>% 
    slice(2:n()) %>% 
    mutate(mp1 = max(K1), 
           mp2 = max(K2), 
           mp12_ratio = mp1/mp2, 
           pw1 = which(K1 == max(K1)), 
           pw2 = which(K2 == max(K2)), 
           pw12_diff = pw1 - pw2) -> everything 
  
  everything %>% 
    slice(n()) %>% 
    select(mp12_ratio, pw12_diff) -> test_sim
  
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

##################################################
# test the effect of psi and chi
sim_params <- expand.grid(chi = seq(0,1,by =0.02), 
                             psi = seq(0,1,by =0.02))


## load example result
load("pomp_longphi_result/best_model_rsva/res_hhs3_arsv_coinfect.rds")


ref_df<-past_est(res_hhs=res_hhs3_arsv_coinfect)%>%
  select(-psi,-chi)
   
params_design<-sim_params%>% bind_cols(ref_df)

pb <- progress_estimated(nrow(params_design))
tic()
library(future)
library(furrr)
plan(multisession)
res_design<- future_map_dfr(.x = 1:nrow(params_design ), .f = multi_simulate_tss,
                         params_mat = params_design, .progress = TRUE)
toc()



sim_res<-params_design %>% bind_cols(res_design)
ggplot(sim_res, aes(x=chi,y=psi,fill=mp12_ratio))+
  geom_tile()+
  gg.theme+
  labs(y=expression(paste("Propotion of inhibition on co-infection(",psi,")")),
       x=expression(paste("Strength of cross-immunity((",chi,")")),
       fill="peak case ration\n(pathogen1/pathogen2)")+
  guides(colorbar = guide_legend(title.position = "top"))+
  theme(legend.position = "top",
        text = element_text(size = 8),
        legend.text = element_text(size=6))+
  scale_fill_viridis("Peak case ration\n(pathogen1/pathogen2)",discrete=FALSE)->mp_sim


ggplot(sim_res, aes(x=chi,y=psi,fill=pw12_diff))+
  geom_tile()+
  gg.theme+
  labs(y=expression(paste("Propotion of inhibition on co-infection (",psi,")")),
       x=expression(paste("Strength of cross-immunity(",chi,")")),
       fill="Peak week difference\n(pathogen1-pathogen2)")+
  theme(legend.position = "bottom",
        text = element_text(size = 8),
        legend.text = element_text(size=6))+
  guides(colorbar = guide_legend(title.position = "top"))+
  scale_fill_viridis("Peak week difference\n(pathogen1-pathogen2)",discrete=FALSE)->pw_sim


ggarrange(mp_sim +
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  plot.margin = margin(t=1,b = 1) ),
            
          pw_sim + 
            theme(axis.ticks.x = element_blank(),
                  plot.margin = margin(b=1,t = 1) ), 
          nrow = 2)->dynamic_simu








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
  tibble(time = seq(0.5, 3.5, by = 1/52), 
         total1 = NA, 
         total2 = NA, N = inc_data_add$N[3]) %>% 
    make_pomp(time_start_sim = -100) -> pomp_sirs
  
  
  
  # browser()
  
  
  pomp_sirs %>%
    trajectory(params=guess_params, t0=-100,format = "d", method = "ode23") %>% 
    slice(2:n()) %>% 
    mutate(case1=K1, 
           case2=K2) -> everything 
  
  everything %>% 
    slice(n()) %>% 
    select(case1, case2) -> test_sim
  
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

##############################################################
## plot the curve with different psi and chi value


point_params<-data.frame(psi=c(1,1,0.3,0.3),chi=c(1,0.5,1,0.5))
pointparams_design<-point_params%>% bind_cols(ref_df)


pointsim_result<- map_df(1:4, multi_simulate_tss, params = pointparams_design, 
                       give_everything = TRUE, show_progress = FALSE)


p_label <- tibble(
  facet_label = c("a", "b", "c","d"),
  label = c(expression(paste(psi==0,",", chi==0)),expression(paste(psi==0,",", chi==0.5)),
            expression(paste(psi==0.7,",", chi==0)),expression(paste(psi==0.7,",", chi==0.5)))
  
)

pointsim_result%>%
  mutate(case1=case1*pointparams_design$rho1[1],
         case2=case2*pointparams_design$rho2[2])%>%
  mutate(scase1=(case1/inc_data_add$N[3])*1e7,
         scase2=(case2/inc_data_add$N[3])*1e7)%>%
  mutate(facet_label = rep(c("a", "b", "c","d"), each = length(seq(0.5, 3.5-1/52, by = 1/52))))%>%
  select(scase1,scase2,time,facet_label)%>%
  gather(key = "Compartment", value = "Count", -c(time, facet_label))->facet_result

facet_result%>%
  ggplot(aes(x = time+2013, y = Count, colour = Compartment, fill = Compartment)) +
  geom_area(position = position_dodge(width = 0), alpha = 0.4) +
  geom_text(data= p_label,aes(x = 2014.2, y =300 ,label = label),inherit.aes = FALSE,parse = T,size=2.5)+
  labs(x = "Time (weeks) ", 
       y = "Cases / million") +
  scale_colour_manual(name = "", 
                      values = c("#66A61E", "#E6AB02"),
                      labels=c("pathogen1","pathogen2")) +
  scale_fill_manual(name = "", 
                      values = c("#66A61E", "#E6AB02"),
                    labels=c("pathogen1","pathogen2"))+
  facet_wrap(.~ facet_label, scales = "fixed", ncol = 1)+
  gg.theme+
  theme(  strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position = "top",
          text = element_text(size = 8))->point_simu


plot_grid(dynamic_simu,point_simu,
          nrow = 1,
          labels = "AUTO",label_size = 12,
          rel_widths = c(1,1.2))->simu_plot  
          

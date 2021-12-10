# this script runs different functions to make sure that they are running alright

# first load prerequisites
source("./src.R", chdir = TRUE)

## specify the virus and HHS region

# make data ready for pomp
pomp_data_hhs1_arsv <- (
  inc_data %.>% 
  make_data_pomp_ready(., virus_combo = c("RSV", "fluA"), HHS_region = 1)
  )



if(FALSE) {
  
  
  
  pseudo_data <- tibble(time = seq(0, 10, by = 1/52), 
                        total1 = NA, 
                        total2 = NA, N = pomp_data_hhs1_arsv$N[1])
  
  # make a pomp object 
  hhs1_a_rsv_po <- make_pomp(df = pseudo_data, time_start_sim = -100)
  # inspect the compiled pomp object
  # spy(hhs1_a_rsv_po)
  
# test if the integratro works as default
test_traj <- trajectory(object = hhs1_a_rsv_po, format = "d", method = "ode23")

plot_comp <- (
  test_traj %.>% 
    slice(., 2:n()) %.>% 
    select(., -`.id`) %.>% 
    gather(., "comp", "count", -time) %.>% 
    ggplot(., aes(x = time, y = count)) +
    geom_line()+
    facet_wrap(.~comp, scales = "free")
)
}

# loading parameter constraints 


res <- (
  DE_traj_match(df = pomp_data_hhs1_arsv, 
                param_constraints = neutral_param_constraints, 
                params = rp_vals_def,
                ode_control = list(method = "ode23"), 
                hypo_name = "neutral", 
                hhs_reg = 1, 
                tot1_name = "RSV", 
                tot2_name = "fluA")
)

if(res$total2 == "RSV") message("Code itegration complete!")

save(res, file = "res_hhs1_arsv.Rdata")



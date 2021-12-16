# first load prerequisites
source("./src.R", chdir = TRUE)

# make data ready for pomp
pomp_data_hhs1_arsv <- (
  inc_data %.>% 
    make_data_pomp_ready(., virus_combo = c("RSV", "fluA"), HHS_region = 1)
)

hhs1_a_rsv_po <- make_pomp(df = pomp_data_hhs1_arsv, time_start_sim = -100)


as.data.frame(res_hhs1_arsv_neutral$DEobj$optim$bestmem)%>%
  cbind(row.names(.))%>%
  pivot_wider(names_from = `row.names(.)`, values_from = `res_hhs1_arsv_neutral$DEobj$optim$bestmem`)%>%
  mutate(phi1=365/30, phi2=365/30,psi=1, chi =1)->param_vector


as.data.frame(res_hhs1_arsv_coinfect$DEobj$optim$bestmem)%>%
  cbind(row.names(.))%>%
  pivot_wider(names_from = `row.names(.)`, values_from = `res_hhs1_arsv_coinfect$DEobj$optim$bestmem`)%>%
  mutate(phi1=365/30, phi2=365/30)->param_vector

rp_vals <- c(R01=param_vector$R01, gamma1=365./9, w1=1/1,
             R02=param_vector$R02, gamma2=365./3, w2=1/1,
             amplitude1 = param_vector$amplitude1,
             amplitude2 = param_vector$amplitude2,
             tpeak1 = param_vector$tpeak1,
             tpeak2 = param_vector$tpeak2,
             phi1=param_vector$phi1, phi2=param_vector$phi2,
             psi=param_vector$psi, chi=param_vector$chi,
             eta1=365., eta2=365.,
             rho1=param_vector$rho1, rho2=param_vector$rho2,
             sigmaSE=0.0000,
             pop= pomp_data_hhs1_arsv$N[1],
             mu=1/80)




test_traj <- trajectory(object = hhs1_a_rsv_po, params = rp_vals, format = "d", method = "ode23")

test_traj %>% 
  slice(., 2:n()) %>% 
  #select(., -`.id`) %.>% 
  select(time, K1, K2) %>%
  mutate(scases1 = param_vector$rho1*K1, 
         scases2 = param_vector$rho2*K2)%>%
  select(time, scases1, scases2) %>%
  gather(key = "type", value = "cases", -time)->test_traj_data 

test_traj_data %>%
  ggplot(aes(x = time, y = cases)) +
  geom_line()+
  facet_wrap(.~type, scales = "free")

pomp_data_hhs1_arsv%>%
  mutate(scases1=total1,
         scases2=total2)%>%
  select(time, scases1, scases2) %>%
  gather(key = "type", value = "obscases", -time) %>%
  full_join(test_traj_data, by =c ("time" = "time","type"="type"))->fit_withdata

fit_withdata%>%
  ggplot(aes(x = time, y = cases))+
  geom_line()+
  facet_wrap(.~type, scales = "free")+
  geom_line(aes(y=obscases), color = "red")



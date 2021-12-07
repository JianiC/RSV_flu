source("./setup.R")
source("./util_funs.R")

# load the data csv

data <- read.csv("./pomp_HHSregion_case.csv")

demog <- read.csv("./pop_HHSRegion.csv") %.>% 
  mutate(., HHS_REGION = HHS_region) %.>%
  select(., 4,3) 


# tidy inc data
inc_data <- (
  data %.>% 
    rbind(., 
          tibble(date = data$date %.>% min(.) - 1/52,
                 time = data$time %.>% min(.) - 1/52, 
                 RSVpos =  NA, 
                 fluApos = NA, 
                 fluBpos = NA, 
                 HHS_REGION = 1:10)) %.>% 
    right_join(., demog, by = "HHS_REGION") %.>% 
    arrange(., time) %.>% 
    gather(., key = "virus", value = "cases", -c(time, date, HHS_REGION, pop_ave_sum)) %.>% 
    mutate(.,
           N = pop_ave_sum,
           virus = str_remove(virus, "pos")) %.>% 
    select(., -pop_ave_sum)
  )

# save(inc_data, file = "inc_data.rds")





  
  





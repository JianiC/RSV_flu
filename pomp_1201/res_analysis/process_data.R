source("./setup.R")
source("./util_funs.R")
library(lubridate)
# load the data csv

data <- read.csv("../original_data/pomp_HHSregion_case.csv")
data <- read.csv("../original_data/pomp_HHSregion_case_new.csv")
demog <- read.csv("../original_data/pop_HHSRegion.csv") %.>% 
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
time_add<-c(decimal_date(ymd("2014-01-04")),decimal_date(ymd("2015-01-03")),decimal_date(ymd("2016-01-02")))

## add missing date point to inc data



# save(inc_data, file = "inc_data.rds")

time_add<-c(decimal_date(ymd("2014-01-04")),decimal_date(ymd("2015-01-03")),
            decimal_date(ymd("2016-01-02")),decimal_date(ymd("2018-01-06"))
                         )

inc_data_add <- (
  data %.>% 
    rbind(., 
          tibble(date = data$date %.>% min(.) - 1/52,
                 time = data$time %.>% min(.) - 1/52, 
                 RSVpos =  NA, 
                 fluApos = NA, 
                 fluBpos = NA, 
                 HHS_REGION = 1:10)) %.>%
    rbind(.,
          tibble(date=time_add[1],
                 time=time_add[1]-2011,
                 RSVpos =  NA, 
                 fluApos = NA, 
                 fluBpos = NA, 
                 HHS_REGION = 1:10))%.>%
    rbind(.,
          tibble(date=time_add[2],
                 time=time_add[2]-2011,
                 RSVpos =  NA, 
                 fluApos = NA, 
                 fluBpos = NA, 
                 HHS_REGION = 1:10))%.>%
    rbind(.,
          tibble(date=time_add[3],
                 time=time_add[3]-2011,
                 RSVpos =  NA, 
                 fluApos = NA, 
                 fluBpos = NA, 
                 HHS_REGION = 1:10))%.>%
    right_join(., demog, by = "HHS_REGION") %.>% 
    arrange(., time) %.>% 
    gather(., key = "virus", value = "cases", -c(time, date, HHS_REGION, pop_ave_sum)) %.>% 
    mutate(.,
           N = pop_ave_sum,
           virus = str_remove(virus, "pos")) %.>% 
    select(., -pop_ave_sum)
)


save(inc_data_add, file = "inc_data_add2.rds")

data_perdic <- read.csv("../original_data/pomp_HHSregion_case_2018.csv")

inc_data_perdic <- (
  data_perdic %.>% 
    rbind(., 
          tibble(date = data$date %.>% min(.) - 1/52,
                 time = data$time %.>% min(.) - 1/52, 
                 RSVpos =  NA, 
                 fluApos = NA, 
                 fluBpos = NA, 
                 HHS_REGION = 1:10)) %.>%
    rbind(.,
          tibble(date=time_add[1],
                 time=time_add[1]-2011,
                 RSVpos =  NA, 
                 fluApos = NA, 
                 fluBpos = NA, 
                 HHS_REGION = 1:10))%.>%
    rbind(.,
          tibble(date=time_add[2],
                 time=time_add[2]-2011,
                 RSVpos =  NA, 
                 fluApos = NA, 
                 fluBpos = NA, 
                 HHS_REGION = 1:10))%.>%
    rbind(.,
          tibble(date=time_add[3],
                 time=time_add[3]-2011,
                 RSVpos =  NA, 
                 fluApos = NA, 
                 fluBpos = NA, 
                 HHS_REGION = 1:10))%.>%
    rbind(.,
          tibble(date=time_add[4],
                 time=time_add[4]-2011,
                 RSVpos =  NA, 
                 fluApos = NA, 
                 fluBpos = NA, 
                 HHS_REGION = 1:10))%.>%
    right_join(., demog, by = "HHS_REGION") %.>% 
    arrange(., time) %.>% 
    gather(., key = "virus", value = "cases", -c(time, date, HHS_REGION, pop_ave_sum)) %.>% 
    mutate(.,
           N = pop_ave_sum,
           virus = str_remove(virus, "pos")) %.>% 
    select(., -pop_ave_sum)
)


save(inc_data_perdic, file = "inc_data_perdic.rds")


  
inc_data_1618 <- (
  data_perdic %.>%
    filter(.,data_perdic>=2015.5) %.>%
    rbind(., 
          tibble(date = data_perdic$date %.>% min(.) - 1/52,
                 time = data_perdic$time %.>% min(.) - 1/52, 
                 RSVpos =  NA, 
                 fluApos = NA, 
                 fluBpos = NA, 
                 HHS_REGION = 1:10)) %.>%
    rbind(.,
          tibble(date=time_add[3],
                 time=time_add[3]-2011,
                 RSVpos =  NA, 
                 fluApos = NA, 
                 fluBpos = NA, 
                 HHS_REGION = 1:10))%.>%
    rbind(.,
          tibble(date=time_add[4],
                 time=time_add[4]-2011,
                 RSVpos =  NA, 
                 fluApos = NA, 
                 fluBpos = NA, 
                 HHS_REGION = 1:10))%.>%
    right_join(., demog, by = "HHS_REGION") %.>% 
    arrange(., time) %.>% 
    gather(., key = "virus", value = "cases", -c(time, date, HHS_REGION, pop_ave_sum)) %.>% 
    mutate(.,
           N = pop_ave_sum,
           virus = str_remove(virus, "pos")) %.>% 
    select(., -pop_ave_sum)
)


save(inc_data_1618 , file = "inc_data_1618.rds")



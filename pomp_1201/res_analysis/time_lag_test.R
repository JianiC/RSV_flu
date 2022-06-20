## time lagged cross-correlation test
Spend <- c(5, 3, 6, 5, 8, 9, 10, 17, 12, 11, 10, 9)
Income <- c(25, 29, 22, 34, 22, 28, 29, 31, 34, 45, 45, 40)
ccf(Spend, Income)
test<-ccf(Spend, Income)
test$lag
test$acf
test$n.used
lag2.plot (soi, rec, 10)

source("./fit_functions.R", chdir = TRUE) 
library(ggpubr)
library("lubridate")
RSV_flu_surveliance <-read.csv("RSV_flu_combine.csv")

pop<-read.csv("../../demographic-data/nst-est2019-alldata.csv")
region_state<-read.csv("../../demographic-data/state_region.csv")

RSV_flu_national<-read_csv("RSV_flu_combine_national.csv")

pop%>%
  left_join(region_state,by=c("US.STATE"="State_name"))%>%
  group_by(HHS_region)%>%
  summarise(across(6:27, sum))%>%
  filter(HHS_region !=0)%>%
  select(HHS_region,starts_with("POPEST"))%>%
  gather(year,pop,POPESTIMATE2010:POPESTIMATE2019,factor_key=TRUE)%>%
  mutate(year=as.numeric(gsub("^POPESTIMATE", "",year)))->pop_hhsregion

RSV_flu_surveliance%>% 
  mutate(year=year(RepWeekDate.y))%>%
  right_join(pop_hhsregion, by=c("HHS_REGION"="HHS_region","year"="year"))%>%
  mutate(RSVpos_norm=(RSVpos/pop)*100000,
         fluApos_norm=(fluApos/pop)*100000,
         fluBpos_norm=(fluBpos/pop)*100000)%>%
  select(RSVpos_norm,fluApos_norm,fluBpos_norm,HHS_REGION,date)->RSV_flu_surveliance_norm

## test with HHS region1 RSV_fluA
RSV_flu_surveliance_norm1<-RSV_flu_surveliance_norm%>%filter(as.numeric(HHS_REGION)==1)

test<-ccf(RSV_flu_surveliance_norm1['RSVpos_norm'], RSV_flu_surveliance_norm1['fluApos_norm'],30)
#remotes::install_github("nickpoison/astsa/astsa_build")
library(astsa)
lag2.plot(as.numeric(RSV_flu_surveliance_norm1['RSVpos_norm']), as.numeric(RSV_flu_surveliance_norm1['fluApos_norm']),5)

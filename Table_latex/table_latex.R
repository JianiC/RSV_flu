install.packages("kableExtra")
library(kableExtra)
library(dplyr)
## Table 1: Hypothesis table
hypo<-read.table("Table_latex/hypothesis.txt",sep = "\t",header = T)
knitr::kable(hypo,"latex",booktabs = TRUE)




#####################################################################################
## Table 2: parameter estimates in HHS region1

res_arsv<-read.csv("pomp_1201/res_analysis/results/result_res_arsv_longphhi_check.csv")
res_brsv<-read.csv("pomp_1201/res_analysis/results/result_res_brsv_longphhi_check.csv")
arsv_comp<-read.csv("pomp_1201/res_analysis/results/arsv_modelcompare.csv")
brsv_comp<-read.csv("pomp_1201/res_analysis/results/brsv_modelcompare.csv")

res<-rbind(res_arsv,res_brsv,arsv_comp,brsv_comp)

## rename column and parameter

res_all<-res%>%rename("Flu subtype"="Pathogen2",
             "HHS region"="HHSregion",
             "Parameter"="variable",
             "No-interaction"="neutral",
             "Inhibition on co-infection"="psi",
             "Cross protection"="chi",
             "Co-infection and cross-protection"="co_infect")%>%
  distinct()%>%
  mutate(`No-interaction`=ifelse(startsWith(Parameter,"rho"),formatC(`No-interaction`,format = "e",digits = 2),round(`No-interaction`,digits = 2)))%>%
  mutate(`Inhibition on co-infection`=ifelse(startsWith(Parameter,"rho"),formatC(`Inhibition on co-infection`,format = "e",digits = 2),round(`Inhibition on co-infection`,digits = 2)))%>%
  mutate(`Cross protection`=ifelse(startsWith(Parameter,"rho"),formatC(`Cross protection`,format = "e",digits = 2),round(`Cross protection`,digits = 2)))%>%
  mutate(`Co-infection and cross-protection`=ifelse(startsWith(Parameter,"rho"),formatC(`Co-infection and cross-protection`,format = "e",digits = 2),round(`Co-infection and cross-protection`,digits = 2)))%>%
  left_join(data.frame(Parameter = parameter_order),  
            .,
            by = "Parameter")%>%
  arrange(`Flu subtype`,`HHS region`)%>%
  mutate(Parameter = case_when(Parameter=="R01"~"R_0^{RSV}",
                              Parameter=="R02"~"R_0^{FLu}",
                              Parameter=="amplitude1"~"b_{RSV}",
                              Parameter=="amplitude2"~"b_{FLu}",
                              Parameter=="tpeak1"~"t_0^{RSV}",
                              Parameter=="tpeak2"~"t_0^{FLu}",
                              Parameter=="rho1"~"rho_{RSV}",
                              Parameter=="rho2"~"rho_{Flu}",
                              Parameter=="psi"~"psi",
                              Parameter=="chi"~"chi",
                              Parameter=="loglik"~"loglik",
                              Parameter=="AIC"~"AIC",
                              Parameter=="deltaAIC"~"deltaAIC",
                              Parameter=="R2"~"R^2",
                              Parameter=="RMSE"~"RMSE"))

parameter_order<-c("R01","R02","amplitude1","amplitude2","tpeak1","tpeak2","rho1","rho2","psi","chi","loglik","AIC","deltaAIC","R2","RMSE")  
res_HHS1<-res_all%>% filter(`HHS region`==1)%>%
  select(-`Flu subtype`,-`HHS region`)
  
knitr::kable(res_HHS1, "latex",booktabs=T) ->res_HHS1_latex 

res_HHS1_latex%>%
  pack_rows("RSV - Flu A",1,15)%>%
  pack_rows("RSV - Flu B",16,30)     ->t_res_HHS1
###########################################################################
## Table 3:AIC and goodness of fit among four epidemiological hypotheses among 10 HHS regions

model_compare_arsv<-read.csv("pomp_1201/res_analysis/results/arsv_modelcompare.csv")
model_compare_brsv<-read.csv("pomp_1201/res_analysis/results/brsv_modelcompare.csv")
parameter_order2<-c("deltaAIC","R2","RMSE")

rbind(model_compare_arsv,model_compare_brsv)%>%
  rename("Flu subtype"="Pathogen2",
            "HHS region"="HHSregion",
            "Parameter"="variable",
            "No-interaction"="neutral",
            "Inhibition on co-infection"="psi",
            "Cross protection"="chi",
            "Co-infection and cross-protection"="co_infect")%>%
  mutate(`No-interaction`=ifelse(startsWith(Parameter,"rho"),formatC(`No-interaction`,format = "e",digits = 2),round(`No-interaction`,digits = 2)))%>%
  mutate(`Inhibition on co-infection`=ifelse(startsWith(Parameter,"rho"),formatC(`Inhibition on co-infection`,format = "e",digits = 2),round(`Inhibition on co-infection`,digits = 2)))%>%
  mutate(`Cross protection`=ifelse(startsWith(Parameter,"rho"),formatC(`Cross protection`,format = "e",digits = 2),round(`Cross protection`,digits = 2)))%>%
  mutate(`Co-infection and cross-protection`=ifelse(startsWith(Parameter,"rho"),formatC(`Co-infection and cross-protection`,format = "e",digits = 2),round(`Co-infection and cross-protection`,digits = 2)))%>%
  left_join(data.frame(Parameter = parameter_order2),  
            .,
            by = "Parameter")%>%
  arrange(`Flu subtype`,`HHS region`)%>%
  select(`HHS region`,`Parameter`,`No-interaction`,`Inhibition on co-infection`,`Cross protection`,`Co-infection and cross-protection`)->model_compare



knitr::kable(model_compare,"latex",booktabs = TRUE)->model_compare_latex
model_compare_latex %>% 
  pack_rows("RSV - Flu A",1,15,latex_gap_space="2em")%>%
  pack_rows("RSV - Flu B",16,30,latex_gap_space="2em")%>%
  collapse_rows(columns = 1, latex_hline = "major", valign = "middle")  

#######################################################################
## Table 4: Parameters value used in POMP model


pomp_par<-read.csv("Table_latex/pomp_param.csv")
knitr::kable(pomp_par,"latex",booktabs = TRUE)


  
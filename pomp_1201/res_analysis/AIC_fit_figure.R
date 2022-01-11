result_res_arsv<-read.csv("./figures/result_res_arsv_longphhi.csv")
result_res_brsv<-read.csv("./figures/result_res_brsv_longphhi.csv")

hypo_levels=c("neutral","psi","chi","co_infect")

rbind(result_res_arsv,result_res_brsv) %>%
  filter(variable=="AIC")%>%
  mutate(
    psi=psi-neutral,
    chi=chi-neutral,
    co_infect=co_infect-neutral,
    neutral=neutral-neutral)%>%
  gather(hypothesis, AIC_diff, neutral:co_infect, factor_key=TRUE)%>%
  ggplot(aes(x=HHSregion,y=AIC_diff,color=hypothesis))+
  geom_point(size = 3)+
  facet_wrap(~Pathogen2,scales="free_y",nrow = 1)+
  theme_classic()+
  scale_x_continuous(breaks = seq(1, 11, by = 1))+
  scale_colour_brewer(palette = "Dark2",name = "Hypothesis", 
                      labels = c("no-interaction", "inhibition of co-infection", "cross-immunity","inhibition of co-ifection + cross-immunity"))+
  theme(legend.position = "bottom")+
  ylab("AIC difference compared with no-interaction hypothesis")->AIC_diff
 






## plot trajectory fit
inc_data_fit<-inc_data_perdic%>%
  mutate(HHSregion=HHS_REGION)%>%
  pivot_wider(names_from = virus, values_from = cases)%>%
  drop_na()%>%
  mutate(method=case_when(
    time<=6.5 ~ "Model fitting",
    time>6.5 ~ "Prediction"
  ))

traj_fit_arsv<-rbind(traj_fit_arsv_neutral,traj_fit_arsv_psi,traj_fit_arsv_chi,traj_fit_arsv_coinfect)
traj_fit_arsv_coinfect%>%
  mutate(hypothesis="co_infect")
hypo_levels=c("neutral","psi","chi","co_infect")
patho_order<-c("RSV","fluA")
traj_fit_arsv%>%
  mutate(hypothesis=factor(hypothesis,levels=hypo_levels))%>%
  mutate(pathogen=factor(type,levels=patho_order))%>%
  mutate(method=case_when(
    time<=6.5 ~ "Model fitting",
    time>6.5 ~ "Prediction"
  ))%>%
  ggplot(aes(x=time+2011,y=cases))+
  geom_area(data=inc_data_fit,aes(x=date,y=RSV),fill="brown",alpha=0.4)+
  
  geom_area(data=inc_data_fit,aes(x=date,y=fluA),fill="gray",alpha=0.6)+
  geom_line(aes(color=hypothesis,linetype=pathogen),size=0.5)+
  facet_grid(HHSregion~method,scales="free_x",space="free")+
  
  
  geom_area(data=inc_data_fit,aes(x=date,y=fluA),fill="gray",alpha=0.4)+
  geom_line(aes(color=hypothesis,linetype=type))+
  facet_grid(HHSregion~method,scales= "free", space="free_x")+
  scale_colour_brewer(palette = "Dark2")+
  theme_classic()+
  xlab("time(weeks)")+
  theme(legend.position="bottom")+
  scale_colour_brewer(palette = "Dark2",name = "Hypothesis", 
                      labels = c("no-interaction", "inhibition of co-infection", "cross-immunity","inhibition of co-ifection + cross-immunity")) ->arsv_trajfit
#+
# 

## plot trajectory fit
inc_data_fit<-inc_data_perdic%>%
  mutate(HHSregion=HHS_REGION)%>%
  pivot_wider(names_from = virus, values_from = cases)%>%
  drop_na()%>%
  mutate(method=case_when(
    time<=6.5 ~ "Fit",
    time>6.5 ~ "predict"
  ))

traj_fit_brsv<-rbind(traj_fit_brsv_neutral,traj_fit_brsv_coinfect,traj_fit_brsv_psi,traj_fit_brsv_chi)

traj_fit_brsv_coinfect%>%
  mutate(hypothesis="co_infect")
hypo_levels=c("neutral","psi","chi","co_infect")
patho_order<-c("RSV","fluB")

inc_data_fit<-inc_data_perdic%>%
  mutate(HHSregion=HHS_REGION)%>%
  pivot_wider(names_from = virus, values_from = cases)%>%
  drop_na()%>%
  mutate(method=case_when(
    time<=6.5 ~ "Model fitting",
    time>6.5 ~ "Prediction"
  ))


traj_fit_brsv%>%
  mutate(hypothesis=factor(hypothesis,levels=hypo_levels))%>%
  mutate(pathogen=factor(type,levels=patho_order))%>%
  mutate(method=case_when(
    time<=6.5 ~ "Model fitting",
    time>6.5 ~ "Prediction"
  ))%>%

  ggplot(aes(x=time+2011,y=cases))+
  geom_line(aes(color=hypothesis,linetype=pathogen ),alpha=0.5)+
  
  geom_area(data=inc_data_fit,aes(x=date,y=fluB),fill="gray",alpha=0.6)+
  geom_area(data=inc_data_fit,aes(x=date,y=RSV),fill="brown",alpha=0.4)+
  facet_grid(HHSregion~method,scales="free",space="free_x")+
  scale_colour_brewer(palette = "Dark2")+
  theme_classic()+
  xlab("time(weeks)")+
  theme(legend.position="bottom")+
  scale_colour_brewer(palette = "Dark2",name = "Hypothesis", 
                      labels = c("no-interaction", "inhibition of co-infection", "cross-immunity","inhibition of co-ifection + cross-immunity"))->brsv_trajfit


library(ggpubr)
ggarrange(
  arsv_trajfit, brsv_trajfit, labels = c("B", "C"),
  nrow=1,
  common.legend = TRUE)->p_trajfit

ggarrange(
  AIC_diff,p_trajfit,
  nrow=2,
  heights = c(1, 2.5)
  
)

## load csv
res_arsv_psi_CI<-read.csv("./res_hhs_CI/res_ares_psilog_CI.csv")


res_arsv_chi_CI<-read.csv("./res_hhs_CI/res_ares_chilog_CI.csv")


res_arsv_chi_CI%>%
  mutate(`parameter value`=chi,
         parameter="chi")%>%
  select(-chi)->res_arsv_chi_CI

res_arsv_psi_CI%>%
  mutate(`parameter value`=psi,
         parameter="psi")%>%
  select(-psi)->res_arsv_psi_CI

rbind(res_arsv_chi_CI,res_arsv_psi_CI)%>%
  ggplot(aes(x= 1- `parameter value`,y=loglik))+
  geom_point(size =0.5)+
  geom_line()+
  #geom_smooth(se=FALSE, col='black')+
  facet_grid(HHS_region~parameter,scales="free_y")+
  geom_hline(aes(yintercept=loglik_cutoff),linetype="dashed",color="blue")+
  geom_vline(aes(xintercept=1-upper),linetype="dashed",color="blue")+
  geom_vline(aes(xintercept =1- lower),linetype="dashed",color="blue")+
  theme_bw()+
  xlab("parameter value")->arsv_CI
  


res_brsv_psi_CI<-read.csv("./res_hhs_CI/res_bresv_psilog_CI.csv")

res_brsv_chi_CI<-read.csv("./res_hhs_CI/res_bresv_chilog_CI.csv")

res_brsv_chi_CI%>%
  mutate(`parameter value`=chi,
         parameter="chi")%>%
  select(-chi)->res_brsv_chi_CI

res_brsv_psi_CI%>%
  mutate(`parameter value`=psi,
         parameter="psi")%>%
  select(-psi)->res_brsv_psi_CI

rbind(res_brsv_chi_CI,res_brsv_psi_CI)%>%
  ggplot(aes(x=1-`parameter value`,y=loglik))+
  geom_point(size =0.5)+
  geom_line()+
  #geom_smooth(se=FALSE, col='black')+
  facet_grid(HHS_region~parameter,scales="free_y")+
  geom_hline(aes(yintercept=loglik_cutoff),linetype="dashed",color="blue")+
  geom_vline(aes(xintercept=1-upper),linetype="dashed",color="blue")+
  geom_vline(aes(xintercept =1- lower),linetype="dashed",color="blue")+
  theme_bw()+
  xlab("parameter value")->brsv_CI

library(ggpubr)  
ggarrange(
  arsv_CI, brsv_CI, labels = c("A", "B"),
  nrow=1)

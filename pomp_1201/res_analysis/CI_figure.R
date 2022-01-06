## load csv
res_brsv_psi_CI<-read.csv("./res_hhs_CI/res_bresv_psilog_CI.csv")
res_brsv_psi_CI%>%
  ggplot(aes(x=psi,y=loglik))+
  geom_point(size =1)+
  geom_line()+
  facet_wrap( ~HHS_region,scales="free_y")+
  geom_hline(aes(yintercept=loglik_cutoff))+
  theme_bw()
  
  
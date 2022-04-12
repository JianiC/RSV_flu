## loglik profile
## check step 1
source("./fit_functions.R", chdir = TRUE) 




loglik_CI_chi<-function(region,loglik){
  loglik%>%
    filter(HHS_region==region)->region_loglik
  approx<-approx(region_loglik$chi,region_loglik$loglik,n=200)
  df_approx<-data.frame(chi=approx$x,
                        loglik=approx$y)
  maxloglik<-max( df_approx$loglik)
  mle<-region_loglik$chi[which.max(region_loglik$loglik)]
  #mle<-df_approx$chi[df_approx$chi == max(df_approx$loglik)]
  cutoff <- maxloglik-qchisq(p=0.95,df=1)/2
  lower<-range(subset(df_approx,loglik>cutoff)$chi)[1]
  upper<-range(subset(df_approx,loglik>cutoff)$chi)[2]
  res<-data.frame("HHS_region"=region,
                  "maxloglik"=maxloglik,
                  "95cutoff"=cutoff,
                  "mle"=mle,
                  "CIlow"=lower,
                  "CIhigh"=upper
                  )
  return(res)
  
}
arsv_chi_loglik1<-read.csv("./res_loglik/res_arsv_chi_loglik.csv")
arsv_psi_loglik1<-read.csv("./res_loglik/res_arsv_psi_loglik.csv")


arsv_CI_loglik<-data.frame()
for (i in 1:10){
  hhs_arsv_chi_CI<-loglik_CI_chi(i,arsv_chi_loglik1)
  arsv_CI_loglik<-rbind(arsv_CI_loglik,hhs_arsv_chi_CI)
}

arsv_chi_loglik1%>%
  ggplot(aes(x= 1-chi,y=loglik))+
  geom_point()+
  geom_line()+
  #geom_smooth(se =F )+
  facet_grid(HHS_region ~ ., scales="free_y")+
  theme_bw()+
  geom_hline(data=arsv_CI_loglik,aes(yintercept = X95cutoff))+
  geom_vline(data=arsv_CI_loglik,aes(xintercept = 1-CIlow))+
  geom_vline(data=arsv_CI_loglik,aes(xintercept =1- CIhigh))->res_arsv_chi_CI


arsv_chi_loglik1%>%
  ggplot(aes(x= psi,y=loglik))+
  geom_point()+
  geom_line()+
  facet_grid(HHS_region ~ ., scales="free_y")+
  theme_bw()->res_arsv_psi_CI

arsv_chi_loglik<-arsv_chi_loglik1%>%
  rename(parm_value=chi)%>%
  mutate(hypo_value=1-parm_value)%>%
  mutate(parm="chi",
         hypo="Strengh of cross-immunity")
  
hypo_levels=c("Propotion of inhibition to be co-infected","Strengh of cross-immunity")

arsv_psi_loglik1%>%
  rename(parm_value=psi)%>%
  mutate(hypo_value=1-parm_value)%>%
  mutate(parm="psi",
         hypo="Propotion of inhibition to be co-infected")%>%
  rbind(arsv_chi_loglik)%>%
  mutate(hypo=factor(hypo,levels=hypo_levels))%>%
  ggplot(aes(x=hypo_value,y=loglik))+
  geom_point()+
  geom_line(aes(colour=hypo))+
  facet_grid(HHS_region~hypo,scales="free_y")+
  theme_bw()+
  scale_colour_brewer(palette = "Dark2")+
  xlab("")+
  theme(legend.position="none")->arsv_loglik


library(ggpubr)  
ggarrange(
  res_arsv_chi_CI, res_arsv_psi_CI, labels = c("chi", "psi"),
  nrow=1)




brsv_chi_loglik1<-read.csv("./res_loglik/res_brsv_chi_loglik.csv")

brsv_psi_loglik1<-read.csv("./res_loglik/res_brsv_psi_loglik.csv")

hypo_levels=c("Propotion of inhibition to be co-infected","Strengh of cross-immunity")

brsv_chi_loglik<-brsv_chi_loglik1%>%
  rename(parm_value=chi)%>%
  mutate(hypo_value=1-parm_value)%>%
  mutate(parm="chi",
         hypo="Strengh of cross-immunity")

brsv_psi_loglik1%>%
  rename(parm_value=psi)%>%
  mutate(hypo_value=1-parm_value)%>%
  mutate(parm="psi",
         hypo="Propotion of inhibition to be co-infected")%>%
  rbind(brsv_chi_loglik)%>%
  mutate(hypo=factor(hypo,levels=hypo_levels))%>%
  ggplot(aes(x=hypo_value,y=loglik))+
  geom_point()+
  geom_line(aes(colour=hypo))+
  facet_grid(HHS_region~hypo,scales="free_y")+
  theme_bw()+
  scale_colour_brewer(palette = "Dark2")+
  xlab("")+
  theme(legend.position="none")->brsv_loglik

ggarrange(
  arsv_loglik, brsv_loglik, labels = c("A", "B"),
  nrow=1)->p_loglik

## quick plot

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
                  "CIhigh"=upper,
                  "parm"="chi",
                  hypo="Strengh of cross-immunity"
  )
  return(res)
  
}

hypo_levels=c("Propotion of inhibition to be co-infected","Strengh of cross-immunity")
arsv_chi_loglik <-read.csv("arsv_best_loglik/arsv_best_chi_loglik.csv")
region<-c(1,2,3,4,5,6,7,9,10)
arsv_chi_CI_loglik<-data.frame()
for (i in region){
  hhs_arsv_chi_CI<-loglik_CI_chi(i,arsv_chi_loglik)
  arsv_chi_CI_loglik<-rbind(arsv_chi_CI_loglik,hhs_arsv_chi_CI)
}

arsv_psi_loglik <-read.csv("arsv_best_loglik/arsv_best_psi_loglik.csv")

arsv_chi_loglik<-arsv_chi_loglik%>%
  rename(parm_value=chi)%>%
  mutate(hypo_value=1-parm_value)%>%
  mutate(parm="chi",
         hypo="Strengh of cross-immunity")


arsv_psi_loglik%>%
  rename(parm_value=psi)%>%
  mutate(hypo_value=1-parm_value)%>%
  mutate(parm="psi",
         hypo="Propotion of inhibition to be co-infected")%>%
  rbind(arsv_chi_loglik)%>%
  mutate(hypo=factor(hypo,levels=hypo_levels))%>%
  ggplot(aes(x=hypo_value,y=loglik))+
  geom_point(aes(colour=parm))+
  geom_line(aes(colour=parm))+
  geom_hline(data=arsv_chi_CI_loglik,aes(yintercept=X95cutoff,colour=parm),linetype="dashed")+
  geom_vline(data=arsv_chi_CI_loglik,aes(xintercept=1-CIlow,colour=parm),linetype="dashed")+
  geom_vline(data=arsv_chi_CI_loglik,aes(xintercept=1-CIhigh,colour=parm),linetype="dashed")+
  facet_wrap(HHS_region ~ ., scales="free")+
  theme_classic()+
  xlab("Parameter value")+
  scale_colour_brewer(palette = "Dark2",name="Parameters",
                      labels=c("Strength of cross-immunity",
                                                "Propotion of inhibition of co-infection"))->arsv_loglik
  
  
  

###############################################

brsv_chi_loglik <-read.csv("brsv_best_loglik/brsv_loglik_chi.csv")
brsv_chi_CI_loglik<-data.frame()
for (i in 1:10){
  hhs_brsv_chi_CI<-loglik_CI_chi(i,brsv_chi_loglik)
  brsv_chi_CI_loglik<-rbind(brsv_chi_CI_loglik,hhs_brsv_chi_CI)
}

brsv_chi_loglik<-brsv_chi_loglik%>%
  rename(parm_value=chi)%>%
  mutate(hypo_value=1-parm_value)%>%
  mutate(parm="chi",
         hypo="Strengh of cross-immunity")


brsv_psi_loglik <-read.csv("brsv_best_loglik/brsv_loglik_psi.csv")
brsv_psi_loglik%>%
  rename(parm_value=psi)%>%
  mutate(hypo_value=1-parm_value)%>%
  mutate(parm="psi",
         hypo="Propotion of inhibition to be co-infected")%>%
  rbind(brsv_chi_loglik)%>%
  mutate(hypo=factor(hypo,levels=hypo_levels))%>%
  ggplot(aes(x=hypo_value,y=loglik))+
  geom_point(aes(colour=parm))+
  geom_line(aes(colour=parm))+
  geom_hline(data=brsv_chi_CI_loglik,aes(yintercept=X95cutoff,colour=parm),linetype="dashed")+
  geom_vline(data=brsv_chi_CI_loglik,aes(xintercept=1-CIlow,colour=parm),linetype="dashed")+
  geom_vline(data=brsv_chi_CI_loglik,aes(xintercept=1-CIhigh,colour=parm),linetype="dashed")+
  facet_wrap(HHS_region ~ ., scales="free")+
  theme_classic()+
  xlab("Parameter value")+
  scale_colour_brewer(palette = "Dark2",name="Parameters",
                      labels=c("Strength of cross-immunity",
                               "Propotion of inhibition of co-infection"))->brsv_loglik

ggarrange(
  arsv_loglik, brsv_loglik, labels = c("A", "B"),
  nrow=2,
  common.legend = TRUE
)->p_loglik

arsv_chi_CI_loglik%>%
  filter(HHS_region==1)->HHS1_arsv_chi_CI_loglik
  
arsv_psi_loglik%>%
  rename(parm_value=psi)%>%
  mutate(hypo_value=1-parm_value)%>%
  mutate(parm="psi",
         hypo="Propotion of inhibition to be co-infected")%>%
  rbind(arsv_chi_loglik)%>%
  mutate(hypo=factor(hypo,levels=hypo_levels))%>%
  filter(HHS_region==1)%>%
  ggplot(aes(x=hypo_value,y=loglik))+
  geom_point(aes(colour=parm))+
  geom_line(aes(colour=parm))+
  geom_hline(data=HHS1_arsv_chi_CI_loglik,aes(yintercept=X95cutoff,colour=parm),linetype="dashed")+
  geom_vline(data=HHS1_arsv_chi_CI_loglik,aes(xintercept=1-CIlow,colour=parm),linetype="dashed")+
  geom_vline(data=HHS1_arsv_chi_CI_loglik,aes(xintercept=1-CIhigh,colour=parm),linetype="dashed")+
  facet_wrap(parm ~ .)+
  theme_classic()+
  xlab("Parameter value")+
  ylab("loglikelihood")+
  scale_colour_brewer(palette = "Dark2",name="Parameters",
                      labels=c("Strength of cross-immunity",
                               "Propotion of inhibition of co-infection"))+
  gg.theme +
  theme(
        legend.position = "none") +
  guides(linetype = guide_legend(ncol = 1, order = 1), 
         colour = guide_legend(ncol = 1)) +
  theme(legend.spacing.y = unit(0.1, "lines"),
        legend.key = element_blank(), 
        legend.background = element_blank(),
        strip.text.x = element_blank())+
  scale_x_continuous(
    labels = label_number(accuracy = 0.1),
    breaks = seq(0,1,0.25))->HHS1_arsv_loglik

brsv_chi_CI_loglik%>%
  filter(HHS_region==1)->HHS1_brsv_chi_CI_loglik

brsv_psi_loglik%>%
  rename(parm_value=psi)%>%
  mutate(hypo_value=1-parm_value)%>%
  mutate(parm="psi",
         hypo="Propotion of inhibition to be co-infected")%>%
  rbind(brsv_chi_loglik)%>%
  mutate(hypo=factor(hypo,levels=hypo_levels))%>%
  filter(HHS_region==1)%>%
  drop_na()%>%
  ggplot(aes(x=hypo_value,y=loglik))+
  geom_point(aes(colour=parm))+
  geom_line(aes(colour=parm))+
  geom_hline(data=HHS1_brsv_chi_CI_loglik,aes(yintercept=X95cutoff,colour=parm),linetype="dashed")+
  geom_vline(data=HHS1_brsv_chi_CI_loglik,aes(xintercept=1-CIlow,colour=parm),linetype="dashed")+
  geom_vline(data=HHS1_brsv_chi_CI_loglik,aes(xintercept=1-CIhigh,colour=parm),linetype="dashed")+
  facet_wrap(parm ~ .)+
  theme_classic()+
  xlab("Parameter value")+
  ylab("Loglikelihood")+
  scale_colour_brewer(palette = "Dark2",name="Parameters",
                      labels=c("Strength of cross-immunity",
                               "Propotion of inhibition of co-infection"))+
  gg.theme +
  theme(
        legend.position = "none") +
  guides(linetype = guide_legend(ncol = 1, order = 1), 
         colour = guide_legend(ncol = 1)) +
  theme(legend.spacing.y = unit(0.1, "lines"),
        legend.key = element_blank(), 
        legend.background = element_blank(),
        strip.text.x = element_blank())+
  scale_x_continuous(
    labels = label_number(accuracy = 0.1),
    breaks = seq(0,1,0.25))->HHS1_brsv_loglik


legend <- get_legend(
  HHS1_arsv_loglik + 
    guides(color = guide_legend(nrow = 2),) +
    theme(legend.position = "bottom",
          legend.title=element_blank(),
          text = element_text(size = 8))
)


plot_grid( HHS1_arsv_loglik+theme(text = element_text(size = 8)), 
           HHS1_brsv_loglik+theme(text = element_text(size = 8)),
           nrow = 2,
        
           labels = "AUTO", label_size = 12)->HHS1_rsv_loglik

plot_grid(legend,
          HHS1_rsv_loglik,
          ncol=1,
          rel_heights = c(.1, 1))->loglik_CI_repre


ggarrange(
  HHS1_arsv_loglik, HHS1_brsv_loglik, labels = c("A", "B"),
  nrow=2,
  common.legend = TRUE
)->p_loglik_HHS1
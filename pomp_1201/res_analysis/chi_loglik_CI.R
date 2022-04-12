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

arsv_chi_loglik <-read.csv("arsv_best_loglik/arsv_best_chi_loglik.csv")
region<-c(1,2,3,4,5,6,7,9,10)
arsv_chi_CI_loglik<-data.frame()
for (i in region){
  hhs_arsv_chi_CI<-loglik_CI_chi(i,arsv_chi_loglik)
  arsv_chi_CI_loglik<-rbind(arsv_chi_CI_loglik,hhs_arsv_chi_CI)
}

arsv_chi_loglik%>%
  mutate(parmvalue=1-chi)%>%
  ggplot(aes(x=parmvalue,y=loglik))+
  geom_point(aes())+
  geom_line(aes())+
  geom_point(data=arsv_chi_CI_loglik,aes(x=1-mle,y=maxloglik, color=parm),shape=1,size=3)+
  geom_hline(data=arsv_chi_CI_loglik,aes(yintercept=X95cutoff,colour=parm),linetype="dashed")+
  geom_vline(data=arsv_chi_CI_loglik,aes(xintercept=1-CIlow,colour=parm),linetype="dashed")+
  geom_vline(data=arsv_chi_CI_loglik,aes(xintercept=1-CIhigh,colour=parm),linetype="dashed")+
  facet_wrap(HHS_region ~ ., scales="free")+
  theme_classic()+
  theme(legend.position = "none")+
  xlab("Strength of cross-immunity")+
  scale_colour_brewer(palette = "Dark2")->p_arsv_chi_CI_loglik

###############################################################################

brsv_chi_loglik <-read.csv("brsv_best_loglik/brsv_loglik_chi.csv")
region<-c(1,2,3,4,5,6,7,9,10)
brsv_chi_CI_loglik<-data.frame()
for (i in region){
  hhs_brsv_chi_CI<-loglik_CI_chi(i,brsv_chi_loglik)
  brsv_chi_CI_loglik<-rbind(brsv_chi_CI_loglik,hhs_brsv_chi_CI)
}

brsv_chi_loglik%>%
  mutate(parmvalue=1-chi)%>%
  ggplot(aes(x=parmvalue,y=loglik))+
  geom_point(aes())+
  geom_line(aes())+
  geom_point(data=brsv_chi_CI_loglik,aes(x=1-mle,y=maxloglik, color=parm),shape=1,size=3)+
  geom_hline(data=brsv_chi_CI_loglik,aes(yintercept=X95cutoff,colour=parm),linetype="dashed")+
  geom_vline(data=brsv_chi_CI_loglik,aes(xintercept=1-CIlow,colour=parm),linetype="dashed")+
  geom_vline(data=brsv_chi_CI_loglik,aes(xintercept=1-CIhigh,colour=parm),linetype="dashed")+
  facet_wrap(HHS_region ~ ., scales="free")+
  theme_classic()+
  xlab("Strength of cross-immunity")+
  theme(legend.position = "none")+
  scale_colour_brewer(palette = "Dark2")->p_brsv_chi_CI_loglik

ggarrange(
  p_arsv_chi_CI_loglik, p_brsv_chi_CI_loglik, labels = c("A", "B"),
  nrow=2,
  legend = "none"
  #common.legend = TRUE
  )->p_loglik


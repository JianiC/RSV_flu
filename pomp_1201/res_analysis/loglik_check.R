## get loglikelihood for a estimate
source("./par_CI_funs.R", chdir = TRUE) 

res_hhs1_arsv_psiloglik<-par_loglik(res_hhs1_arsv_coinfect,"psi")
res_hhs1_arsv_psilog_CI<-CI_loglik_psi(res_hhs1_arsv_psiloglik)
save(res_hhs1_arsv_psilog_CI,file = "./res_hhs_CI/res_hhs1_arsv_psilog_CI.rds")

res_hhs2_arsv_psiloglik<-par_loglik(res_hhs2_arsv_coinfect,"psi")
res_hhs2_arsv_psilog_CI<-CI_loglik_psi(res_hhs2_arsv_psiloglik)
save(res_hhs2_arsv_psilog_CI,file = "./res_hhs_CI/res_hhs2_arsv_psilog_CI.rds")

res_hhs3_arsv_psiloglik<-par_loglik(res_hhs3_arsv_coinfect,"psi")
res_hhs3_arsv_psilog_CI<-CI_loglik_psi(res_hhs3_arsv_psiloglik)
save(res_hhs3_arsv_psilog_CI,file = "./res_hhs_CI/res_hhs3_arsv_psilog_CI.rds")

res_hhs4_arsv_psiloglik<-par_loglik(res_hhs4_arsv_coinfect,"psi")
res_hhs4_arsv_psilog_CI<-CI_loglik_psi(res_hhs4_arsv_psiloglik)
save(res_hhs4_arsv_psilog_CI,file = "./res_hhs_CI/res_hhs4_arsv_psilog_CI.rds")

res_hhs5_arsv_psiloglik<-par_loglik(res_hhs5_arsv_coinfect,"psi")
res_hhs5_arsv_psilog_CI<-CI_loglik_psi(res_hhs5_arsv_psiloglik)
save(res_hhs5_arsv_psilog_CI,file = "./res_hhs_CI/res_hhs5_arsv_psilog_CI.rds")

res_hhs6_arsv_psiloglik<-par_loglik(res_hhs1_arsv_coinfect,"psi")
res_hhs6_arsv_psilog_CI<-CI_loglik_psi(res_hhs6_arsv_psiloglik)
save(res_hhs6_arsv_psilog_CI,file = "./res_hhs_CI/res_hhs6_arsv_psilog_CI.rds")

res_hhs7_arsv_psiloglik<-par_loglik(res_hhs7_arsv_coinfect,"psi")
res_hhs7_arsv_psilog_CI<-CI_loglik_psi(res_hhs7_arsv_psiloglik)
save(res_hhs7_arsv_psilog_CI,file = "./res_hhs_CI/res_hhs7_arsv_psilog_CI.rds")

res_hhs8_arsv_psiloglik<-par_loglik(res_hhs8_arsv_coinfect,"psi")
res_hhs8_arsv_psilog_CI<-CI_loglik_psi(res_hhs8_arsv_psiloglik)
save(res_hhs8_arsv_psilog_CI,file = "./res_hhs_CI/res_hhs8_arsv_psilog_CI.rds")

res_hhs9_arsv_psiloglik<-par_loglik(res_hhs9_arsv_coinfect,"psi")
res_hhs9_arsv_psilog_CI<-CI_loglik_psi(res_hhs9_arsv_psiloglik)
save(res_hhs9_arsv_psilog_CI,file = "./res_hhs_CI/res_hhs9_arsv_psilog_CI.rds")

res_hhs10_arsv_psiloglik<-par_loglik(res_hhs10_arsv_coinfect,"psi")
res_hhs10_arsv_psilog_CI<-CI_loglik_psi(res_hhs10_arsv_psiloglik)
save(res_hhs10_arsv_psilog_CI,file = "./res_hhs_CI/res_hhs10_arsv_psilog_CI.rds")


##############################################################################

## loglik profile for RSV and fluA  chi estiamtes

##############################################################################
res_hhs1_arsv_chiloglik<-par_loglik(res_hhs1_arsv_coinfect,"chi")
res_hhs1_arsv_chilog_CI<-CI_loglik_chi(res_hhs1_arsv_chiloglik)
save(res_hhs1_arsv_chilog_CI,file = "./res_hhs_CI/res_hhs1_arsv_chilog_CI.rds")

res_hhs2_arsv_chiloglik<-par_loglik(res_hhs2_arsv_coinfect,"chi")
res_hhs2_arsv_chilog_CI<-CI_loglik_chi(res_hhs2_arsv_chiloglik)
save(res_hhs2_arsv_chilog_CI,file = "./res_hhs_CI/res_hhs2_arsv_chilog_CI.rds")

res_hhs3_arsv_chiloglik<-par_loglik(res_hhs3_arsv_coinfect,"chi")
res_hhs3_arsv_chilog_CI<-CI_loglik_chi(res_hhs3_arsv_chiloglik)
save(res_hhs3_arsv_chilog_CI,file = "./res_hhs_CI/res_hhs3_arsv_chilog_CI.rds")

res_hhs4_arsv_chiloglik<-par_loglik(res_hhs4_arsv_coinfect,"chi")
res_hhs4_arsv_chilog_CI<-CI_loglik_chi(res_hhs4_arsv_chiloglik)
save(res_hhs4_arsv_chilog_CI,file = "./res_hhs_CI/res_hhs4_arsv_chilog_CI.rds")

res_hhs5_arsv_chiloglik<-par_loglik(res_hhs5_arsv_coinfect,"chi")
res_hhs5_arsv_chilog_CI<-CI_loglik_chi(res_hhs5_arsv_chiloglik)
save(res_hhs5_arsv_chilog_CI,file = "./res_hhs_CI/res_hhs5_arsv_chilog_CI.rds")

res_hhs6_arsv_chiloglik<-par_loglik(res_hhs1_arsv_coinfect,"chi")
res_hhs6_arsv_chilog_CI<-CI_loglik_chi(res_hhs6_arsv_chiloglik)
save(res_hhs6_arsv_chilog_CI,file = "./res_hhs_CI/res_hhs6_arsv_chilog_CI.rds")

res_hhs7_arsv_chiloglik<-par_loglik(res_hhs7_arsv_coinfect,"chi")
res_hhs7_arsv_chilog_CI<-CI_loglik_chi(res_hhs7_arsv_chiloglik)
save(res_hhs7_arsv_chilog_CI,file = "./res_hhs_CI/res_hhs7_arsv_chilog_CI.rds")

res_hhs8_arsv_chiloglik<-par_loglik(res_hhs8_arsv_coinfect,"chi")
res_hhs8_arsv_chilog_CI<-CI_loglik_chi(res_hhs8_arsv_chiloglik)
save(res_hhs8_arsv_chilog_CI,file = "./res_hhs_CI/res_hhs8_arsv_chilog_CI.rds")

res_hhs9_arsv_chiloglik<-par_loglik(res_hhs9_arsv_coinfect,"chi")
res_hhs9_arsv_chilog_CI<-CI_loglik_chi(res_hhs9_arsv_chiloglik)
save(res_hhs9_arsv_chilog_CI,file = "./res_hhs_CI/res_hhs9_arsv_chilog_CI.rds")

res_hhs10_arsv_chiloglik<-par_loglik(res_hhs10_arsv_coinfect,"chi")
res_hhs10_arsv_chilog_CI<-CI_loglik_chi(res_hhs10_arsv_chiloglik)
save(res_hhs10_arsv_chilog_CI,file = "./res_hhs_CI/res_hhs10_arsv_chilog_CI.rds")

######################################################################################
## Check the loglik profile functions - res_hhs3_arsv "chi"
######################################################################################
res_hhs3_arsv_coinfect$total2

# 1. loglik profile function
res_arsv_hhs3_pomp_data_hhs <- (
  inc_data_add %.>% 
    make_data_pomp_ready(., virus_combo =c("RSV","fluA"), HHS_region = res_hhs3_arsv_coinfect$HHS)
)


res_arsv_hhs3_fit_par <- get_rp_vals(data = inc_data_add,res_hhs3_arsv_coinfect)
loglik(fit_par = res_arsv_hhs3_fit_par,pomp_data = res_arsv_hhs3_pomp_data_hhs )
res_hhs3_arsv_coinfect$DEobj$optim$bestval
#######################################################################################

## 2. check simulation process
result<-data.frame()
for (i in seq(0.3, 0.5, by = 0.001)) {
  
  res_arsv_hhs3_fit_par ["chi"]<-i
  loglik_value <-loglik(fit_par = res_arsv_hhs3_fit_par,pomp_data = res_arsv_hhs3_pomp_data_hhs)
  
  out<-list("chi"=i,"loglik"=-loglik_value)
  print(out)
  result<-rbind(result,out)
  
}

#######################################################################################
maxloglik <- max(result$loglik)
plot(loglik~chi,result,type="l")
cutoff <- maxloglik-qchisq(p=0.95,df=1)/2
range(subset(result,loglik>cutoff))

## narrow down the simulation space of chi

result2<-data.frame()
for (i in seq(0.4135, 0.4145, by = 0.00005)) {
  
  res_arsv_hhs3_fit_par ["chi"]<-i
  loglik_value <-loglik(fit_par = res_arsv_hhs3_fit_par,pomp_data = res_arsv_hhs3_pomp_data_hhs)
  
  out<-list("chi"=i,"loglik"=-loglik_value)
  print(out)
  result2<-rbind(result2,out)
  
}

maxloglik <- max(result2$loglik)
plot(loglik~chi,result2,type="l")
cutoff <- maxloglik-qchisq(p=0.95,df=1)/2
range(subset(result2,loglik>cutoff)$chi)


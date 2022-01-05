## get loglikelihood for a estimate
source("./fit_functions.R", chdir = TRUE) 
source("./par_CI_funs.R", chdir = TRUE) 

#open $HOME/.Renviron
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

res_hhs6_arsv_psiloglik<-par_loglik(res_hhs6_arsv_coinfect,"psi")
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
###
chi
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

res_hhs6_arsv_chiloglik<-par_loglik(res_hhs6_arsv_coinfect,"chi")
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













maxloglik <- max(test$loglik)
plot(loglik~psi,test,type="l",ylim=maxloglik+c(-10,0))
cutoff <- maxloglik-qchisq(p=0.95,df=1)/2
abline(h=c(0,cutoff))
abline(v=range(subset(test,loglik>cutoff)$psi),lty=2)

range(subset(test,loglik>cutoff)$psi)[1]

loglik(fit_par,pomp_data = pomp_data_hhs1_arsv )

fit_par<-get_rp_vals(data=inc_data_add,res_hhs1_arsv_coinfect)

## get loglikelihood for a estimate
source("./fit_functions.R", chdir = TRUE) 
source("./par_CI_funs.R", chdir = TRUE) 

#open $HOME/.Renviron

res_hhs1_arsv_psiloglik<-par_loglik(res_hhs1_arsv_coinfect,"psi")
res_hhs1_arsv_psilog_CI<-CI_loglik_psi(res_hhs1_arsv_psiloglik)
write.csv(res_hhs1_arsv_psilog_CI,"./res_hhs_CI/res_hhs1_arsv_psilog_CI.csv",row.names = FALSE)

res_hhs2_arsv_psiloglik<-par_loglik(res_hhs2_arsv_coinfect,"psi")
res_hhs2_arsv_psilog_CI<-CI_loglik_psi(res_hhs2_arsv_psiloglik)
write.csv(res_hhs2_arsv_psilog_CI,"./res_hhs_CI/res_hhs2_arsv_psilog_CI.csv",row.names = FALSE)

res_hhs3_arsv_psiloglik<-par_loglik(res_hhs3_arsv_coinfect,"psi")
res_hhs3_arsv_psilog_CI<-CI_loglik_psi(res_hhs3_arsv_psiloglik)
write.csv(res_hhs3_arsv_psilog_CI,"./res_hhs_CI/res_hhs3_arsv_psilog_CI.csv",row.names = FALSE)


res_hhs4_arsv_psiloglik<-par_loglik(res_hhs4_arsv_coinfect,"psi")
res_hhs4_arsv_psilog_CI<-CI_loglik_psi(res_hhs4_arsv_psiloglik)
write.csv(res_hhs4_arsv_psilog_CI,"./res_hhs_CI/res_hhs4_arsv_psilog_CI.csv",row.names = FALSE)

res_hhs5_arsv_psiloglik<-par_loglik(res_hhs5_arsv_coinfect,"psi")
res_hhs5_arsv_psilog_CI<-CI_loglik_psi(res_hhs5_arsv_psiloglik)
write.csv(res_hhs5_arsv_psilog_CI,"./res_hhs_CI/res_hhs5_arsv_psilog_CI.csv",row.names = FALSE)

res_hhs6_arsv_psiloglik<-par_loglik(res_hhs6_arsv_coinfect,"psi")
res_hhs6_arsv_psilog_CI<-CI_loglik_psi(res_hhs6_arsv_psiloglik)
write.csv(res_hhs6_arsv_psilog_CI,"./res_hhs_CI/res_hhs6_arsv_psilog_CI.csv",row.names = FALSE)

res_hhs7_arsv_psiloglik<-par_loglik(res_hhs7_arsv_coinfect,"psi")
res_hhs7_arsv_psilog_CI<-CI_loglik_psi(res_hhs7_arsv_psiloglik)
write.csv(res_hhs7_arsv_psilog_CI,"./res_hhs_CI/res_hhs7_arsv_psilog_CI.csv",row.names = FALSE)

res_hhs8_arsv_psiloglik<-par_loglik(res_hhs8_arsv_coinfect,"psi")
res_hhs8_arsv_psilog_CI<-CI_loglik_psi(res_hhs8_arsv_psiloglik)
write.csv(res_hhs8_arsv_psilog_CI,"./res_hhs_CI/res_hhs8_arsv_psilog_CI.csv",row.names = FALSE)

res_hhs9_arsv_psiloglik<-par_loglik(res_hhs9_arsv_coinfect,"psi")
res_hhs9_arsv_psilog_CI<-CI_loglik_psi(res_hhs9_arsv_psiloglik)
write.csv(res_hhs9_arsv_psilog_CI,"./res_hhs_CI/res_hhs9_arsv_psilog_CI.csv",row.names = FALSE)

res_hhs10_arsv_psiloglik<-par_loglik(res_hhs10_arsv_coinfect,"psi")
res_hhs10_arsv_psilog_CI<-CI_loglik_psi(res_hhs10_arsv_psiloglik)
write.csv(res_hhs10_arsv_psilog_CI,"./res_hhs_CI/res_hhs10_arsv_psilog_CI.csv",row.names = FALSE)
###

res_hhs1_arsv_chiloglik<-par_loglik(res_hhs1_arsv_coinfect,"chi")
res_hhs1_arsv_chilog_CI<-CI_loglik_chi(res_hhs1_arsv_chiloglik)
write.csv(res_hhs1_arsv_chilog_CI,"./res_hhs_CI/res_hhs1_arsv_chilog_CI.csv",row.names = FALSE)

res_hhs2_arsv_chiloglik<-par_loglik(res_hhs2_arsv_coinfect,"chi")
res_hhs2_arsv_chilog_CI<-CI_loglik_chi(res_hhs2_arsv_chiloglik)
write.csv(res_hhs2_arsv_chilog_CI,"./res_hhs_CI/res_hhs2_arsv_chilog_CI.csv",row.names = FALSE)

res_hhs3_arsv_chiloglik<-par_loglik(res_hhs3_arsv_coinfect,"chi")
res_hhs3_arsv_chilog_CI<-CI_loglik_chi(res_hhs3_arsv_chiloglik)
write.csv(res_hhs3_arsv_chilog_CI,"./res_hhs_CI/res_hhs3_arsv_chilog_CI.csv",row.names = FALSE)


res_hhs4_arsv_chiloglik<-par_loglik(res_hhs4_arsv_coinfect,"chi")
res_hhs4_arsv_chilog_CI<-CI_loglik_chi(res_hhs4_arsv_chiloglik)
write.csv(res_hhs4_arsv_chilog_CI,"./res_hhs_CI/res_hhs4_arsv_chilog_CI.csv",row.names = FALSE)

res_hhs5_arsv_chiloglik<-par_loglik(res_hhs5_arsv_coinfect,"chi")
res_hhs5_arsv_chilog_CI<-CI_loglik_chi(res_hhs5_arsv_chiloglik)
write.csv(res_hhs5_arsv_chilog_CI,"./res_hhs_CI/res_hhs5_arsv_chilog_CI.csv",row.names = FALSE)

res_hhs6_arsv_chiloglik<-par_loglik(res_hhs6_arsv_coinfect,"chi")
res_hhs6_arsv_chilog_CI<-CI_loglik_chi(res_hhs6_arsv_chiloglik)
write.csv(res_hhs6_arsv_chilog_CI,"./res_hhs_CI/res_hhs6_arsv_chilog_CI.csv",row.names = FALSE)

res_hhs7_arsv_chiloglik<-par_loglik(res_hhs7_arsv_coinfect,"chi")
res_hhs7_arsv_chilog_CI<-CI_loglik_chi(res_hhs7_arsv_chiloglik)
write.csv(res_hhs7_arsv_chilog_CI,"./res_hhs_CI/res_hhs7_arsv_chilog_CI.csv",row.names = FALSE)

res_hhs8_arsv_chiloglik<-par_loglik(res_hhs8_arsv_coinfect,"chi")
res_hhs8_arsv_chilog_CI<-CI_loglik_chi(res_hhs8_arsv_chiloglik)
write.csv(res_hhs8_arsv_chilog_CI,"./res_hhs_CI/res_hhs8_arsv_chilog_CI.csv",row.names = FALSE)

res_hhs9_arsv_chiloglik<-par_loglik(res_hhs9_arsv_coinfect,"chi")
res_hhs9_arsv_chilog_CI<-CI_loglik_chi(res_hhs9_arsv_chiloglik)
write.csv(res_hhs9_arsv_chilog_CI,"./res_hhs_CI/res_hhs9_arsv_chilog_CI.csv",row.names = FALSE)

res_hhs10_arsv_chiloglik<-par_loglik(res_hhs10_arsv_coinfect,"chi")
res_hhs10_arsv_chilog_CI<-CI_loglik_chi(res_hhs10_arsv_chiloglik)
write.csv(res_hhs10_arsv_chilog_CI,"./res_hhs_CI/res_hhs10_arsv_chilog_CI.csv",row.names = FALSE)

################interaction with fluB


res_hhs1_brsv_psiloglik<-par_loglik(res_hhs1_brsv_coinfect,"psi")
res_hhs1_brsv_psilog_CI<-CI_loglik_psi(res_hhs1_brsv_psiloglik)
write.csv(res_hhs1_brsv_psilog_CI,"./res_hhs_CI/res_hhs2_brsv_psilog_CI.csv",row.names = FALSE)

res_hhs2_brsv_psiloglik<-par_loglik(res_hhs2_brsv_coinfect,"psi")
res_hhs2_brsv_psilog_CI<-CI_loglik_psi(res_hhs2_brsv_psiloglik)
write.csv(res_hhs2_brsv_psilog_CI,"./res_hhs_CI/res_hhs2_brsv_psilog_CI.csv",row.names = FALSE)

res_hhs3_brsv_psiloglik<-par_loglik(res_hhs3_brsv_coinfect,"psi")
res_hhs3_brsv_psilog_CI<-CI_loglik_psi(res_hhs3_brsv_psiloglik)
write.csv(res_hhs3_brsv_psilog_CI,"./res_hhs_CI/res_hhs3_brsv_psilog_CI.csv",row.names = FALSE)


res_hhs4_brsv_psiloglik<-par_loglik(res_hhs4_brsv_coinfect,"psi")
res_hhs4_brsv_psilog_CI<-CI_loglik_psi(res_hhs4_brsv_psiloglik)
write.csv(res_hhs4_brsv_psilog_CI,"./res_hhs_CI/res_hhs4_brsv_psilog_CI.csv",row.names = FALSE)

res_hhs5_brsv_psiloglik<-par_loglik(res_hhs5_brsv_coinfect,"psi")
res_hhs5_brsv_psilog_CI<-CI_loglik_psi(res_hhs5_brsv_psiloglik)
write.csv(res_hhs5_brsv_psilog_CI,"./res_hhs_CI/res_hhs5_brsv_psilog_CI.csv",row.names = FALSE)

res_hhs6_brsv_psiloglik<-par_loglik(res_hhs6_brsv_coinfect,"psi")
res_hhs6_brsv_psilog_CI<-CI_loglik_psi(res_hhs6_brsv_psiloglik)
write.csv(res_hhs6_brsv_psilog_CI,"./res_hhs_CI/res_hhs6_brsv_psilog_CI.csv",row.names = FALSE)

res_hhs7_brsv_psiloglik<-par_loglik(res_hhs7_brsv_coinfect,"psi")
res_hhs7_brsv_psilog_CI<-CI_loglik_psi(res_hhs7_brsv_psiloglik)
write.csv(res_hhs7_brsv_psilog_CI,"./res_hhs_CI/res_hhs7_brsv_psilog_CI.csv",row.names = FALSE)

res_hhs8_brsv_psiloglik<-par_loglik(res_hhs8_brsv_coinfect,"psi")
res_hhs8_brsv_psilog_CI<-CI_loglik_psi(res_hhs8_brsv_psiloglik)
write.csv(res_hhs8_brsv_psilog_CI,"./res_hhs_CI/res_hhs8_brsv_psilog_CI.csv",row.names = FALSE)

res_hhs9_brsv_psiloglik<-par_loglik(res_hhs9_brsv_coinfect,"psi")
res_hhs9_brsv_psilog_CI<-CI_loglik_psi(res_hhs9_brsv_psiloglik)
write.csv(res_hhs9_brsv_psilog_CI,"./res_hhs_CI/res_hhs9_brsv_psilog_CI.csv",row.names = FALSE)

res_hhs10_brsv_psiloglik<-par_loglik(res_hhs10_brsv_coinfect,"psi")
res_hhs10_brsv_psilog_CI<-CI_loglik_psi(res_hhs10_brsv_psiloglik)
write.csv(res_hhs10_brsv_psilog_CI,"./res_hhs_CI/res_hhs10_brsv_psilog_CI.csv",row.names = FALSE)
###

res_hhs1_brsv_chiloglik<-par_loglik(res_hhs1_brsv_coinfect,"chi")
res_hhs1_brsv_chilog_CI<-CI_loglik_chi(res_hhs1_brsv_chiloglik)
write.csv(res_hhs1_brsv_chilog_CI,"./res_hhs_CI/res_hhs1_brsv_chilog_CI.csv",row.names = FALSE)

res_hhs2_brsv_chiloglik<-par_loglik(res_hhs2_brsv_coinfect,"chi")
res_hhs2_brsv_chilog_CI<-CI_loglik_chi(res_hhs2_brsv_chiloglik)
write.csv(res_hhs2_brsv_chilog_CI,"./res_hhs_CI/res_hhs2_brsv_chilog_CI.csv",row.names = FALSE)

res_hhs3_brsv_chiloglik<-par_loglik(res_hhs3_brsv_coinfect,"chi")
res_hhs3_brsv_chilog_CI<-CI_loglik_chi(res_hhs3_brsv_chiloglik)
write.csv(res_hhs3_brsv_chilog_CI,"./res_hhs_CI/res_hhs3_brsv_chilog_CI.csv",row.names = FALSE)


res_hhs4_brsv_chiloglik<-par_loglik(res_hhs4_brsv_coinfect,"chi")
res_hhs4_brsv_chilog_CI<-CI_loglik_chi(res_hhs4_brsv_chiloglik)
write.csv(res_hhs4_brsv_chilog_CI,"./res_hhs_CI/res_hhs4_brsv_chilog_CI.csv",row.names = FALSE)

res_hhs5_brsv_chiloglik<-par_loglik(res_hhs5_brsv_coinfect,"chi")
res_hhs5_brsv_chilog_CI<-CI_loglik_chi(res_hhs5_brsv_chiloglik)
write.csv(res_hhs5_brsv_chilog_CI,"./res_hhs_CI/res_hhs5_brsv_chilog_CI.csv",row.names = FALSE)

res_hhs6_brsv_chiloglik<-par_loglik(res_hhs6_brsv_coinfect,"chi")
res_hhs6_brsv_chilog_CI<-CI_loglik_chi(res_hhs6_brsv_chiloglik)
write.csv(res_hhs6_brsv_chilog_CI,"./res_hhs_CI/res_hhs6_brsv_chilog_CI.csv",row.names = FALSE)

res_hhs7_brsv_chiloglik<-par_loglik(res_hhs7_brsv_coinfect,"chi")
res_hhs7_brsv_chilog_CI<-CI_loglik_chi(res_hhs7_brsv_chiloglik)
write.csv(res_hhs7_brsv_chilog_CI,"./res_hhs_CI/res_hhs7_brsv_chilog_CI.csv",row.names = FALSE)

res_hhs8_brsv_chiloglik<-par_loglik(res_hhs8_brsv_coinfect,"chi")
res_hhs8_brsv_chilog_CI<-CI_loglik_chi(res_hhs8_brsv_chiloglik)
write.csv(res_hhs8_brsv_chilog_CI,"./res_hhs_CI/res_hhs8_brsv_chilog_CI.csv",row.names = FALSE)

res_hhs9_brsv_chiloglik<-par_loglik(res_hhs9_brsv_coinfect,"chi")
res_hhs9_brsv_chilog_CI<-CI_loglik_chi(res_hhs9_brsv_chiloglik)
write.csv(res_hhs9_brsv_chilog_CI,"./res_hhs_CI/res_hhs9_brsv_chilog_CI.csv",row.names = FALSE)

res_hhs10_brsv_chiloglik<-par_loglik(res_hhs10_brsv_coinfect,"chi")
res_hhs10_brsv_chilog_CI<-CI_loglik_chi(res_hhs10_brsv_chiloglik)
write.csv(res_hhs10_brsv_chilog_CI,"./res_hhs_CI/res_hhs10_brsv_chilog_CI.csv",row.names = FALSE)
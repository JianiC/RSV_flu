## check  mle
## FluA HHS3

res_arsv_hhs3_pomp_data_hhs <- (
  inc_data_add %.>% 
    make_data_pomp_ready(., virus_combo =c("RSV","fluA"), HHS_region = res_hhs3_arsv_coinfect$HHS)
)

## optimal estimate 
res_arsv_hhs3_fit_par <- get_rp_vals(data = inc_data_add,res_hhs3_arsv_coinfect)
res_arsv_hhs3_fit_par

loglik(fit_par = res_arsv_hhs3_fit_par,pomp_data = res_arsv_hhs3_pomp_data_hhs )->fit_loglik
fit_loglik
result<-data.frame()
for (i in seq(0, 1, by = 0.1)) {
  
  res_arsv_hhs3_fit_par ["psi"]<-i
  loglik_value <-loglik(fit_par = res_arsv_hhs3_fit_par,pomp_data = res_arsv_hhs3_pomp_data_hhs)
  
  out<-list("psi"=i,"loglik"=-loglik_value)
  #print(out)
  result<-rbind(result,out)
  
}
result

## check hhs6
res_arsv_hhs6_pomp_data_hhs <- (
  inc_data_add %.>% 
    make_data_pomp_ready(., virus_combo =c("RSV","fluA"), HHS_region = res_hhs6_arsv_coinfect$HHS)
)

## optimal estimate 
res_arsv_hhs6_fit_par <- get_rp_vals(data = inc_data_add,res_hhs6_arsv_coinfect)
res_arsv_hhs6_fit_par
res_hhs6_arsv_coinfect$DEobj$optim$bestval

loglik(fit_par = res_arsv_hhs6_fit_par,pomp_data = res_arsv_hhs6_pomp_data_hhs )->fit_loglik
fit_loglik
result<-data.frame()
for (i in seq(0, 1, by = 0.1)) {
  
  res_arsv_hhs6_fit_par ["psi"]<-i
  loglik_value <-loglik(fit_par = res_arsv_hhs6_fit_par,pomp_data = res_arsv_hhs6_pomp_data_hhs)
  
  out<-list("psi"=i,"loglik"=-loglik_value)
  #print(out)
  result<-rbind(result,out)
  
}
result
##############
## check HHS7
res_arsv_hhs7_pomp_data_hhs <- (
  inc_data_add %.>% 
    make_data_pomp_ready(., virus_combo =c("RSV","fluA"), HHS_region = res_hhs7_arsv_coinfect$HHS)
)

## optimal estimate 
res_arsv_hhs7_fit_par <- get_rp_vals(data = inc_data_add,res_hhs7_arsv_coinfect)
res_arsv_hhs7_fit_par
res_hhs7_arsv_coinfect$DEobj$optim$bestval

loglik(fit_par = res_arsv_hhs7_fit_par,pomp_data = res_arsv_hhs7_pomp_data_hhs )->fit_loglik
fit_loglik
result<-data.frame()
for (i in seq(0, 1, by = 0.1)) {
  
  res_arsv_hhs7_fit_par ["psi"]<-i
  loglik_value <-loglik(fit_par = res_arsv_hhs7_fit_par,pomp_data = res_arsv_hhs7_pomp_data_hhs)
  
  out<-list("psi"=i,"loglik"=-loglik_value)
  #print(out)
  result<-rbind(result,out)
  
}



result
###

res_arsv_hhs4_pomp_data_hhs <- (
  inc_data_add %.>% 
    make_data_pomp_ready(., virus_combo =c("RSV","fluA"), HHS_region = res_hhs4_arsv_coinfect$HHS)
)

## optimal estimate 
res_arsv_hhs4_fit_par <- get_rp_vals(data = inc_data_add,res_hhs4_arsv_coinfect)
res_arsv_hhs4_fit_par
res_hhs4_arsv_coinfect$DEobj$optim$bestval

loglik(fit_par = res_arsv_hhs4_fit_par,pomp_data = res_arsv_hhs4_pomp_data_hhs )->fit_loglik
fit_loglik
result<-data.frame()
for (i in seq(0, 1, by = 0.1)) {
  
  res_arsv_hhs4_fit_par ["psi"]<-i
  loglik_value <-loglik(fit_par = res_arsv_hhs4_fit_par,pomp_data = res_arsv_hhs4_pomp_data_hhs)
  
  out<-list("psi"=i,"loglik"=-loglik_value)
  #print(out)
  result<-rbind(result,out)
  
}
## HHS4

res_arsv_hhs4_pomp_data_hhs <- (
  inc_data_add %.>% 
    make_data_pomp_ready(., virus_combo =c("RSV","fluA"), HHS_region = res_hhs4_arsv_chi$HHS)
)

## optimal estimate 
res_arsv_hhs4_fit_par <- get_rp_vals(data = inc_data_add,res_hhs4_arsv_chi)
res_arsv_hhs4_fit_par
res_hhs4_arsv_chi$DEobj$optim$bestval

loglik(fit_par = res_arsv_hhs4_fit_par,pomp_data = res_arsv_hhs4_pomp_data_hhs )->fit_loglik
fit_loglik
result<-data.frame()
for (i in seq(0, 1, by = 0.1)) {
  
  res_arsv_hhs4_fit_par ["psi"]<-i
  loglik_value <-loglik(fit_par = res_arsv_hhs4_fit_par,pomp_data = res_arsv_hhs4_pomp_data_hhs)
  
  out<-list("psi"=i,"loglik"=-loglik_value)
  #print(out)
  result<-rbind(result,out)
  
}

res_arsv_hhs5_pomp_data_hhs <- (
  inc_data_add %.>% 
    make_data_pomp_ready(., virus_combo =c("RSV","fluA"), HHS_region = res_hhs5_arsv_coinfect$HHS)
)

## optimal estimate 
res_arsv_hhs5_fit_par <- get_rp_vals(data = inc_data_add,res_hhs5_arsv_coinfect)
res_arsv_hhs5_fit_par
res_hhs5_arsv_coinfect$DEobj$optim$bestval
res_hhs5_arsv_chi$DEobj$optim$bestval
loglik(fit_par = res_arsv_hhs5_fit_par,pomp_data = res_arsv_hhs5_pomp_data_hhs )->fit_loglik
fit_loglik
result<-data.frame()
for (i in seq(0, 1, by = 0.1)) {
  
  res_arsv_hhs5_fit_par ["psi"]<-i
  loglik_value <-loglik(fit_par = res_arsv_hhs5_fit_par,pomp_data = res_arsv_hhs5_pomp_data_hhs)
  
  out<-list("psi"=i,"loglik"=-loglik_value)
  #print(out)
  result<-rbind(result,out)
  
}


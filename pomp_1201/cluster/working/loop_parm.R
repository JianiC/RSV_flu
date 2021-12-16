## get parameter estimates
getparam(res_hhs=res_hhs6_arsv_coinfect)
## write a loop to rbind all esitimates

res_arsv_coinfect_list<-list(
               res_hhs1_arsv_coinfect,res_hhs2_arsv_coinfect,res_hhs3_arsv_coinfect,
               res_hhs4_arsv_coinfect,res_hhs5_arsv_coinfect,res_hhs6_arsv_coinfect,
               res_hhs7_arsv_coinfect,res_hhs8_arsv_coinfect,res_hhs9_arsv_coinfect,res_hhs10_arsv_coinfect)
est_res_arsv_coinfect<-get_est_all(res_arsv_coinfect_list)
rm(res_arsv_coinfect_list,res_hhs1_arsv_coinfect,res_hhs2_arsv_coinfect,res_hhs3_arsv_coinfect,
      res_hhs4_arsv_coinfect,res_hhs5_arsv_coinfect,res_hhs6_arsv_coinfect,
      res_hhs7_arsv_coinfect,res_hhs8_arsv_coinfect,res_hhs9_arsv_coinfect,res_hhs10_arsv_coinfect)
######
res_arsv_neutral_list<-list(res_hhs1_arsv_neutral,res_hhs2_arsv_neutral,res_hhs3_arsv_neutral,
                            res_hhs4_arsv_neutral,res_hhs5_arsv_neutral,res_hhs6_arsv_neutral,
                            res_hhs7_arsv_neutral,res_hhs8_arsv_neutral,res_hhs9_arsv_neutral,res_hhs10_arsv_neutral)
est_res_arsv_neutral<-get_est_all(res_arsv_neutral_list)
rm(res_hhs1_arsv_neutral,res_hhs2_arsv_neutral,res_hhs3_arsv_neutral,
   res_hhs4_arsv_neutral,res_hhs5_arsv_neutral,res_hhs6_arsv_neutral,
   res_hhs7_arsv_neutral,res_hhs8_arsv_neutral,res_hhs9_arsv_neutral,res_hhs10_arsv_neutral,res_arsv_neutral_list)

######
res_arsv_psi_list<-list(res_hhs1_arsv_psi,res_hhs2_arsv_psi,res_hhs3_arsv_psi,
                        res_hhs4_arsv_psi,res_hhs5_arsv_psi,
                        res_hhs7_arsv_psi,res_hhs8_arsv_psi,res_hhs9_arsv_psi,res_hhs10_arsv_psi)
est_res_arsv_psi<-get_est_all(res_arsv_psi_list)
rm(res_hhs1_arsv_psi,res_hhs2_arsv_psi,res_hhs3_arsv_psi,
   res_hhs4_arsv_psi,res_hhs5_arsv_psi,res_hhs6_arsv_psi,
   res_hhs7_arsv_psi,res_hhs8_arsv_psi,res_hhs9_arsv_psi,res_hhs10_arsv_psi,res_arsv_psi_list)
###,

res_arsv_chi_list<-list(res_hhs1_arsv_chi,res_hhs3_arsv_chi,res_hhs5_arsv_chi,
                        res_hhs6_arsv_chi,res_hhs7_arsv_chi,res_hhs9_arsv_chi,res_hhs10_arsv_chi)
est_res_arsv_chi<-get_est_all(res_arsv_chi_list)
rm(res_hhs1_arsv_chi,res_hhs2_arsv_chi,res_hhs3_arsv_chi,
   res_hhs4_arsv_chi,res_hhs5_arsv_chi,res_hhs6_arsv_chi,
   res_hhs7_arsv_chi,res_hhs8_arsv_chi,res_hhs9_arsv_chi,res_hhs10_arsv_chi,res_arsv_chi_list)

rbind(est_res_arsv_neutral,est_res_arsv_coinfect,est_res_arsv_psi,est_res_arsv_chi)%>%
  ggplot(aes(x=region,y=loglik))

res_hhs2_arsv_chi
res_hhs4_arsv_chi
res_hhs8_arsv_chi,
res_brsv_coinfect_list<-list(
  res_hhs1_brsv_coinfect,res_hhs2_brsv_coinfect,res_hhs3_brsv_coinfect,
  res_hhs4_brsv_coinfect,res_hhs5_brsv_coinfect,res_hhs6_brsv_coinfect,
  res_hhs7_brsv_coinfect,res_hhs8_brsv_coinfect,res_hhs9_brsv_coinfect,res_hhs10_brsv_coinfect)
est_res_brsv_coinfect<-get_est_all(res_brsv_coinfect_list)
rm(res_brsv_coinfect_list,res_hhs1_brsv_coinfect,res_hhs2_brsv_coinfect,res_hhs3_brsv_coinfect,
   res_hhs4_brsv_coinfect,res_hhs5_brsv_coinfect,res_hhs6_brsv_coinfect,
   res_hhs7_brsv_coinfect,res_hhs8_brsv_coinfect,res_hhs9_brsv_coinfect,res_hhs10_brsv_coinfect)
######
res_brsv_neutral_list<-list(res_hhs1_brsv_neutral,res_hhs2_brsv_neutral,res_hhs3_brsv_neutral,
                            res_hhs4_brsv_neutral,res_hhs5_brsv_neutral,res_hhs6_brsv_neutral,
                            res_hhs7_brsv_neutral,res_hhs8_brsv_neutral,res_hhs9_brsv_neutral,res_hhs10_brsv_neutral)
est_res_brsv_neutral<-get_est_all(res_brsv_neutral_list)
rm(res_hhs1_brsv_neutral,res_hhs2_brsv_neutral,res_hhs3_brsv_neutral,
   res_hhs4_brsv_neutral,res_hhs5_brsv_neutral,
   res_hhs7_brsv_neutral,res_hhs8_brsv_neutral,res_hhs9_brsv_neutral,res_hhs10_brsv_neutral,res_brsv_neutral_list)

######
res_brsv_psi_list<-list(res_hhs1_brsv_psi,res_hhs2_brsv_psi,res_hhs3_brsv_psi,
                        res_hhs4_brsv_psi,res_hhs5_brsv_psi,
                        res_hhs7_brsv_psi,res_hhs8_brsv_psi,res_hhs9_brsv_psi,res_hhs10_brsv_psi)
est_res_brsv_psi<-get_est_all(res_brsv_psi_list)
rm(res_hhs1_brsv_psi,res_hhs2_brsv_psi,res_hhs3_brsv_psi,
   res_hhs4_brsv_psi,res_hhs5_brsv_psi,res_hhs6_brsv_psi,
   res_hhs7_brsv_psi,res_hhs8_brsv_psi,res_hhs9_brsv_psi,res_hhs10_brsv_psi,res_brsv_psi_list)
###,

res_brsv_chi_list<-list(res_hhs1_brsv_chi,res_hhs1_brsv_chi,res_hhs3_brsv_chi,res_hhs4_brsv_chi,res_hhs7_brsv_chi,res_hhs8_brsv_chi,res_hhs9_brsv_chi,res_hhs10_brsv_chi)
est_res_brsv_chi<-get_est_all(res_brsv_chi_list)
rm(res_hhs1_brsv_chi,res_hhs2_brsv_chi,res_hhs3_brsv_chi,
   res_hhs4_brsv_chi,res_hhs5_brsv_chi,res_hhs6_brsv_chi,
   res_hhs7_brsv_chi,res_hhs8_brsv_chi,res_hhs9_brsv_chi,res_hhs10_brsv_chi,res_brsv_chi_list)

rbind(est_res_brsv_neutral,est_res_brsv_coinfect,est_res_brsv_psi,est_res_brsv_chi)%>%
  ggplot(aes(x=region,y=loglik))
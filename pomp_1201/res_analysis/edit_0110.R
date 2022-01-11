load("pomp_longphi_result/res_arsv_chi/res_hhs4_arsv_chi.rds")
res_hhs4_arsv_chi$DEobj$optim$bestval<-12603.05
save(res_hhs4_arsv_chi, file = "pomp_longphi_result/res_arsv_chi/res_hhs4_arsv_chi.rds")
load("pomp_longphi_result/res_arsv_coinfect/res_hhs4_arsv_coinfect.rds")
res_hhs4_arsv_coinfect$DEobj$optim$bestval<-12591.42
save(res_hhs4_arsv_coinfect, file = "pomp_longphi_result/res_arsv_coinfect/res_hhs4_arsv_coinfect.rds")
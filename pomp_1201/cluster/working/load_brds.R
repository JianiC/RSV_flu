res_brsv_coinfect<-list.files(path="pomp_result_1213/res_brsv_coinfect/",pattern=".rds")
for (i in 1:length(res_brsv_coinfect)){
  load(paste("pomp_result_1213/res_brsv_coinfect/",res_brsv_coinfect[i],sep=""))
  print(paste0("load ",res_brsv_coinfect[i]))
}


res_brsv_neutral<-list.files(path="pomp_result_1213/res_brsv_neutral/",pattern=".rds")
for (i in 1:length(res_brsv_neutral)){
  load(paste("pomp_result_1213/res_brsv_neutral/",res_brsv_neutral[i],sep=""))
  print(paste0("load ",res_brsv_neutral[i]))
}


res_brsv_psi<-list.files(path="pomp_result_1213/res_brsv_psi/",pattern=".rds")
for (i in 1:length(res_brsv_psi)){
  load(paste("pomp_result_1213/res_brsv_psi/",res_brsv_psi[i],sep=""))
  print(paste0("load ",res_brsv_psi[i]))
}


res_brsv_chi<-list.files(path="pomp_result_1213/res_brsv_chi/",pattern=".rds")
for (i in 1:length(res_brsv_chi)){
  load(paste("pomp_result_1213/res_brsv_chi/",res_brsv_chi[i],sep=""))
  print(paste0("load ",res_brsv_chi[i]))
}

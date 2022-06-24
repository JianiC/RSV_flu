## genetic interaction
## load skyride data
library(ggplot2)
library(ggpubr)
library(ggtree)
library(treeio)
library(scales)
##

## load skyride data
## RSV skyride
ON_skyride<-read.table('Genetic_analysis/RSV_flu_Beast/RSV/gisaid_ON1_north_100901_190430_skyride.txt', skip = 1, header =TRUE, sep ='\t')
BA_skyride<-read.table('Genetic_analysis/RSV_flu_Beast/RSV/gisaid_BA_north_100901_190430_skyride.txt', skip = 1, header =TRUE, sep ='\t')

H3_skyride<-read.table('Genetic_analysis/RSV_flu_Beast/Flu/gisaid_H3_100901_190401_skyride.txt', skip = 1, header =TRUE, sep ='\t')
H1_skyride<-read.table('Genetic_analysis/RSV_flu_Beast/Flu/gisaid_H1_100901_190430_skyride.txt', skip = 1, header =TRUE, sep ='\t')

Vic_skyride<-read.table('Genetic_analysis/RSV_flu_Beast/Flu/gisaid_Vic_100901_190430_skyride.txt', skip = 1, header =TRUE, sep ='\t')
Yam_skyride<-read.table('Genetic_analysis/RSV_flu_Beast/Flu/gisaid_Yam_100901_190430_skyride.txt', skip = 1, header =TRUE, sep ='\t')

## approaximate skyride
skyride_approx<-function(skyride,timepoint){
  approx_out<-approx(skyride$Time,skyride$Median,xout=timepoint)
  return(approx_out)
}



## specify x out
skyride_timpoint<-seq(2010, 2019, by = 0.1)
ON_skyride_approx<-skyride_approx(ON_skyride,skyride_timpoint)
BA_skyride_approx<-skyride_approx(BA_skyride,skyride_timpoint)

H3_skyride_approx<-skyride_approx(H3_skyride,skyride_timpoint)
H1_skyride_approx<-skyride_approx(H1_skyride,skyride_timpoint)

Vic_skyride_approx<-skyride_approx(Vic_skyride,skyride_timpoint)
Yam_skyride_approx<-skyride_approx(Yam_skyride,skyride_timpoint)

## test
plot(H1_skyride_approx$x,H1_skyride_approx$y)
#################################################################
## test interaction: skyride vs skyride
#################################################################


## start with pearson correlation
## check normial distribution
set.seed(0)
shapiro.test(ON_skyride_approx$y)
shapiro.test(BA_skyride_approx$y)
shapiro.test(H3_skyride_approx$y)
shapiro.test(H1_skyride_approx$y)
shapiro.test(Vic_skyride_approx$y)
shapiro.test(Yam_skyride_approx$y)

## pairwise skyrdie correlation test

skyride_corr<-function(skyride_approx1,skyride_approx2){
  res<-cor.test(skyride_approx1,skyride_approx2,method="pearson")
  out<-data.frame(method="pearson",
                  corre=res$estimate,
                  CI_low=res$conf.int[1],
                  CI_high=res$conf.int[2])
  return(out)
}

## store the skyride approx to a list 
idx_skyride<-list(ON=ON_skyride_approx$y,BA=BA_skyride_approx$y,H3=H3_skyride_approx$y,H1=H1_skyride_approx$y,
                  Vic=Vic_skyride_approx$y,Yam=Yam_skyride_approx$y)

## loop through all index
pathogen_idx=combn(seq(1,length(idx_skyride)), 2, FUN = NULL, simplify = TRUE)

## loop through all list to caluculate the coorelation
df_skyride_cor<-data.frame()
for(i in seq(1, ncol(pathogen_idx))) {

  id1=pathogen_idx[1,i]
  id2=pathogen_idx[2,i]
  corr_out <- skyride_corr(idx_skyride[[id1]],idx_skyride[[id2]])
  corr_out$virus1 =names(idx_skyride[id1])
  corr_out$virus2 =names(idx_skyride[id2])

  ## concatenate loop result together
  df_skyride_cor<-rbind(df_skyride_cor,corr_out)
}

df_skyride_cor
write.csv(df_skyride_cor,"Genetic_analysis/skyride_cor.csv", row.names = FALSE)
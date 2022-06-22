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

skyride_corr<-function(skyride_approx1,skyride_approx2,label1,label2){
  cor<-cor.test(skyride_approx1$y,skyride_approx2$y,method="pearon")
  out<-data.frame(virus1=label1,
                  virus2=label2,
                  method="pearson",
                  corre=res$estimate,
                  CI_low=res$conf.int[1],
                  CI_high=res$conf.int[2])
  return(out)
  
  
}
skyride_cor1<-

# Script to visualize csv output file(export from Tracer )from BEAST Bayesian Skyride coalescent analysis
# Input: csv file (remove first line with	GMRF Skyride: xxxxxxxx.log)
# Output: ggplot of skyride wiht 95% HPD 


library(ggplot2)
library(lubridate)


dat <- read.table('Genetic_data/usa_12-14/Raxml/BEAST/rsva_12_14_v2/skyride_plot.txt', skip = 1, header =TRUE, sep ='\t')

dat$date <- as.Date(paste(date_decimal(dat$Time)),"%Y-%m-%d")

p1 <- ggplot(dat, aes(y=Mean, x=date)) + geom_line() + 
  geom_ribbon(aes(ymin=Lower, ymax=Upper, alpha = 0.2)) +
  theme(legend.position = "bottom") + labs(alpha="95% hpd") + 
  labs(title = "Skyride coalescent") + 
  scale_y_continuous(name ="Log mean pop (N)", trans = 'log10') 


p1
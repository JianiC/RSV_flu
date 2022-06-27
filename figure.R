## Generate Figure in PDF format
## Color palettes, themes, margins, labels, etc.

## PNAS specs
## For details see https://blog.pnas.org/digitalart.pdf
## max height = 9in / 22.5cm
## knitr expects inches by default
## scale up by *1.56 here, down by 0.64 in final tex
#figwidth <- c(8.7, 11.4, 17.8)/2.54
figwidth <- c(13.6, 17.8, 27.7)/2.54
## Intended target dims
outwidth <- c('8.7cm', '11.4cm', '17.8cm')

## cowplot::plot_grid(label_size=)
lab.size = 113.6

## ggplot: geom_text size
lab.size.small=4
## for grid.text, match cowplot::plot_grid labels
gp.label <- gpar(fonsize=lab.size, fontface='bold')

## use consistent theme throughout
gg.font.base <-12
gg.theme <- (
  theme_bw(base_size=gg.font.base)
  + theme(
    legend.background=element_blank(),
    panel.grid.minor.x=element_blank(),
    panel.grid.minor.y=element_blank()
  )
  # + ...
)

##################################################################
## figure 1 data figure

source("datafigure.R", chdir = TRUE) 


pdf(file = "Figure1_datafigure_cor.pdf",
    width = 8.7/2.54,
    height = 15/2.54)
#plot(data_figure)
plot(data_cor_figure)

dev.off()

##################################################################
## figure2 Skyride and phylogeny

source("plot_tree_pact.R", chdir = TRUE) 

pdf(file = "Figure2_phylogeny.pdf",
    width = 17.8/2.54,
    height = 13/2.54)
plot(p_skyrdie_pact)

dev.off()

##################################################################
## figure3 Pomp model simulaiton plot

source("dynamic_simulate.R", chdir = TRUE) 
source("point_simulate.R", chdir = TRUE) 

pdf(file = "Figure3_simulation.pdf",
    width = 9/2.54,
    height = 12/2.54)
plot(simu_plot)

dev.off()




##################################################################
## figure4 Fit with different hypothesis ~ example with HHS region 1

source("hypo_fit_plot.R")

pdf(file = "Figure4_hypofit.pdf",
    width = 9/2.54,
    height = 12/2.54)
plot(p_fit_hypo_repre)

dev.off()

##################################################################
## figure5 loglike profile plot  ~ example with HHS region 1
source("loglik_plot_0222.R")
pdf(file = "Figure5_loglik.pdf",
    width = 8.7/2.54,
    height = 10/2.54)
plot(loglik_CI_repre)

dev.off()


###############################################################
## supplemmentary Fig2
source("plot_tree_pact.R", chdir = TRUE) 

pdf(file = "S2Fig_pactmeam.pdf",
    width = 9/2.54,
    height = 13/2.54)
plot(pact_mean)

dev.off()




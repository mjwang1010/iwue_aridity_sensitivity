# Settings for all PROs

rm(list = ls())
# Packages------------

library(tidyverse)
library(patchwork)
library(trend)
library(boot)

inpath1 <- 'Path for input files'
setwd(inpath1)

# Functions ---------------------------------------------------------------

# Bootstrapping
boot_mean_fun <- function(data, index) {
  return(mean(data[index], na.rm = T))
}
boot_ci_fun <- function(data, repnum = 1000, type = 'up') {
  # CI for mean values
  bootobj <- boot(data, boot_mean_fun, R = repnum)
  bootci <- boot.ci(bootobj, type = 'norm')
  bound <- bootci$normal
  # 95 CI upper or low
  if (type == 'up') {
    return(bound[, 3])
  } else if (type == 'low') {
    return(bound[, 2])
  }
}

# Plotting ---------------

a.color <- '#E64B35FF'
g.color <- '#4DBBD5FF'
a.label <- 'Angiosperm'
g.label <- 'Gymnosperm'

all.color <- '#969696'

theme_char <- theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.tag = element_text(size = 7.5),
    plot.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7.5),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6.5)
  )
labs_geo_char <- labs(x = expression('Longitude' ~ (degree * 'E')),
                      y = expression('Latitude' ~ (degree * 'N')))

labs_ts_char <- labs(x = 'Year', 
                     y = expression(S[italic(W)]), 
                     fill = NULL, 
                     color = NULL)

limits_ts_char <- lims(x = c(1951, 2010))

line_ts_char <- geom_hline(
  yintercept = 0, 
  linetype = 'dashed'
)

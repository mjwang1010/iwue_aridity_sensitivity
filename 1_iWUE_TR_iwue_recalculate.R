# Reclaculate the iWUE according to the recommended CO2 data

rm(list = ls())

##load libraries
library(doBy)
library(data.table)
library(reshape2)
library(plyr)
library(tidyverse)

options("scipen"=100, "digits"=4)

# Load in collated input data ---------------------------------------------

setwd('Path for the input file folder')

# Original stable carbon isotope data
# NOTE: the data is already contained in the final file
# 'TR_iwue_newco2_processed_iwue_full_iav_20nyrs_all.csv'
x <- read_csv('The input file name') 


##partition based on data type (will recombine later after calculation of common Delta)
xu <- subset(x, datatype=="u") #uncorreced d13Cplant
xc <- subset(x, datatype=="c") #corrected d13C
xt <- subset(x, datatype=="t") #'triangle' Delta (capital Delta)
xw <- subset(x, datatype=="w") #iwUE
xi <- subset(x, datatype=="i") #ci


##### y=uncorrected d13Cplant -----
#####     (1) compute Delta

delta_uncor_fun <- function(d13ca, d13c) {
  
  delta <- (d13ca - d13c) / (1 + 0.001*d13c)
  
  return(delta)
}
xu$Delta <- delta_uncor_fun(xu$d13Catm, xu$y)
xu$Delta_new <- delta_uncor_fun(xu$d13Catm_recom_interpolated, xu$y)


##### y=corrected d13Cplant -----
#####     (1) de-correct using d13Catm
#####     (2) compute Delta from decorrected value

delta_cor_fun <- function(d13ca, d13ccor) {
  
  dconst <- -6.4 # (PIN d13ca)
  d13c <- d13ccor + d13ca - dconst
  delta <- (d13ca - d13c) / (1 + 0.001*d13c)
  
  return(delta)
}

xc$Delta <- delta_cor_fun(xc$d13Catm, xc$y) 
xc$Delta_new <- delta_cor_fun(xc$d13Catm_recom_interpolated, xc$y)


##### y=WUEi-----
#####     (1) back out Delta as a + (b-a)*(1 - 1.6*WUEi/ca)

delta_iwue_fun <- function(ca, iwue) {
  
  a <- 4.4 
  b <- 27
  
  # Use the simple model to re-estimate the delta
  delta <- a + (b - a)*(1 - 1.6*iwue/ca)
  
  return(delta)
}

delta_iwue_fun2 <- function(ca, iwue) {
  
  # For 'study-34'
  a <- 4.4 
  b <- 29
  
  # Use the simple model to re-estimate the delta
  delta <- a + (b - a)*(1 - 1.6*iwue/ca)
  
  return(delta)
}

#study #34 used b = 29
xw_sub1 <- xw %>% 
  filter(study != 34)
xw_sub2 <- xw %>% 
  filter(study == 34)

xw_sub1$Delta = delta_iwue_fun(xw_sub1$ca, xw_sub1$y)
xw_sub1$Delta_new <- delta_iwue_fun(xw_sub1$ca_recom, xw_sub1$y)

xw_sub2$Delta = delta_iwue_fun2(xw_sub2$ca, xw_sub2$y)
xw_sub2$Delta_new <- delta_iwue_fun2(xw_sub2$ca_recom, xw_sub2$y)

# Re-combine the two
xw <- rbind(xw_sub1, xw_sub2)


##### y=ci-----
#####     (1) back out Delta as a + (b-a)*ci/ca

delta_ci_fun <- function(ca, ci) {
  
  a <- 4.4
  b <- 27
  
  # Use the simple model to re-estimate the delta
  delta <- a + (b - a)*ci/ca
  
  return(delta)
}

xi$Delta <- delta_ci_fun(xi$ca, xi$y)
xi$Delta_new <- delta_ci_fun(xi$ca_recom, xi$y)


##### y=Delta-----
#####     (1) keep as is
xt$Delta <- xt$y
xt$Delta_new = xt$y


######################################################
######################################################
#####
##### //done deprocessing data (all now Delta on a common basis)
##### //recombine and save

##recombine dfs 
x_new <- Reduce(function(x,y) merge(x,y,all=TRUE), list(xu, xc, xt, xw, xi))


# Calculations for iWUE ---------------------------------------------------


# (1) estimate A/ca needed for 'full' isotope model -----------------------
# Assume doubling ratio (DR) = 1.45 from Keeling et al 2017

A_ca_fun <- function(ca) {
  
  doublingratio=1.45
  beta=log(doublingratio)/log(2)-1
  A_ca280 = 9/280 #assume pre-industrial A=9 
  
  A_ca <- A_ca280*((ca/280)^beta)
  
  return(A_ca)
}

x_new$A_ca <- A_ca_fun(x_new$ca)
x_new$A_ca_new <- A_ca_fun(x_new$ca_recom)


# (2) compute ci/ca from Delta, and iWUE from ci/ca ---------------------------


ci_ca_fun <- function(delta, A_ca, ca) {
  
  # Set parameters for full isotope model
  
  gm = 0.2 #mesophyll conductance (mol m-2 s-1)
  f = 12 #photorespiration discrimination term (per mille)
  gammastar = 43 #photorespiratory CO2 compensation point (ppm)
  a = 4.4 #diffusive discrimination (per mille)
  b = 30 #Rubisco discrimination (per mille)
  am = 1.8 # (per mille)
  
  ci_ca = (delta - a + (b - am)*A_ca/gm + f*gammastar/ca)/(b - a)
  
  return(ci_ca)
}

x_new$ci_ca <- ci_ca_fun(x_new$Delta, x_new$A_ca, x_new$ca)
x_new$ci_ca_new <- ci_ca_fun(x_new$Delta_new, x_new$A_ca_new, x_new$ca_recom)

x_new$iWUE <- x_new$ca * (1 - x_new$ci_ca) / 1.6
x_new$iWUE_new <- x_new$ca_recom * (1 - x_new$ci_ca_new) / 1.6

# Using simple isotope model
a <- 4.4
b_simple <- 27

x_new$ci_ca_simple <- (x_new$Delta - a)/(b_simple - a)
x_new$ci_ca_simple_new <- (x_new$Delta_new - a)/(b_simple - a)

x_new$iWUE_simple <- x_new$ca * (1 - x_new$ci_ca_simple) / 1.6
x_new$iWUE_simple_new <- x_new$ca_recom * (1 - x_new$ci_ca_simple_new) / 1.6

#Save results -----------------------------------------------------------


write_csv(x_new, 
          path = 'The output file name')

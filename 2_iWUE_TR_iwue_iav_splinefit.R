# Spline method to estimate the interannual variations of the iWUE
# Or exponential function to fit the trend
# Then, de-trend

rm(list = ls())
library(tidyverse)
library(dplR)

# File paths and names ----------------------------------------------------

infilepath1 <- 'Path of the input file folder'
infilename1 <- 'Path of the input file of the iWUE dataset'

outfilepath1 <- infilepath1 # Output file folder

# Load data ---------------------------------------------------------------

# iWUE data
# NOTE: the data is already contained in the final file
# 'TR_iwue_newco2_processed_iwue_full_iav_20nyrs_all.csv'
df1 <- read_csv(infilename1)


# Estimate the interannual trend ------------------------------------------

trend_calc_fun <- function(x, xyear, flag) {
  
  # Input - 
  # x - target variables (y)
  # xyear - the years (x)
  
  df <- data.frame(y = x, 
                   t = xyear)
  df1 <- df[complete.cases(x), ]
  # Reorder
  inds_order <- order(df1$t)
  df1 <- df1[inds_order, ]
  
  if (flag == 'E') {
    # Exponential fitting
    # y = a0 + a1*exp((x - 1850) / a2)
    expmod <- nls(y ~ a0 + a1 * exp((t - 1850) * a2),
                  data = df1,
                  start = list(a0 = min(df1$y), a1 = 0.1*mean(df1$y), a2 = 0.01))
    # Fitting values
    dffit <- data.frame(t = df$t,
                        yfit = predict(expmod, df$t))
    
  }
  
  return(dffit$yfit)
  
}

# For each tree 
tree_name <- unique(df1$tree)
num_tree <- length(unique(df1$tree))

df_datatype <- data.frame() # Save trees with more than one datatype
df_trendiav1 <- data.frame() # trend and interannual variability

spaniyear <- 1901 # Ranges of the years
spanfyear <- max(df1$year) # 2015

for (i in 1:num_tree) {
  
  # Find whether there are multiple data types
  dftemp1 <- df1 %>% 
    filter(tree == tree_name[i])
  
  datatypetemp1 <- unique(dftemp1$datatype)
  num_datatypetemp1 <- length(datatypetemp1)
  if (num_datatypetemp1 > 1) {
    
    dftemp1 <- dftemp1 %>% 
      filter(datatype == 'u')
  }

  # Reorder
  inds_order <- order(dftemp1$year)
  dftemp1 <- dftemp1[inds_order, ]
  
# Trend detection ---------------------------------------------------------

  
  # High-pass filtering ----
  
  # Parameter settings for cubic smoothing spline fitting
  nyrs <- 20 # Original setting in the reference
  
  # Combine data in duplicated years
  dftemp1_uniyear <- dftemp1 %>% 
    group_by(year) %>% 
    summarize(iWUE_new_uni = mean(iWUE_new, na.rm = T))
    # summarise(iWUE_simple_new_uni = mean(iWUE_simple_new, na.rm = T)) # iWUE simple model
  
  data_detrend1 <- dplR::detrend(
    data.frame(y = dftemp1_uniyear$iWUE_new_uni, row.names = dftemp1_uniyear$year),
    # data.frame(y = dftemp1_uniyear$iWUE_simple_new_uni, row.names = dftemp1_uniyear$year),
    method = 'Spline',
    nyrs = nyrs, # The setting in the reference (e.g., Saurer et al., 2014)
    f = 0.5,
    difference = T, 
    return.info = T)
  
  df_detrend1 <- data.frame(year = dftemp1_uniyear$year, 
                            iWUE_new_trend = data_detrend1$curves$y, 
                            iWUE_new_iav = data_detrend1$series$y)
  
  # Set the ranges for years and count the data
  # from 1901 to 2015
  df_detrend1_sub1 <- df_detrend1 %>% 
    filter(year >= spaniyear & year <= spanfyear) %>% 
    mutate(nobs = length(year))
  
  # Connect with the original data frame
  dftemp2 <- dftemp1 %>% 
    inner_join(df_detrend1_sub1, by = 'year')
  
  df_trendiav1 <- rbind(df_trendiav1, dftemp2)
}


# Save the data -----------------------------------------------------------

write_csv(df_trendiav1, 
          file = paste0(outfilepath1, 'tree_ring_iwue_newco2_processed_iwue_full_iav_', nyrs, 
                        'nyrs_all.csv'))

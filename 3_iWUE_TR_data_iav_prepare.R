# Prepare data for sensitivity analysis
# Connect IAV of iWUE with SPEI 
# and estimate IAV of SPEI

source('Path for iWUE_setting.R')
infilepath1 <- 'Path of the input file folder'
outfilepath1 <- infilepath1

# Load data ---------------------------------------------------------------

# Get iWUE IAV
df_iwue <- read_csv(paste0('TR_iwue_newco2_processed_iwue_full_iav_20nyrs_all.csv'))

# Get growing-season SPEI data 
df2 <- read_csv(paste0('TR_spei6mon_gs3.csv'))

# Estimate IAV ------------------------------------------------------------

iav_estimate_fun <- function(year, vecx) {
  # Estimate the IAV for one vector
  # Return the IAV vector
  
  df1 <- data.frame(
    year = year, 
    x = vecx
  )
  df_caltrend <- df1 %>% 
    drop_na()
  
  if (dim(df_caltrend)[1] == 0) {
    print(paste0('No available data during the period'))
    x_iav <- as.numeric(rep(NA, length(df1$x)))
    return(x_iav)
  } 
  
  slope <- sens.slope(as.ts(df_caltrend$x)) # Get slope
  interc <- median(df_caltrend$x - slope$estimates * df_caltrend$year) # Get interception
  
  # Estimate the IAV of the AI
  x_iav <- df1$x - (slope$estimates * df1$year + interc)
  return(x_iav)
  
}

# Connect different datasets ----------------------------------------------

df_iwue_var <- df_iwue %>% 
  group_by(tree, year) %>% 
  left_join(df2 %>% 
              dplyr::select(-all_of(c('lat', 'lon'))))
# Rename the variable
df_iwue_var <- df_iwue_var %>% 
  rename(spei6mon_mean = SPEI_mean)
# Get the IAV
df_iwue_var_iav <- df_iwue_var %>% 
  group_by(tree) %>%
  mutate(
    spei6mon_iav = iav_estimate_fun(year, spei6mon_mean)
  )

# Save results ------------------------------------------------------------

write_csv(
  df_iwue_var_iav, 
  file = paste0(outfilepath1, 'TR_iwueiav_v3.csv')
)

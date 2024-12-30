# Estimate differences in sensitivity between the first and last periods
# For SW_SPEI
# Input: sensitivity in each moving window

source('Path for iWUE_setting.R')
infilepath1 <- 'Path of the input file folder'
outfilepath1 <- infilepath1

# Load data ---------------------------------------------------------------

# Tree clades
df_iav <- read_csv(paste0(infilepath1, 'TR_iwueiav_v3.csv'))
# Sensitivities in moving windows
winsize <- 22
var_iav <- 'spei6moniav_gs3_'
df_winsens <- read_csv(paste0(infilepath1, 
                              'TR_iwueiav_slr_', 
                              var_iav, 
                              'winenvi_win', 
                              winsize, 
                              'sen.csv'))

# Climatic conditions at sites
df_cond <- read_csv(paste0(infilepath1, 
                           'TR_climate_condition_geoinfo_cru_1951_2015.csv'))

# Add clades to the data
df_clade <- df_iav %>% 
  dplyr::select(c(tree, group)) %>% 
  distinct()
df_winsens_clade <- df_winsens %>% 
  left_join(df_clade, 
            by = c('tree' = 'tree'))

worldmap <- map_data("world")

# Difference --------------------------------------------------------------

# from year 1951 to year 2010
start_period_year <- 1951
end_period_year <- 2010
# 'win_year' center of the moving window
start_win_year <- start_period_year+winsize/2
end_win_year <- end_period_year-(winsize/2-1)

df_win_twoperiod <- df_winsens_clade %>% 
  filter(win_year %in% c(start_win_year, end_win_year))

# Same samples in both periods
df_numwindiff <- df_win_twoperiod %>% 
  drop_na(sens) %>% 
  group_by(tree) %>% 
  summarize(num_win_year = n()) %>% 
  filter(num_win_year == 2)

# Estimate the difference

var_windiff_fun <- function(df, tree_name, start_win_year, var_name) {
  # Estimate the difference between the first and last periods
  # df: long tables, data from two periods
  # tree_name: tree names which have data in the two periods
  # start_win_year: indicates the first period
  # var_name: choose which variables to estimate the difference
  
  df_windiff <- df %>% 
    filter(tree %in% tree_name) %>% 
    # Rename years
    mutate(
      win_year = if_else(win_year == start_win_year, 'Early_period', 'Late_period')
    ) %>% 
    dplyr::select(all_of(c('tree', 'win_year', var_name, 'group'))) %>% 
    pivot_wider(
      names_from = win_year, 
      values_from = var_name
    )
  # Estimate the difference
  df_windiff <- df_windiff %>% 
    mutate(
      diff = Late_period - Early_period
    ) 
  return(df_windiff)
}

# Sensitivity
df_sens_windiff <-
  var_windiff_fun(df_win_twoperiod,
                  df_numwindiff$tree,
                  start_win_year,
                  'sens')

# Correlation coefficients
df_rcor_windiff <-
  var_windiff_fun(df_win_twoperiod,
                  df_numwindiff$tree,
                  start_win_year,
                  'Rcor')

# Difference in other variables between the two periods
# Precipitation
df_pre_windiff <- var_windiff_fun(df_win_twoperiod, df_numwindiff$tree, 
                                  start_win_year, 'PRE_ave')
df_pre_windiff <- df_pre_windiff %>% 
  dplyr::select(tree, diff) %>% 
  rename(prediff = diff)

# Aridity index
df_ai_windiff <- var_windiff_fun(df_win_twoperiod, df_numwindiff$tree, 
                                 start_win_year, 'AI_ave')
df_ai_windiff <- df_ai_windiff %>% 
  dplyr::select(tree, diff) %>% 
  rename(aidiff = diff)

# Temperature
df_tmp_windiff <- var_windiff_fun(df_win_twoperiod, df_numwindiff$tree, 
                                  start_win_year, 'TMP_ave')
df_tmp_windiff <- df_tmp_windiff %>% 
  dplyr::select(tree, diff) %>% 
  rename(tmpdiff = diff)

# Potential factors -------------------------------------------------------

# Climatic conditions
df_windiff_conddiff <- df_sens_windiff %>% 
  left_join(df_cond) %>% 
  left_join(df_pre_windiff) %>% 
  left_join(df_ai_windiff) %>% 
  left_join(df_tmp_windiff)

df_rcor_windiff_conddiff <- df_rcor_windiff %>% 
  left_join(df_cond) %>% 
  left_join(df_pre_windiff) %>% 
  left_join(df_ai_windiff) %>% 
  left_join(df_tmp_windiff)


# Save data ---------------------------------------------------------------

# Difference in sensitivity to 6month-scale SPEI
write_csv(
  df_windiff_conddiff, 
  file = paste0(outfilepath1, 
                'TR_iwue_sensdiff_', 
                var_iav, 
                'win', winsize, 'sens_', 
                start_period_year, '_', end_period_year, '.csv')
)



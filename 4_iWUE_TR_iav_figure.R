# IAV of anomalies
# Figure 1 in MS 2023-04

source('Path for iWUE_setting.R')

infilepath1 <- 'Path of the input file folder'
outfilepath1 <- infilepath1

# Load data ---------------------------------------------------------------

# IAV in iWUE and AI or other variables
df_iav <- read_csv(paste0(infilepath1, 'TR_iwueiav_v3.csv'))

# Prepare data ------------------------------------------------------------

# Interannual variation in anomalies
var_name <- 'spei6mon_iav'
df_iav_mean <- df_iav %>% 
  drop_na(iWUE_new_iav, as.name(var_name)) %>% 
  group_by(year) %>% 
  summarise(
    iWUE_iav_mean = mean(iWUE_new_iav, na.rm = T), 
    var_iav_mean = mean(get(as.name(var_name)), na.rm = T), 
    iWUE_iav_mean_up = boot_ci_fun(iWUE_new_iav, type = 'up'), 
    iWUE_iav_mean_low = boot_ci_fun(iWUE_new_iav, type = 'low'), 
    var_iav_mean_up = boot_ci_fun(get(as.name(var_name)), type = 'up'), 
    var_iav_mean_low = boot_ci_fun(get(as.name(var_name)), type = 'low'), 
    iav_n = n()
  )


# from year 1951 to year 2010
start_period_year <- 1951
end_period_year <- 2010

df_iav_mean_sub <- df_iav_mean %>% 
  filter(year>=start_period_year & 
           year<=end_period_year)

# Equal sample sizes -----------

df_iavnum <-  df_iav %>%
  filter(year >= start_period_year &
           year <= end_period_year) %>%
  drop_na(iWUE_new_iav, as.name(var_name)) %>%
  group_by(tree) %>%
  summarise(num_iav = n()) %>%
  # Set "==" if there is an equal sample size
  filter(num_iav > 40 )

# 
df_iav_mean_limitsample <- df_iav %>% 
  filter(year >= start_period_year &
           year <= end_period_year) %>%
  filter(tree %in% df_iavnum$tree) %>% 
  drop_na(iWUE_new_iav, as.name(var_name)) %>% 
  group_by(year) %>% 
  summarise(
    iWUE_iav_mean = mean(iWUE_new_iav, na.rm = T), 
    var_iav_mean = mean(get(as.name(var_name)), na.rm = T), 
    iWUE_iav_mean_up = boot_ci_fun(iWUE_new_iav, type = 'up'), 
    iWUE_iav_mean_low = boot_ci_fun(iWUE_new_iav, type = 'low'), 
    var_iav_mean_up = boot_ci_fun(get(as.name(var_name)), type = 'up'), 
    var_iav_mean_low = boot_ci_fun(get(as.name(var_name)), type = 'low'), 
    iav_n = n()
  )

# Calculate R and sensitivity in sliding windows ---------------------------------------------

winsize <- 22
win_startyear <- start_period_year + floor(winsize/2)
win_endyear <- end_period_year - winsize + 1 + floor(winsize/2)
num_win <- win_endyear-win_startyear+1

df_win_sensrcor <- data.frame()
# For each window
for (j in 1:num_win) {
  
  # Set years
  dfwin <- df_iav_mean_sub %>% 
    filter(year >= (win_startyear - floor(winsize/2) - 1 + j) &
             year <= (win_startyear - floor(winsize/2) - 1 + j + winsize-1))
  
  df_regress1 <- dfwin %>% 
    dplyr::select(iWUE_iav_mean, var_iav_mean) %>%
    drop_na()
  
  # There should be enough data in a window
  # How to set the threshold?
  # Larger than half of the window size
  if (dim(df_regress1)[1] < winsize/2) {
    
    df_stat <- data.frame(
      year = NA, 
      sens = NA, 
      pval_lin = NA, 
      sens_cil = NA, 
      sens_ciu = NA, 
      Rcor = NA, 
      pval_rcor = NA
    )
    next
    
  } else {
    
    # Sensitivity
    regress1 <- lm(data = df_regress1, 
                   iWUE_iav_mean ~ var_iav_mean)
    regress1summ <- summary(regress1)
    regress1ci <- confint(regress1) # CI intervals
    # Correlation coefficient
    correlation1 <- cor.test(data = df_regress1, 
                             ~iWUE_iav_mean + var_iav_mean)
    
    df_stat <- data.frame(
      year = win_startyear + j-1, 
      sens = regress1summ$coefficients[, 1][2], 
      pval_lin = regress1summ$coefficients[, 4][2], 
      sens_cil = regress1ci[, 1][2], 
      sens_ciu = regress1ci[, 2][2], 
      Rcor = correlation1$estimate, 
      pval_rcor = correlation1$p.value, 
      Rcor_cil = correlation1$conf.int[1], 
      Rcor_ciu = correlation1$conf.int[2]
    )
    
  }
  
  df_win_sensrcor<- bind_rows(df_win_sensrcor, df_stat)
  
}


# Plotting ----------------------------------------------------------------

# a, Two y axises
if (var_name %in% 'sm0_10cm_iav') {
  secy_name <- '-Soil moisture anomaly'
  secfactor <- 1.5
} else if (var_name == 'spei6mon_iav') {
  secy_name <- "SPEI anomaly"
  secfactor <- 2.1
} else if (var_name == 'sm0_100cm_iav') {
  secy_name <- '-Soil moisture anomaly'
  secfactor <- 0.2
} else {
  secy_name <- ''
  secfactor <- 1
}

f_all_v2 <- ggplot(data = df_iav_mean_limitsample) + 
  geom_line(
    aes(x = year, 
        y = iWUE_iav_mean)
  ) + 
  geom_line(
    aes(x = year, 
        y = -var_iav_mean*secfactor), # Reverse the axis
    color = 'red'
  ) + 
  scale_y_continuous(
    # Features of the first axis
    name = expression(italic(W)~'anomaly ('*mu*'mol'~mol^-1*')'), 
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans = ~-./secfactor, name = secy_name)
  ) + 
  labs(x = 'Year', 
       tag = 'a') + 
  geom_hline(yintercept = 0, 
             linetype = 'dashed') + 
  theme_char

# b, Rcor time series
f_rcor <- ggplot(data = df_win_sensrcor) + 
  geom_line(
    aes(
      x = year, 
      y = Rcor
    )
  ) + 
  # Standard errors
  geom_ribbon(
    aes(x = year, 
        y = Rcor, 
        ymin = Rcor_cil, 
        ymax = Rcor_ciu
    ), 
    fill = 'gray', 
    alpha = 0.5
  ) + 
  labs_ts_char + 
  labs(y = expression(italic(R)), 
       x = "Year", 
       tag = 'b') + 
  limits_ts_char + 
  theme_char

# c, Sensitivity time series
f_sens <- ggplot(data = df_win_sensrcor) + 
  geom_line(
    aes(
      x = year, 
      y = sens
    )
  ) + 
  # confidence intervals
  geom_ribbon(
    aes(x = year, 
        y = sens, 
        ymin = sens_cil, 
        ymax = sens_ciu
    ), 
    fill = 'gray', 
    alpha = 0.5
  ) + 
  labs_ts_char + 
  labs(
    # y = expression(S[italic(W)]),
    y = expression(atop(
      paste(S[italic(W)]),
      paste('(', mu, 'mol ', mol ^ -1, ' per unit SPEI)')
    )), 
    x = "Year",
    tag = 'c'
  ) + 
  limits_ts_char + 
  theme_char

# Combine the figures
figure1_v1 <- f_all_v2 / (f_rcor + f_sens)

# Statistical test --------------------------------------------------------

summary(lm(data = df_win_sensrcor, 
           Rcor~year))
sens.slope(df_win_sensrcor$Rcor)

summary(lm(data = df_win_sensrcor, 
           sens~year))
sens.slope(df_win_sensrcor$sens)


# Save figures ------------------------------------------------------------


ggsave(
  figure1_v1,
  filename = paste0(
    outfilepath1,
    'Figures/',
    'fig1_iwueiav_', 
    var_name, 
    '_.jpeg'
  ),
  units = 'mm',
  dpi = 900,
  width = 165,
  height = 135
)


# Temporal variations in sensitivities to SPEI for all available interval periods
# By using a sliding window (size)
# Get environmental conditions (e.g., mean and SD) in each window

source('Path for iWUE_setting.R')
filepath1 <- 'Path of the input file folder'
metsource <- ''
winsize <- 20 # The moving window size
windataT <- winsize/2 # Threshold for the data number

# Load data ---------------------------------------------------------------

# IAV data including iWUE and SPEI
df1 <- read_csv(paste0(filepath1, 'TR_iwueiav_v3.csv'))

# Sliding window sensitivities --------------------------------------------

# For each tree calculate the sensitivities
tree_name <- unique(df1$tree)
num_tree <- length(tree_name)
spaniyear <- min(df1$year, na.rm = T)  # Initial year of the period
spanfyear <- max(df1$year, na.rm = T)

# The cencter of the first window
startyear <- spaniyear + floor(winsize/2)
endyear <- spanfyear - winsize + 1 + floor(winsize/2)
# How many windows for each tree
num_win <- endyear - startyear + 1

df_winsens <- data.frame()
# 
for (i in 1:num_tree) {
  
  dftemp1 <- df1 %>% 
    filter(tree == tree_name[i])
  
  # For each window calculate the sensitivities
  # estimate environmental conditions
  dftemp1_winsens <- data.frame(
    win_year = startyear:endyear, 
    tree = tree_name[i], 
    sens = NA, 
    sens_cil = NA,  # Lower CI
    sens_ciu = NA,  # Upper CI
    pval_lin = NA,  # p value from linear regression
    Rcor = NA,  # correlation coefficient
    pval_rcor = NA,  # p value from correlation
    # 
    TMP_ave = NA, 
    TMP_sd = NA,
    PRE_ave = NA, 
    PRE_sd = NA, 
    PET_ave = NA, 
    PET_sd = NA, 
    AI_ave = NA, 
    AI_sd = NA, 
    iWUE_ave = NA
  )
  
  for (j in 1:num_win) {
    
    # Set years
    dftemp2 <- dftemp1 %>% 
      filter(year >= (startyear - floor(winsize/2) - 1 + j) &
               year <= (startyear - floor(winsize/2) - 1 + j + winsize-1))
    
    df_regress1 <- dftemp2 %>% 
      dplyr::select(iWUE_new_iav, spei6mon_iav) %>%
      drop_na()
    
    # There should be enough data in a window
    # How to set the threshold?
    # Larger than half of the window size
    if (dim(df_regress1)[1] < winsize/2) {
      
      next
      
    } else {
      
      # Sensitivity
      regress1 <- lm(data = df_regress1, 
                     iWUE_new_iav ~ spei6mon_iav)
      regress1summ <- summary(regress1)
      regress1ci <- confint(regress1) # CI intervals
      dftemp1_winsens$sens[j] <- regress1summ$coefficients[, 1][2]
      dftemp1_winsens$pval_lin[j] <- regress1summ$coefficients[, 4][2]
      dftemp1_winsens$sens_cil[j] <- regress1ci[, 1][2]
      dftemp1_winsens$sens_ciu[j] <- regress1ci[, 2][2]
      
      # Correlation coefficient
      correlation1 <- cor.test(data = df_regress1, 
                               ~iWUE_new_iav + spei6mon_iav)
      dftemp1_winsens$Rcor[j] <- correlation1$estimate
      dftemp1_winsens$pval_rcor[j] <- correlation1$p.value
      
      # Conditions in each window
      dftemp1_winsens$TMP_ave[j] <- mean(dftemp2$tmp)
      dftemp1_winsens$PRE_ave[j] <- mean(dftemp2$pre)
      dftemp1_winsens$PET_ave[j] <- mean(dftemp2$pet)
      dftemp1_winsens$AI_ave[j] <- mean(dftemp2$ai)
      
      dftemp1_winsens$TMP_sd[j] <- sd(dftemp2$tmp)
      dftemp1_winsens$PRE_sd[j] <- sd(dftemp2$pre)
      dftemp1_winsens$PET_sd[j] <- sd(dftemp2$pet)
      dftemp1_winsens$AI_sd[j] <- sd(dftemp2$ai)
      
      # iWUE means
      dftemp1_winsens$iWUE_ave[j] <- mean(dftemp2$iWUE_new, na.rm = T)
      
    }
    
  }
  
  df_winsens <- rbind(df_winsens, dftemp1_winsens)
  
}

# Save results ------------------------------------------------------------

outfiletag_sens <- 'spei6moniav_gs3_'
write_csv(
  df_winsens,
  file = paste0(
    'Results_revise_part1/', 
    'TR_iwueiav_slr_',
    outfiletag_sens,
    'winenvi_win',
    winsize ,
    'sen.csv'
  )
)

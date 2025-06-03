# Test the parameter sensitivity in calculation iWUE

source('iWUE_setting.R')

# Load data ---------------------------------------------------------------

# carbon isotope data
infilepath1 <- 'input file paths'
infilename1 <- 'tree_ring_iwue_newco2_processed_iwue_full_iav_20nyrs_all.csv'
df_iso <- read_csv(paste0(infilepath1, infilename1))

# Process data ------------------------------------------------------------


# Subset the data frame, choose essential data for sensitivity analysis
# 
df_iso_sub <- df_iso %>% 
  dplyr::select(all_of(c('tree', 'year', 'Delta_new', 'A_ca_new', 'ca_recom', 'iWUE_new')))


# A function to test different 'b' values
ci_ca_fun <- function(delta, A_ca, ca, b) {
  
  # Set parameters for full isotope model
  
  gm = 0.2 #mesophyll conductance (mol m-2 s-1)
  f = 12 #photorespiration discrimination term (per mille)
  gammastar = 43 #photorespiratory CO2 compensation point (ppm)
  a = 4.4 #diffusive discrimination (per mille)
  #b = 30 #Rubisco discrimination (per mille)
  am = 1.8 # (per mille)
  
  ci_ca = (delta - a + (b - am)*A_ca/gm + f*gammastar/ca)/(b - a)
  return(ci_ca)
}

# Sensitivity analysis ----------------------------------------------------------------

b_para <- seq(27, 30, 0.1)
# Get mean values for the year 1951 and 2010
year1 <- 1951
year2 <- 2010
df_iso_year1 <- df_iso_sub %>% filter(year == year1)
df_iso_year2 <- df_iso_sub %>% filter(year == year2)
#
Delta_year1 <- mean(df_iso_year1$Delta_new, na.rm = T)
A_ca_year1 <- mean(df_iso_year1$A_ca_new)
ca_year1 <- mean(df_iso_year1$ca_recom)
Delta_year2 <- mean(df_iso_year2$Delta_new, na.rm = T)
A_ca_year2 <- mean(df_iso_year2$A_ca_new)
ca_year2 <- mean(df_iso_year2$ca_recom)

# Simulate the change in ci/ca
ci_ca_year1 <- ci_ca_fun(Delta_year1, A_ca_year1, ca_year1, b_para)
ci_ca_year2 <- ci_ca_fun(Delta_year2, A_ca_year2, ca_year2, b_para)


# Comparison --------------------------------------------------------------

# Data for one tree
df_iso_filteryear <- df_iso_sub %>% 
  group_by(tree) %>% 
  filter(min(year) <= year1 & max(year) >= year2)
df_iso_onetree <- df_iso_filteryear %>% 
  filter(tree == unique(df_iso_filteryear$tree)[4]) %>% 
  filter(year >= year1 & year <= year2)


# One tree data under different b
b_para_sub <- c(27, 28, 30)
df_compare <- data.frame(year = df_iso_onetree$year)
for (b_one in b_para_sub) {
  #
  df_out <- as.data.frame(ci_ca_fun(df_iso_onetree$Delta_new, df_iso_onetree$A_ca_new, df_iso_onetree$ca_recom, b_one))
  colnames(df_out) <- as.character(paste0('b_', b_one))
  df_compare <- df_compare %>% bind_cols(df_out)
}

# Make a long data frame
df_compare_long <- df_compare %>% 
  pivot_longer(names_to = 'b_para', values_to = 'b_value', 
               cols = c('b_27', 'b_28', 'b_30'))


# Plotting ----------------------------------------------------------------

###
# Change in ci/ca
df_year1 <- data.frame(
  b = b_para, 
  cica = ci_ca_year1, 
  cica_rel = (ci_ca_year1-ci_ca_year1[1]) / ci_ca_year1[1]*100
)
df_year2 <- data.frame(
  b = b_para, 
  cica = ci_ca_year2, 
  cica_rel = (ci_ca_year2-ci_ca_year2[1]) / ci_ca_year2[1]*100
)
#plot(x = b_para, y = ci_ca_year1)
f_year1 <- ggplot(data = df_year1, aes(x = b, y = cica_rel)) + 
  geom_line() + 
  labs(
    x = expression(italic(b)~'(\u2030)'), 
    y = 'Relative change (%)', 
    tag = 'a', 
    title = paste0('Year ', year1)
  ) + 
  scale_y_continuous(limits = c(-10, 0)) + 
  theme_char + 
  geom_hline(yintercept = 0, linetype = 'dashed')
# f_year1

f_year2 <- ggplot(data = df_year2, aes(x = b, y = cica_rel)) + 
  geom_line() + 
  labs(
    x = expression(italic(b)~'(\u2030)'), 
    y = 'Relative change (%)', 
    tag = 'b', 
    title = paste0('Year ', year2)
  ) + 
  scale_y_continuous(limits = c(-10, 0)) + 
  theme_char + 
  geom_hline(yintercept = 0, linetype = 'dashed')

f_years_combine <- f_year1 + f_year2

###

# Time series of data
ggplot(data = df_iso_onetree, aes(x = year, y = iWUE_new)) + 
  geom_point()

# Relationship between values under different b
f_compare1 <- ggplot(data = df_compare, 
       aes(x = b_27, y = b_30)) + 
  geom_point(
    alpha = 0.7,
    stroke = 0.6, 
    size = 6*5/14
  ) + 
  geom_smooth(linewidth = 0.6) + 
  labs(
    x = expression(italic(b)~'='~'27\u2030'), 
    y = expression(italic(b)~'='~'30\u2030'), 
    tag = 'a', 
    title = expression('Comparison of '*italic(c)[italic(i)]*'/'*italic(c)[italic(a)])
  ) + 
  scale_y_continuous(breaks = seq(0.68, 0.76, 0.02)) + 
  theme_char

f_compare2 <- ggplot(data = df_compare, 
       aes(x = b_28, y = b_30)) + 
  geom_point(
    alpha = 0.7,
    stroke = 0.6, 
    size = 6*5/14
  ) + 
  geom_smooth(linewidth = 0.6) + 
  labs(
    x = expression(italic(b)~'='~'28\u2030'), 
    y = expression(italic(b)~'='~'30\u2030'), 
    tag = 'b', 
    title = expression('Comparison of '*italic(c)[italic(i)]*'/'*italic(c)[italic(a)])
  ) + 
  scale_y_continuous(breaks = seq(0.68, 0.76, 0.02)) + 
  theme_char

f_compare_combine <- f_compare1 + f_compare2

###
# Examine other variables
ggplot(data = df_iso_sub, aes(x = year, y = A_ca_new)) + 
  geom_point()


# Save data ---------------------------------------------------------------

# To save data indicating relative changes
df_year1_ylab <- df_year1 %>% add_column(year = 1951)
df_year2_ylab <- df_year2 %>% add_column(year = 2010)
df_parab_rel <- bind_rows(df_year1_ylab, df_year2_ylab)
data_outpath = "out file paths"
write_csv(df_parab_rel, 
          file = paste0(data_outpath, 'SuppFig8/', 'iwue_para_b_supp_fig8.csv'))
write_csv(df_compare, 
          file = paste0(data_outpath, 'SuppFig9/', 'iwue_para_b_supp_fig9.csv'))


# Save figures ------------------------------------------------------------

ggsave(
  filename = paste0(infilepath1, 'Figures/Response_cica_calculation_relchange.jpeg'), 
  f_years_combine, 
  units = 'mm',
  dpi = 900,
  width = 120,
  height = 65
)

ggsave(
  filename = paste0(infilepath1, 'Figures/Response_cica_calculation_compare.jpeg'), 
  f_compare_combine, 
  units = 'mm',
  dpi = 900,
  width = 120,
  height = 65
)

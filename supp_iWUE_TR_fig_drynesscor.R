# Relationships between IAV of SM/VPD and SPEI

source('/iWUE_setting.R')
infilepath1 <- 'Results_revise_part1/'
outfilepath1 <- infilepath1
world_map <- map_data('world')

# Load data ---------------------------------------------------------------

df <- read_csv(paste0(infilepath1, 'TR_iwueiav_v3.csv'))

# Prepare data ------------------------------------------------------------

cor_fun <- function(x1, x2) {
  # Calculate the correlation between two dryness variables
  dfcor <- data.frame(x = x1, y = x2)
  dfcor <- dfcor %>% drop_na()
  if (dim(dfcor)[1] < 5) {
    # Not enough data!
    Rcor <- NA
    pval <- NA
    return(c(Rcor, pval))
  }
  correlation <- cor.test(x1, x2)
  Rcor = correlation$estimate 
  pval = correlation$p.value
  return(c(Rcor, pval))
}

# Mean of sites --------
start_period_year <- 1951
end_period_year <- 2010
df_num <- df %>%
  filter(year >= start_period_year &
           year <= end_period_year) %>%
  drop_na(SoilMoi0_10cm_inst, gs_vpd, spei6mon_mean) %>%
  group_by(tree) %>%
  summarise(num = n())
  # Equal sample size
  # filter(num > 40)

# Note: cannot make sure enough sample size (6?) during the whole period
df_mean <- df %>%
  filter(year >= start_period_year &
           year <= end_period_year) %>%
  filter(tree %in% df_num$tree) %>%
  drop_na(SoilMoi0_10cm_inst, gs_vpd, spei6mon_mean) %>%
  group_by(year) %>%
  summarise(
    sm0_10cm_mean = mean(SoilMoi0_10cm_inst, na.rm = T),
    vpd_mean = mean(gs_vpd, na.rm = T),
    spei_mean = mean(spei6mon_mean, na.rm = T),
  )

# Estimate correlation at each tree --------
df_smcor <- df %>% 
  filter(year >= start_period_year &
           year <= end_period_year) %>%
  group_by(tree) %>% 
  summarize(
    rcor = cor_fun(SoilMoi0_10cm_inst, spei6mon_mean)[1], 
    pval = cor_fun(SoilMoi0_10cm_inst, spei6mon_mean)[2]
  )

df_vpdcor <- df %>% 
  group_by(tree) %>% 
  filter(year >= start_period_year &
           year <= end_period_year) %>%
  summarize(
    rcor = cor_fun(gs_vpd, spei6mon_mean)[1], 
    pval = cor_fun(gs_vpd, spei6mon_mean)[2]
  )

df_vpdsm <- df %>% 
  group_by(tree) %>% 
  filter(year >= start_period_year &
           year <= end_period_year) %>%
  summarize(
    rcor = cor_fun(gs_vpd, SoilMoi0_10cm_inst)[1], 
    pval = cor_fun(gs_vpd, SoilMoi0_10cm_inst)[2]
  )

# IAV
df_smiavcor <- df %>% 
  group_by(tree) %>% 
  summarize(
    rcor = cor_fun(sm0_10cm_iav, spei6mon_iav)[1], 
    pval = cor_fun(sm0_10cm_iav, spei6mon_iav)[2]
  )

df_vpdiavcor <- df %>% 
  group_by(tree) %>% 
  summarize(
    rcor = cor_fun(gs_vpd_iav, spei6mon_iav)[1], 
    pval = cor_fun(gs_vpd_iav, spei6mon_iav)[2]
  )

# Analysis ----------------------------------------------------------------

mean(df_smcor$rcor, na.rm = T)
median(df_smcor$rcor, na.rm = T)
t.test(df_smcor$rcor)

# Relationship Figures -------------------------------------------------------------------------

# Density figures of R --------
color_dens <- '#636363'
# VPD
f_vpd1 <- ggplot(data = df_vpdcor %>% filter(tree %in% df_num$tree)) + 
  geom_density(
    aes(x = rcor), 
    fill = alpha(color_dens, 0.5), 
    color = color_dens, 
    position = 'identity'
  ) + 
  labs(x = expression(italic(R)~"between SPEI and "*italic(D)), 
       y = "Density") + 
  geom_vline(
    xintercept = median(df_vpdcor$rcor, na.rm = T), 
    linetype = 'dashed'
  )+
  theme_char + 
  theme(
    axis.title = element_text(size = 7.5)
  )

# SM
f_sm1 <- ggplot(data = df_smcor %>% filter(tree %in% df_num$tree)) + 
  geom_density(
    aes(x = rcor), 
    fill = alpha(color_dens, 0.5), 
    color = color_dens, 
    position = 'identity'
  ) + 
  labs(x = expression(italic(R)~"between SPEI and SM"), 
       y = "Density", 
       tag = 'b') + 
  geom_vline(
    xintercept = median(df_smcor$rcor, na.rm = T), 
    linetype = 'dashed'
  )+
  theme_char + 
  theme(
    axis.title = element_text(size = 7.5)
  )

# VPD vs SM
f_vpdsm1 <- ggplot(data = df_vpdsm %>% filter(tree %in% df_num$tree)) + 
  geom_density(
    aes(x = rcor), 
    fill = alpha(color_dens, 0.5), 
    color = color_dens, 
    position = 'identity'
  ) + 
  labs(x = expression(italic(R)~"between SM and "*italic(D)), 
       y = "Density", 
       tag = 'd') + 
  geom_vline(
    xintercept = median(df_vpdsm$rcor, na.rm = T), 
    linetype = 'dashed'
  )+
  theme_char + 
  theme(
    axis.title = element_text(size = 7.5)
  )

# Scattering ----------
var_point_size <- 5
var_text_size <- 2.5

if (F) {
# VPD
corvpd <- cor.test(data = df_mean, ~ spei_mean + vpd_mean)
f_vpd2 <- ggplot(data = df_mean,
       aes(x = vpd_mean,
           y = spei_mean)) +
  geom_point(size = var_point_size * 5 / 14) +
  labs(x = expression(italic(D) ~ "(kPa)"),
       y = 'SPEI',
       tag = 'a') +
  geom_smooth(method = 'lm',
              color = 'black',
              linetype = 'dashed') +
  geom_text(
    x = 6,
    y = -0.25,
    label = paste('italic(R)', '==', round(corvpd$estimate, 2),
                  sep = ''),
    size = var_text_size,
    parse = T
  ) +
  geom_text(
    x = 6,
    y = -0.3,
    label = paste('italic(p)', '<~0.001',
                  sep = ''),
    size = var_text_size,
    parse = T
  ) +
  theme_char 
}

# SM
corsm <- cor.test(data = df_mean, ~ spei_mean + sm0_10cm_mean)
f_sm2 <- ggplot(data = df_mean,
       aes(x = sm0_10cm_mean,
           y = spei_mean)) +
  geom_point(size = var_point_size * 5 / 14) +
  labs(x = expression(SM ~ "(kg "*m^-2*")"),
       y = 'SPEI',
       tag = 'a') +
  geom_smooth(method = 'lm',
              color = 'black',
              linetype = 'dashed') +
  geom_text(
    x = 26,
    y = -0.27,
    label = paste('italic(R)', '==', round(corsm$estimate, 2),
                  sep = ''),
    size = var_text_size,
    parse = T
  ) +
  geom_text(
    x = 26,
    y = -0.32,
    label = paste('italic(p)', '<~0.01',
                  sep = ''),
    size = var_text_size,
    parse = T
  ) +
  theme_char 

# SM vs VPD
corvpdsm <- cor.test(data = df_mean, ~ sm0_10cm_mean + vpd_mean)
f_vpdsm2 <- ggplot(data = df_mean,
                 aes(x = vpd_mean,
                     y = sm0_10cm_mean)) +
  geom_point(size = var_point_size * 5 / 14) +
  labs(x = expression(italic(D) ~ "(kPa)"),
       y = expression(SM ~ "(kg "*m^-2*")"),
       tag = 'c') +
  geom_smooth(method = 'lm',
              color = 'black',
              linetype = 'dashed') +
  geom_text(
    x =5.8,
    y = 25,
    label = paste('italic(R)', '==', round(corvpdsm$estimate, 2),
                  sep = ''),
    size = var_text_size,
    parse = T
  ) +
  geom_text(
    x =5.8,
    y = 24.9,
    label = paste('italic(p)', '<~0.05',
                  sep = ''),
    size = var_text_size,
    parse = T
  ) +
  theme_char 

# Combine figures ---------------------------------------------------------

f_combine1 <- (f_sm2+f_sm1) / (f_vpdsm2 + f_vpdsm1)
f_combine2 <- f_sm2+f_sm1

# Save figures ------------------------------------------------------------

ggsave(
  f_combine1, 
  filename = paste0(
    outfilepath1,
    'Figures/',
    'supp_dryness_cor2.jpeg'
  ),
  units = 'mm',
  dpi = 900,
  width = 120,
  height = 120
)

ggsave(
  f_combine2, 
  filename = paste0(
    outfilepath1,
    'Figures/',
    'supp_dryness_cor2_v2.jpeg'
  ),
  units = 'mm',
  dpi = 900,
  width = 120,
  height = 60
)

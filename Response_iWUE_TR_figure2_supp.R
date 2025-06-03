# Use growing season mean of climates (MATGS, MAPGS)
# to examine the relationship between SW and climatic conditions
# Supplements for Figure 2

rm(list = ls())

source('iWUE_setting.R')
library(tidyverse)

# Load data ---------------------------------------------------------------


inputpath <- '/Users/mw/Documents/Proj_Postdoc/Proj_WUETrait/Revision_202212/Results_revise_part1/'
 
# Growing-season annual data
df_cru = read_csv(paste0(inputpath, 'TRMET_data_cru_climate_gs3.csv'))

# Original SW; iWUE sensitivity
df_sw = read_csv(paste0(inputpath, 'TR_iwue_sensdiff_win22sens_1951_2010.csv'))
# Use alternative data
# df_sw = read_csv(paste0('/Users/mw/Documents/Proj_Postdoc/Proj_WUETrait/Proj_WUETrait_back/Results/', 
#                  'TRtree_climate_condition_geoinfo_cru_1951_2015.csv'))


# Growing season mean for each tree ---------------------------------------

tree_id = unique(df_cru$tree)
df_cru_gsmean <- df_cru %>% 
  # Note: the years are not consistent with those for calculating total rain
  # filter(year >= 1971 & year <= 2015) %>% 
  filter(year >= 1951 & year <=2015 ) %>% 
  group_by(tree) %>% 
  summarize(mean_tmp = mean(tmp, na.rm = T), 
            mean_pre = mean(pre, na.rm = T))

df_sw_cru_gsmean <- df_sw %>% 
  left_join(df_cru_gsmean, by = 'tree')


# Statistics --------------------------------------------------------------

cor.test(x = df_sw_cru_gsmean$mat, y = df_sw_cru_gsmean$mean_tmp)
cor.test(x = df_sw_cru_gsmean$map, y = df_sw_cru_gsmean$mean_pre)

# changes in SW and climates
# cor.test(x = df_sw_cru_gsmean$diff, y = df_sw_cru_gsmean$mat)
# cor.test(x = df_sw_cru_gsmean$diff, y = df_sw_cru_gsmean$map)
# 
# cor.test(x = df_sw_cru_gsmean$diff, y = df_sw_cru_gsmean$mean_tmp)
# cor.test(x = df_sw_cru_gsmean$diff, y = df_sw_cru_gsmean$mean_pre)


# Visualization -----------------------------------------------------------

# Correlation for climate variables ---------------------------------------


# Temperature
# For text
tmp_res <- cor.test(x = df_sw_cru_gsmean$mat, y = df_sw_cru_gsmean$mean_tmp)
# 
p_value = round(tmp_res$p.value, 3)
if (p_value != 0) {
  pvalue_text = paste('italic(p)', '==', p_value, sep = '')
} else if (p_value == 0) {
  pvalue_text <- paste('italic(p)', '<', 0.001, sep = '')
}
# Correlation coefficient
R <- round(tmp_res$estimate, 2)
R_text <- paste('italic(R)', '==', R, sep = '')

# Plotting
fig_tmp <- ggplot(data = df_sw_cru_gsmean,
       aes(x = mat,
           y = mean_tmp)) + 
  geom_point(
    alpha = 0.7,
    stroke = 0.6, 
    size = 6*5/14
  ) + 
  geom_smooth(method = 'lm', formula = y ~ x) + 
  # Set axises
  labs(tag = 'a', 
       x = expression('MAT ('*degree*C*')'), 
       y = expression(MAT[GS]~'('*degree*C*')')) + 
  lims(x = c(-15, 30),  
       y = c(-5, 30)) + 
  theme_bw() + 
  theme_char + 
  # Show texts
  annotate("text", x = min(df_sw_cru_gsmean$mat), y = max(df_sw_cru_gsmean$mean_tmp), 
           label = R_text, size = 7*5/14, 
           hjust = 0, vjust = 0.5, parse = T) + 
  annotate("text", x = min(df_sw_cru_gsmean$mat), y = max(df_sw_cru_gsmean$mean_tmp), 
           label = pvalue_text, size = 7*5/14, 
           hjust = 0, vjust = 2, parse = T)


# Precipitation
# Note: growing season length: 4-10, 7 months
pre_res <- cor.test(x = df_sw_cru_gsmean$map, y = df_sw_cru_gsmean$mean_pre)
# 
p_value = round(pre_res$p.value, 3)
if (p_value != 0) {
  pvalue_text = paste('italic(p)', '==', p_value, sep = '')
} else if (p_value == 0) {
  pvalue_text <- paste('italic(p)', '<', 0.001, sep = '')
}
# Correlation coefficient
R <- round(pre_res$estimate, 2)
R_text <- paste('italic(R)', '==', R, sep = '')

fig_pre <- ggplot(data = df_sw_cru_gsmean,
       aes(x = map*12,
           y = mean_pre*7)) + 
  geom_point(
    alpha = 0.7,
    stroke = 0.6, 
    size = 6*5/14
  ) + 
  geom_smooth(method = 'lm', formula = y ~ x) + 
  labs(tag = 'b', 
       x = 'MAP (mm)', 
       y = expression(MAP[GS]~'(mm)')) + 
  theme_bw() + 
  theme_char + 
  # Show texts
  annotate("text", x = min(df_sw_cru_gsmean$map*12), y = max(df_sw_cru_gsmean$mean_pre*7), 
           label = R_text, size = 7*5/14, 
           hjust = 0, vjust = 1.5, parse = T) + 
  annotate("text", x = min(df_sw_cru_gsmean$map*12), y = max(df_sw_cru_gsmean$mean_pre*7), 
           label = pvalue_text, size = 7*5/14, 
           hjust = 0, vjust = 3, parse = T)
fig_pre


# Combine figures

fig_combine <- fig_tmp + fig_pre


# Scattering climate space--------

var_fill_direction <- 1
var_palette <- 'RdYlBu'
var_fill_type <- 'seq'
var_point_size <- 6
f_climspace1 <- ggplot(data = df_sw_cru_gsmean) + 
  geom_point(
    aes(
      x = mean_tmp, 
      y = mean_pre*7,  # Growing season: April - October 
      color = diff, 
      shape = group
    ), 
    alpha = 0.7,
    stroke = 0.6, 
    size = var_point_size*5/14
  ) + 
  labs(
    x = expression(MAT[GS]~'('*degree*C*')'), 
    y = expression(MAP[GS]~'(mm)'), 
    fill = NULL, 
    color = NULL, 
    # tag = 'c'
  ) + 
  scale_shape_manual(
    guide = 'none', 
    values = c('a' = 2, 'g' = 1)
  ) + 
  scale_color_fermenter(
    direction = var_fill_direction, 
    palette = var_palette, 
    type = var_fill_type, 
    breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)
  ) + 
  scale_y_continuous(
    breaks = c(500, 1000, 1500, 2000, 2500, 3000, 3500)
  )+
  theme_char + 
  theme(
    # legend.key.height= unit(3.5, 'mm'),
    legend.key.width= unit(3.5, 'mm')
  )
f_climspace1

# Save figures ------------------------------------------------------------

outfigpath <- inputpath

ggsave(
  fig_combine, 
  filename = paste0(outfigpath, 'Figures/Response_climate_relation_all_gs.jpeg'),
  units = 'mm',
  dpi = 900,
  width = 120,
  height = 65
)

# Scattering
ggsave(
  f_climspace1, 
  filename = paste0(outfigpath, 'Figures/Response_climate_space_gs.jpeg'),
  units = 'mm',
  dpi = 900,
  width = 70,
  height = 55
)

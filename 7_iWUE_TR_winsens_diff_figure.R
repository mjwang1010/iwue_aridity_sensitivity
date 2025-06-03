# Figure 2, difference of SW between the two periods
# Spatial patterns
# Density distribution
# Scattering

source('iWUE_setting.R')
# Path of the input file folder
infilepath1 <- '...'
# Path of the output file folder
outfilepath1 <- '...'
world_map <- map_data('world')

# Load data ---------------------------------------------------------------

df <- read_csv(paste0(infilepath1, 'TR_iwue_sensdiff_win22sens_1951_2010.csv'))

# Prepare data ------------------------------------------------------------


# Statistical analysis ----------------------------------------------------------------

dim(df %>% filter(group == 'a'))[1]
dim(df %>% filter(group == 'g'))[1]

min(df$diff)
max(df$diff)
# Negative change in SW
dim(df %>% filter(diff < 0))[1]
dim(df %>% filter(diff < 0&aidiff < 0))[1]
dim(df %>% filter(diff < 0&aidiff >= 0))[1]
dim(df %>% filter(diff >= 0))[1]
dim(df %>% filter(diff >= 0&aidiff < 0))[1]
dim(df %>% filter(diff >= 0&aidiff >= 0))[1]

cor.test(data = df, 
         ~diff+mat)
cor.test(data = df, 
         ~diff+map)
cor.test(data = df, 
         ~diff+maa)
t.test(df$diff)

# Angiosperms
df_ang <- df %>% 
  filter(group == 'a')
df_ang_dry <- df %>% 
  filter(group == 'a') %>% 
  filter(aidiff < 0 & diff < 0 )
hist(df_ang$diff)
cor.test(df_ang$aidiff, df_ang$diff)
t.test(df_ang$diff)

# Gymnosperms
df_gym <- df %>% filter(group == 'g')
df_gym1 <- df %>% filter(group == 'g') %>% 
  filter(maa <= quantile(df_gym$maa, 0.5, type = 3))
df_gym2 <- df %>% filter(group == 'g') %>% 
  filter(maa > quantile(df_gym$maa, 0.5, type = 3))
df_gym_dry <- df %>% filter(group == 'g') %>% 
  filter(aidiff < 0 & diff < 0 )
hist(df_gym1$diff)
hist(df_gym2$diff)
plot(df_gym$aidiff, df_gym$diff)
cor(df_gym2$aidiff, df_gym2$diff)
t.test(df_gym1$diff, df_gym2$diff)
t.test(df_gym1$diff)
t.test(df_gym2$diff)

# Estimate relative change for gymnosperms
df_gym <- df %>% filter(group == 'g')
df_gym <- df_gym %>% mutate(reldiff = diff/Early_period)
df_gym_earlyneg <- df_gym %>% filter(Early_period <= 0)
reldiff_gym_mean <- mean(df_gym$reldiff, na.rm = T)

mean(df_gym_earlyneg$diff)
mean(df_gym_earlyneg$reldiff)
median(df_gym_earlyneg$diff)
median(df_gym_earlyneg$reldiff)

summary(df_gym_earlyneg$reldiff)

# 
ggplot(data = df_gym_earlyneg) + 
  geom_density(aes(x = diff))

# Figures -----------------------------------------------------------------

# Spatial map of tree sites --------

geo_map_fun1 <- function(df1, var1) {
  
 if (var1 == 'diff') {
    
    var_title <- expression(atop(paste('Change in '*S[italic(W)]),
                            paste('(',mu,'mol ',mol^-1,' per unit SPEI)')))
    var_fill_direction <- 1
    var_palette <- 'RdYlBu'
    var_fill_type <- 'div'
    var_showlimit <- F
    # Consistent with the scattering plot
    var_breaks <- c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)
    # var_breaks <- waiver()
    
  } else if (var1 == 'maa') {
    
    var_title <- 'Aridity index from 1971 to 2015'
    var_fill_direction <- 1
    # var_palette <- "YlOrRd"
    var_palette <- 'RdYlBu'
    var_fill_type <- 'div'
    var_showlimit <- F
    var_breaks <- c(0.05, 0.2, 0.5, 0.65, 1, 2, 3) # Refer to online information
    
  }
  
  map1 <- ggplot() + 
    geom_polygon(
      data = world_map,
      aes(x = long, y = lat, group = group),
      color = NA,
      fill = '#d9d9d9',
      linewidth = 0.05
    ) + 
    geom_point(
      data = df1,
      aes(
        x = lon,
        y = lat,
        shape = group, 
        fill = get(as.name(var1))
      ), 
      size = 4 * 5 / 14, 
      stroke = 0.1, 
      alpha = 0.6
    ) + 
    scale_x_continuous(breaks = seq(-180, 180, 60)) + 
    scale_shape_manual(
      guide = 'none', 
      values = c(24, 21)
    ) + 
    scale_fill_fermenter(
      breaks = var_breaks,
      type = var_fill_type, 
      palette = var_palette,
      direction = var_fill_direction
    ) + 
    guides(
      fill = guide_colorsteps(
        show.limits = var_showlimit, 
        ticks = T, 
        ticks.colour = 'black', 
        barwidth = unit(45, 'mm'), # original 60
        barheight = unit(5, 'points'), 
        title = var_title
      )
    ) + 
    labs(
      x = paste('Longitude (\u00b0E)'), 
      y = paste('Latitude (\u00b0N)')
    ) +  
    coord_fixed() + 
    theme_bw() + 
    theme(
      panel.grid = element_blank(), 
      text = element_text(size = 7), 
      legend.title = element_text(size = 8),
      legend.position = 'bottom', 
      legend.box.spacing = unit(0, 'mm'), 
      plot.title = element_text(hjust = 0.5, size = 8)
    )
  return(map1)
}
f_sensdiff_map <- geo_map_fun1(df, 'diff') + 
  labs(tag = 'a')
f_sensdiff_map

# Density distribution of SW differences-------

color_sensdiff <- '#5e3c99'
# Which clades
df_in <- df %>% 
  filter(group == 'g')
df_in <- df

f_diffdens1 <- ggplot(data = df_in) + 
  geom_density(
    aes(x = diff), 
    fill = alpha(color_sensdiff, 0.5), 
    color = color_sensdiff, 
    position = 'identity'
  ) + 
  labs(x = expression("Change in "*S[italic(W)]), 
       y = "Density", 
       tag = 'b') + 
  geom_vline(
    xintercept = 0, 
    linetype = 'dashed'
  )+
  theme_char + 
  theme(
    axis.title = element_text(size = 8.5)
  )
f_diffdens1


# Scattering--------

var_fill_direction <- 1
var_palette <- 'RdYlBu'
var_fill_type <- 'seq'
var_point_size <- 6
f_climspace1 <- ggplot(data = df) + 
  geom_point(
    aes(
      x = mat, 
      y = map*12, 
      color = diff, 
      shape = group
    ), 
    alpha = 0.7,
    stroke = 0.6, 
    size = var_point_size*5/14
  ) + 
  labs(
    x = expression('MAT ('*degree*C*')'), 
    y = 'MAP (mm)', 
    fill = NULL, 
    color = NULL, 
    tag = 'c'
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
    breaks = c(500, 1000, 1500, 2000, 2500)
  )+
  theme_char + 
  theme(
    # legend.key.height= unit(3.5, 'mm'),
    legend.key.width= unit(3.5, 'mm')
  )
f_climspace1

# Climate feature space, by climate difference
f_climspace2 <- ggplot(data = df) + 
  geom_point(
    aes(
      x = tmpdiff, 
      # y = prediff*12,
      y = aidiff,
      fill = diff, 
      shape = group
    ), 
    alpha = 0.6, 
    size = var_point_size*5/14
  ) + 
  labs(
    x = expression('Change in MAT ('*degree*C*')'), 
    # y = 'Change in MAP (mm)',
    y = 'Change in AI',
    fill = NULL, 
    color = NULL
  ) + 
  scale_shape_manual(
    guide = 'none', 
    values = c('a' = 24, 'g' = 21)
  ) + 
  scale_fill_fermenter(
    direction = var_fill_direction, 
    palette = var_palette, 
    type = var_fill_type, 
    breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5)
  ) + 
  theme_char + 
  theme(
    # legend.key.height= unit(3.5, 'mm'),
    legend.key.width= unit(3.5, 'mm')
  ) + 
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  geom_hline(yintercept = 0, linetype = 'dashed')
f_climspace2


# Combine figures ---------------------------------------------------------

figure2_v1 <- f_sensdiff_map / (f_diffdens1 | f_climspace1) + 
  plot_layout(heights = c(2, 1.5))

# Save figures ------------------------------------------------------------

ggsave(
  figure2_v1,
  filename = paste0(
    outfilepath1,
    'Figures/',
    'fig2_iwue_spei6mon_gs3_sensdiff_v2.jpeg'
  ),
  units = 'mm',
  dpi = 900,
  width = 120,
  height = 120
)

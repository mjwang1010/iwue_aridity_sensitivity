# Figure3, Potential factors controlling SW changes
# AI change

source('Path for iWUE_setting.R')
infilepath1 <- 'Path of the input file folder'
outfilepath1 <- infilepath1

# Load data ---------------------------------------------------------------

df1 <- read_csv(paste0(infilepath1, 'TR_iwue_sensdiff_win22sens_1951_2010.csv'))

# Analysis ----------------------------------------------------------------

max(df1$diff)
min(df1$diff)
max(df1$aidiff)
min(df1$aidiff)

# Fraction
# Gymnosperms
df1_g <- df1 %>% filter(group == 'g')
dim(df1_g %>% filter(diff<0))[1]
dim(df1_g %>% filter(diff<0&aidiff<0))[1]/dim(df1_g)[1]
dim(df1_g %>% filter(diff>0&aidiff<0))[1]/dim(df1_g)[1]
dim(df1_g %>% filter(diff<0&aidiff>0))[1]/dim(df1_g)[1]
dim(df1_g %>% filter(diff>0&aidiff>0))[1]/dim(df1_g)[1]

t.test(df1_g$diff)
t.test(df1_g$aidiff)
dim(df1_g %>% filter(aidiff < 0))[1]

# Angiosperms
df1_a <- df1 %>% filter(group == 'a')
dim(df1_a %>% filter(diff<0))[1]
dim(df1_a %>% filter(diff<0&aidiff<0))[1]/dim(df1_a)[1]
dim(df1_a %>% filter(diff>0&aidiff<0))[1]/dim(df1_a)[1]
dim(df1_a %>% filter(diff<0&aidiff>0))[1]/dim(df1_a)[1]
dim(df1_a %>% filter(diff>0&aidiff>0))[1]/dim(df1_a)[1]

t.test(df1_a$diff)
t.test(df1_a$aidiff)
dim(df1_a %>% filter(aidiff < 0))[1]

# Correlation between changes in SW and AI
cor.test(x = df1_g$diff, 
         y = df1_g$aidiff)

# Figures -----------------------------------------------------------------

# Scattering--------
df1_label <- df1 %>% 
  mutate(
    label = ifelse(diff>=0&aidiff<0, 'type1', 
                   ifelse(diff>=0&aidiff>=0, 'type2', 
                          ifelse(diff<0&aidiff<0, 'type3', 
                                 ifelse(diff<0&aidiff>=0, 'type4', NA))))
  )

lim_x <- c(-0.45, 0.25)
lim_y <- c(-7.5, 12)
break_x <- seq(-0.4, 0.2, 0.1)
break_y <- seq(-6, 10, 2)
point_size <- 5 * 5/14

fig_diff_factor_fun <- function(df, clade = 'g') {
  # Change in AI against change in SW
  
  if (clade == 'g') {
    shape_char <- 21
    tag_char <- 'c'
  } else if (clade == 'a') {
    shape_char <- 24
    tag_char <- 'b'
  }
  
  f <- ggplot(data = df %>% filter(group == clade)) + 
    geom_point(
      aes(x = aidiff, 
          y = diff, 
          color = label), 
      shape = shape_char, 
      size = point_size, 
      show.legend = F
    ) + 
    labs(
      x = 'Change in AI', 
      # y = expression('Change in '*S[italic(W)]), 
      y = expression(atop(paste('Change in '*S[italic(W)]),
                          paste('(',mu,'mol ',mol^-1,' per unit SPEI)'))), 
      tag = tag_char
    ) + 
    scale_x_continuous(
      limits = lim_x, 
      breaks = break_x
    ) + 
    scale_y_continuous(
      limits = lim_y, 
      breaks = break_y
    ) + 
    scale_color_manual(
      values = c('type1' = '#018571', 
                 'type2' = '#80cdc1', 
                 'type3' = '#a6611a', 
                 'type4' = '#dfc27d')
    ) + 
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    geom_hline(yintercept = 0, linetype = 'dashed') + 
    theme_char
  
  return(f)
}

f_diff_gym <- fig_diff_factor_fun(df1_label, 'g')
f_diff_gym

f_diff_ang <- fig_diff_factor_fun(df1_label, 'a')
f_diff_ang

# Box plot -----------
# Change in SW
color_ang <- '#d01c8b'
color_gym <- '#4dac26'

f_box1 <- ggplot(data = df1, 
       aes(x = group, y = diff)) + 
  geom_boxplot(
    aes(fill = group, 
        color = group),
    width = 0.2, 
    outlier.size = 1.5 * 5 / 14
  ) + 
  stat_summary(
    aes(group = group, 
        color = group),
    fun = mean,
    size = 3.5*5/14, 
    shape  = 2, 
    geom = 'point'
  ) + 
  scale_fill_manual(
    values = alpha(c(color_ang, color_gym), 0.5),
    guide = 'none'
  ) + 
  scale_x_discrete(labels = c('Angiosperms', 'Gymnosperms')) + 
  scale_color_manual(values = c(color_ang, color_gym), 
                     guide = 'none') + 
  labs(
    tag = 'a', 
    x = NULL,
    # y = expression('Change in '*S[italic(W)])
    y = expression(atop(paste('Change in '*S[italic(W)]),
                    paste('(',mu,'mol ',mol^-1,' per unit SPEI)')))) + 
  # Add mean points
  theme_bw() + 
  theme_char + 
  theme(axis.text = element_text(size = 5)) + 
  geom_hline(yintercept = 0, linetype = 'dashed')

# Combine figures ---------------------------------------------------------

f_combine1 <- f_box1 + f_diff_ang + f_diff_gym + 
  plot_layout(widths = c(1, 1.5, 1.5))

# Save figures ------------------------------------------------------------

ggsave(
  f_combine1,
  filename = paste0(
    outfilepath1,
    'Figures/',
    'fig3_iwue_sensdiff_factor_aidiff.jpeg'
  ),
  units = 'mm',
  dpi = 900,
  width = 160,
  height = 50
)


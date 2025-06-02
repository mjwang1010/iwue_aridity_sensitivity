# CMIP6 results analysis
# Data visualization and statistical analysis


source('/Users/mw/Documents/Proj_Postdoc/Proj_WUETrait/Revision_202212/iWUE_setting.R')


# Load results ------------------------------------------------------------

infilepath = '/Users/mw/Documents/Proj_Postdoc/Proj_WUETrait/Revision_202212/Results_revise_part1/CMIPOutput/Postprocess/'

co2type = 'meanco2'

filenamelist = list.files(path = infilepath, pattern = paste0('.*', co2type))
# Read all SW from different scenarios
df_sim <- data.frame()
for (filename in filenamelist) {
  df_temp = read_csv(paste0(infilepath, filename))
  df_sim <- bind_rows(df_sim, df_temp)
}

if (FALSE) {
  # Save the table
  write_csv(
    df_sim, 
    file = paste0(infilepath, 'model_sim_sens_multiscen.csv')
  )
}


# Data visualization ------------------------------------------------------

# Create a box plot
colors <- c('gray', '#2c7bb6', '#abd9e9', '#ffffbf', '#fdae61')
mean_historical <- mean(df_sim[df_sim$experiment=='historical', 'sw']$sw, na.rm = T)
plot <- ggplot(data = df_sim, aes(x = experiment, y = sw, fill = experiment)) +
  geom_boxplot(color = "black", 
               #outlier.shape = NA, 
               width = 0.4) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 0.75) + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(limits = c(-4, -1)) + 
  labs(x = NULL, 
       # y: soil moisture sensitivity
       y = expression(atop(paste(S[italic(W)]),
                           #paste('(',mu,'mol ',mol^-1,' per unit SPEI)')
                           )), 
       title = expression("Comparisons of "~S[italic(W)]~"in Different Scenarios")) + 
  geom_hline(yintercept = mean_historical, linetype = "dashed", color = "black") +
  theme_char + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "none", 
    plot.title = element_text(size = 6.5)
  )

plot

# Save the figure with specified dimensions
outfigpath <- '/Users/mw/Documents/Proj_Postdoc/Proj_WUETrait/Revision_202212/Results_revise_part1/Figures/'
ggsave(filename = paste0(outfigpath, 'fig_iwue_model_smsens_future_scen_v3.jpeg'), 
       plot, 
       dpi = 600, 
       width = 80, 
       height = 60, 
       units = "mm")

# Statistical analysis ----------------------------------------------------

# Compare historical to future scenarios
# hist vs. ssp126
df_test <- df_sim %>% filter(experiment %in% c('historical', 'ssp126'))
t.test(formula = sw~experiment, data = df_test)

df_sim %>% 
  group_by(experiment) %>% 
  summarize(swmean = mean(sw, na.rm = T), 
            n_mod = n())

df_sim %>% 
  group_by(experiment) %>% 
  summarize(
    swmean = mean(sw, na.rm = T), 
    n_mod = n(), 
    swsd = sd(sw, na.rm = T), 
    swse = sd(sw, na.rm = T)/n_mod
  )

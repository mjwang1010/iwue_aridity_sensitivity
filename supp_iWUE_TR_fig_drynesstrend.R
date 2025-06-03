# The dryness trends

source('/iWUE_setting.R')
infilepath1 <- 'Results_revise_part1/'
outfilepath1 <- infilepath1
world_map <- map_data('world')

# Load data ---------------------------------------------------------------

df <- read_csv(paste0(infilepath1, 'TR_iwueiav_v3.csv'))

# Prepare data ------------------------------------------------------------

start_period_year <- 1951
end_period_year <- 2010
# Annual mean from all available sites
df_mean <- df %>%
  filter(year >= start_period_year &
           year <= end_period_year) %>%
  drop_na(pre, pet, spei6mon_mean, ai) %>%
  group_by(year) %>%
  summarise(
    pre_mean = mean(pre, na.rm = T),
    pet_mean = mean(pet, na.rm = T),
    spei_mean = mean(spei6mon_mean, na.rm = T),
    ai_mean = mean(ai, na.rm = T), 
    nsite = n()
  )

# Combine precipitation and PET
df_meanpre <- df_mean %>%
  dplyr::select(all_of(c('pre_mean', 'year'))) %>%
  mutate(name = 'Precipitation',
         # Annual
         pre_mean = pre_mean * 12) %>%
  rename(value = pre_mean)

df_meanpet <- df_mean %>%
  dplyr::select(all_of(c('pet_mean', 'year'))) %>%
  mutate(name = 'PET',
         # Annual
         pet_mean = pet_mean * 365) %>%
  rename(value = pet_mean)

df_doublep <- bind_rows(df_meanpre, df_meanpet)

# Figures -----------------------------------------------------------------


# Time series
fig_ts_fun <- function(df, varname) {
  # Time series of dryness 
  
  if (varname %in% c('pre_mean', 'pet_mean')) {
    y_char <- 'Precipitation or PET (mm)'
  } else if (varname == 'spei_mean') {
    y_char <- 'SPEI'
    tag_char <- 'a'
  } else if (varname == 'ai_mean') {
    y_char <- 'AI'
    tag_char <- 'b'
  }
  
  f <- ggplot(data = df_mean,
              aes(x = year,
                  y = get(as.name(varname)))) +
    geom_line() +
    labs(x = 'Year', y = y_char) +
    theme_char
  
  return(f)
}

# fig_ts_fun(df_mean, 'ai_mean')
# Estimate the trends
sensres <- sens.slope(df_mean$spei_mean)
f_ts <- fig_ts_fun(df_mean, 'spei_mean') + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  geom_smooth(method = 'lm', color = 'black') + 
  geom_text(
    x = 1960,
    y = -0.25,
    label = paste('italic(p)', '<~0.001',
                  sep = ''),
    size = 3,
    parse = T
  )
  

if (F) {
  var_name <- 'value'
  ggplot(data = df_doublep,
         aes(
           x = year,
           y = get(as.name(var_name)),
           color = name
         )) +
    geom_line() +
    geom_smooth(method = 'lm') +
    scale_color_manual(values = c('Precipitation' = '#5e3c99',
                                  'PET' = '#e66101')) +
    labs(x = 'Year', y = 'Precipitation or PET (mm)',
         color = NULL) +
    theme_char
}

# Save figures ------------------------------------------------------------

ggsave(
  f_ts, 
  filename = paste0(
    outfilepath1,
    'Figures/',
    'supp_spei_trend.jpeg'
  ),
  units = 'mm',
  dpi = 900,
  width = 100,
  height = 50
)

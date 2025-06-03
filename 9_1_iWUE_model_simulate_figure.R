# Figure 4
# Combine different figures
# Use different data to form different figures
# Then, combine the figures

source('/Users/mw/Documents/Proj_Postdoc/Proj_WUETrait/Revision_202212/iWUE_setting.R')


# Set parameters ----------------------------------------------------------

thetaw_value <- 0.05 * 100  # from m3 m-3 to kg m-2
thetac_value <- 0.8 * 100
D <- seq(0.5, 3.5, 0.1)

# D using (near) realistic values based on site values
# D <- seq(0.5, 0.7, 0.005)

SM <- seq(50, 10, length.out = length(D)) # Inversely related to D

# SM using (near) realistic values based on site values
# SMmin <- 23
# SMmax <- 26
# SM <- seq(SMmax, SMmin, length.out = length(D))


# the Beta function parameter
c_value <- 0.5

Ca <- seq(300, 400, 5)
cainit <- Ca[1]

g1g <- 2.35 # Gymnosperm from Lin et al., 2015/De Kauwe et al., 2015
g1a <- 4.45 # Angiosperm from De kauwe et al., 2015

# Assume a linear change in g1 as ca increases
# Data from Gardner et al., 2022
# eg1: 2; ag1: 1.4 (kPa^0.5)
# eca: 703; aca: 343 (ppm)
eg1 <- 2
ag1 <- 1.4
eca <- 703
aca <- 343
# Slope. Use the same slope
k_g1_ca <- (eg1-ag1) / (eca-aca)
# Interception
i_g1_ca <- eg1 - k_g1_ca*eca


# Settings for plotting
# point_size <- 4.5*5/14
point_size <- 3.5*5/14
axis_text_size <- 5.5


# Functions defined -------------------------------------------------------


beta1_fun <- function(M, thetaw, thetac, c=0.7) {
  
  # Refer to Sabot et al., 2020
  # M - variable indicating moisture (e.g., precipitation or soil moisture)
  # thetaw - wilting point
  # thetac - field capacity
  # c - constant parameter, 1 or 0.5
  
  beta <- ((M-thetaw) / (thetac-thetaw))^c
  return(beta)
  
}

w_sm_fun <- function(D, ca, g1, M, 
                     thetaw = thetaw_value, thetac = thetac_value) {
  # sensitivity of W to SM
  
  beta <- beta1_fun(M, thetaw, thetac, c=c_value)
  dwdsm <- -D^0.5*ca*c_value*g1*beta / (1.6*(D^0.5 + g1*beta)^2*(M - thetaw))
  return(dwdsm)
}




# Factor simulation data prep ---------------------------------------------

# Make a long vector
Dlong <- rep(D, times = length(Ca))
Calong <- rep(Ca, each = length(D))
SMlong <- rep(SM, times = length(Ca))

# Gymnosperms
wsens_sm1 <- w_sm_fun(Dlong, Calong, g1g, SMlong, 
                      thetac = thetac_value, thetaw = thetaw_value)

wsens_sm1_vpdconst <- w_sm_fun(2, Calong, g1g, SMlong, 
                               thetac = thetac_value, thetaw = thetaw_value)

# Angiosperms
wsens_sm1_a <- w_sm_fun(Dlong, Calong, g1a, SMlong, 
                        thetac = thetac_value, thetaw = thetaw_value)

df_wsens_sm <- data.frame(
  VPD = Dlong,
  # VPD = 2, 
  Ca = Calong, 
  SM = SMlong, 
  Sens = wsens_sm1
)


# Factor simulation plots -------------------------------------------------

# Create Heatmaps
# x = SM
barlimt_sim <- c(-5, 0)
palet_sim <- 'OrRd'

barlimt_real <- c(-2,-1) # For near-real data
palet_real <- 'YlOrBr'

f_sens_sm <- ggplot(data = df_wsens_sm) +
  geom_tile(aes(x = SM,
                y = Ca,
                fill = Sens)) +
  scale_fill_distiller(direction = -1,
                       palette = palet_sim,
                       limits = barlimt_sim
                       # limits = c(-6, 0)
                       ) +
                       labs(
                         fill = expression(S[italic(W)]),
                         x = expression(SM ~ '(kg ' * m ^ -2 * ')'),
                         y = expression(italic(C[a]) ~ '(' * mu * mol ~ mol ^ -1 * ')'),
                         tag = 'a'
                       ) +
                         theme_char +
                         theme(legend.key.width = unit(2.5, 'mm'))
f_sens_sm



# Future simulation data prep ---------------------------------------------

infilepath = '/Users/mw/Documents/Proj_Postdoc/Proj_WUETrait/Revision_202212/Results_revise_part1/CMIPOutput/Postprocess/'
infilename1 <- 'model_sim_sens_multiscen.csv'
df_sim <- read_csv(paste0(infilepath, infilename1))


# Future simulation plots -------------------------------------------------

# Create a box plot
colors <- c('gray', '#2c7bb6', '#abd9e9', '#ffffbf', '#fdae61')
mean_historical <- mean(df_sim[df_sim$experiment=='historical', 'sw']$sw, na.rm = T)
f_fsim <- ggplot(data = df_sim, aes(x = experiment, y = sw, fill = experiment)) +
  geom_boxplot(color = "black", 
               #outlier.shape = NA, 
               width = 0.4) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 0.75) + 
  scale_fill_manual(values = colors) + 
  scale_y_continuous(limits = c(-4, -1)) + 
  labs(x = NULL, 
       y = expression(paste(S[italic(W)])), 
       #title = expression("Comparisons of "~S[italic(W)]~"in Different Scenarios"), 
       tag = 'b'
       ) + 
  geom_hline(yintercept = mean_historical, linetype = "dashed", color = "black") +
  theme_char + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "none", 
    plot.title = element_text(size = 6.5)
  )
f_fsim


# Combine figures ---------------------------------------------------------

f_combine1 <- f_sens_sm + f_fsim + plot_layout(widths = c(1, 1.5))

if (F) {
outfigpath <- '/Users/mw/Documents/Proj_Postdoc/Proj_WUETrait/Revision_202212/Results_revise_part1/Figures/'
ggsave(
  f_combine1, 
  filename = paste0(outfigpath, 'fig4_iwue_model_smsens_v3_combine_future_scen.jpeg'), 
  dpi = 600, 
  unit = 'mm', 
  width = 130, 
  height = 65
)
}

# To save as PDF
outfilepath1 = '/Users/mw/Documents/Proj_Postdoc/Proj_WUETrait/Revision_202212/Docs_revise_part1/NatGeo_submit/AIP/Figures_AIP/'
ggsave(
  f_combine1,
  filename = paste0(
    outfilepath1,
    'Figure4.pdf'
  ),
  units = 'mm',
  width = 130,
  height = 65
)


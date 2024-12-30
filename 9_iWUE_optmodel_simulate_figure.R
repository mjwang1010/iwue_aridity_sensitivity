# Figure4, Simulate iWUE variations according to theoretical models
# Simulate sensitivity variations under different Ca, D, and SM

source('Path for iWUE_setting.R')
filepath1 <- 'Path of the input file folder'
outfilepath1 <- infilepath1

# Model functions ---------------------------------------------------------

# Model1 from Medlyn et al., 2011

beta1_fun <- function(M, thetaw, thetac, c=0.7) {
  
  # Refer to Sabot et al., 2020
  # M - variable indicating moisture (e.g., precipitation or soil moisture)
  # thetaw - wilting point
  # thetac - field capacity
  # c - constant parameter, 1 or 0.5
  
  beta <- ((M-thetaw) / (thetac-thetaw))^c
  return(beta)
  
}

m1_fun <- function(D, Ca, g1) {
  # iWUE model
  # D - VPD, kPa
  # Ca - CO2 concentrations, umol mol-1
  # g1 - slope, kPa^0.5
  
  W = D^0.5 * Ca * (D^0.5 + g1)^(-1) / 1.6
  # iWUE
  return(W)
}

m1v2_fun <- function(D, Ca, g1, beta) {
  # iWUE model considering soil moisture impact
  # D - VPD, kPa
  # Ca - CO2 concentrations, umol mol-1
  # g1 - slope, kPa^0.5
  # beta - estimated from a beta function
  
  W = D^0.5 * Ca * (D^0.5 + g1*beta)^(-1) / 1.6
  
  # iWUE
  return(W)
}

g1_model <- function(psil, ca = 320, PFT = 'a') {
  
  # psil - leaf water potential
  # ca - co2 concentration
  # PFT - plant functional types
  
  # Typical values from Manzoni 2011
  ca_star <- 400 # umol mol-1
  gamma_star <- 43 # umol mol-1; Bernacchi 2001
  
  if (PFT == 'a') {
    # deciduous; wet condition
    lambda_max <- 2912 # umol mol-1
    beta <- 0.78
    psil_max <- -1.85 # MPa
    
  } else if (PFT == 'g') {
    # conifer; wet condition
    lambda_max <- 6944 # umol mol-1
    beta <- 0.98
    psil_max <- -1.40 # MPa
    
  }
  
  lambda <- lambda_max * ca / ca_star * exp(-beta * (psil - psil_max)^2) # umol mol-1
  gamma_star_kPa <- gamma_star * 101.3 * 10^(-6) # from umol mol-1 to kPa
  g1 <- (3/1.6 * gamma_star_kPa/(lambda*10^(-6)))^0.5
  
}

# Theoretical sensitivity
w_ca_fun <- function(D, g1, beta) {
  # sensitivity of W to Ca

  dwdca <- D^0.5 / 1.6 / (D^0.5 + g1*beta)
  return(dwdca)
}

w_vpd_fun <- function(ca, g1, beta, D) {
  # sensitivity of W to VPD
  
  dwdvpd <- ca*g1*beta/(2*D^0.5) / (1.6*(D^0.5 + g1*beta)^2)
  return(dwdvpd)
}

w_sm_fun <- function(D, ca, g1, M, thetaw, thetac) {
  # sensitivity of W to SM
  
  beta <- beta1_fun(M, thetaw, thetac)
  dwdsm <- -D^0.5*ca*c*g1*beta / (1.6*(D^0.5 + g1*beta)^2*(M - thetaw))
  return(dwdsm)
}

# Sensitivity model --------------------------------------------

# Case 1: wet condition (high soil moisture)
beta <- 1
g1 <- 2.35
D <- seq(0.5, 2.5, 0.1)
Ca <- seq(300, 400, 5)
# Make a long vector
Dlong <- rep(D, times = length(Ca))
Calong <- rep(Ca, each = length(D))
# Sensitivity to VPD
wsens_vpd <- w_vpd_fun(Calong, g1, beta, Dlong)
df_wsens_vpd <- data.frame(
  VPD = Dlong, 
  Ca = Calong, 
  Sens = wsens_vpd
)

# Case 2: dry condition (low soil moisture)
beta <- 0.5
wsens_vpd2 <- w_vpd_fun(Calong, g1, beta, Dlong)
df_wsens_vpd2 <- data.frame(
  VPD = Dlong, 
  Ca = Calong, 
  Sens = wsens_vpd2
)

# Sensitivity to soil moisture
thetaw <- 0.05 * 100  # from m3 m-3 to kg m-2
thetac <- 0.8 * 100
c <- 0.5
SM <- seq(0.7, 0.1, length.out = length(D))
SM <- SM*100
SMlong <- rep(SM, times = length(Ca))
wsens_sm1 <- w_sm_fun(Dlong, Calong, g1, SMlong, 
                      thetac = thetac, thetaw = thetaw)
df_wsens_sm1 <- data.frame(
  VPD = Dlong, 
  Ca = Calong, 
  SM = SMlong, 
  Sens = wsens_sm1
)

# Case study --------------------------------------------------------------
# A quick test
thetaw <- 0.05 * 100  # from m3 m-3 to kg m-2
thetac <- 0.8 * 100
c <- 0.5
SM <- seq(0.7, 0.1, length.out = length(D))
SM <- SM*100
SMmean <- mean(SM, na.rm = T)
D <- seq(0.5, 2.5, 0.1)
Dmean <- mean(D, na.rm = T)
Ca <- c(312, 388)
g1 <- 2.35

wsens_sm1 <- w_sm_fun(Dmean, Ca, g1, SMmean, 
                      thetac = thetac, thetaw = thetaw)
df_wsens_sm1 <- data.frame(
  VPD = Dmean, 
  Ca = Ca, 
  SM = SMmean, 
  Sens = wsens_sm1
)
wsens_change <- wsens_sm1[2] - wsens_sm1[1]
wsens_change_rel <- wsens_change/wsens_sm1[1]

# Sensitivity figures ------------------------------------------------------

# VPD sensitivity
var_breaks <- seq(10, 40, 5)
f_sens1 <- ggplot(data = df_wsens_vpd) + 
  geom_tile(
    aes(
      x = VPD, 
      y = Ca, 
      fill = Sens
    )
  ) + 
  scale_fill_distiller(
    direction = 1, 
    palette = 'PuBu', 
    breaks = var_breaks
  ) + 
  labs(fill = expression(S[italic(W)]), 
       x = expression(italic(D)~'(kPa)'), 
       y = expression(italic(C[a])~'('*mu*mol~mol^-1*')')) + 
  theme_char + 
  theme(
    legend.key.width= unit(3.5, 'mm')
  )

f_sens2 <- ggplot(data = df_wsens_vpd2) + 
  geom_tile(
    aes(
      x = VPD, 
      y = Ca, 
      fill = Sens
    )
  ) + 
  scale_fill_distiller(
    direction = 1, 
    palette = 'PuBu', 
    breaks = var_breaks
  ) + 
  labs(fill = expression(S[italic(W)]), 
       x = 'VPD (kPa)', 
       y = expression(italic(C[a])~'('*mu*mol~mol^-1*')')) + 
  theme_char + 
  theme(
    legend.key.width= unit(3.5, 'mm')
  )

# SM Sensitivity
f_sens_sm1 <- ggplot(data = df_wsens_sm1) + 
  geom_tile(
    aes(
      x = SM, 
      y = Ca, 
      fill = Sens
    )
  ) + 
  scale_fill_distiller(
    direction = -1, 
    palette = 'OrRd',
  ) + 
  labs(fill = expression(S[italic(W)]), 
       x = expression(SM~'(kg '*m^-2*')'), 
       y = expression(italic(C[a])~'('*mu*mol~mol^-1*')'), 
       tag = 'a') + 
  theme_char + 
  theme(
    legend.key.width= unit(3.5, 'mm')
  )
f_sens_sm1

f_sens_sm2 <- ggplot(data = df_wsens_sm1) + 
  geom_tile(
    aes(
      x = VPD, 
      y = Ca, 
      fill = Sens
    )
  ) + 
  scale_fill_distiller(
    direction = -1, 
    palette = 'BuPu',
  ) + 
  labs(fill = expression(S[italic(W)]), 
       x = expression(italic(D)~'(kPa)'), 
       y = expression(italic(C[a])~'('*mu*mol~mol^-1*')'), 
       tag = 'b') + 
  theme_char + 
  theme(
    legend.key.width= unit(3.5, 'mm')
  )

# Combine figures ---------------------------------------------------------

f_combine1 <- f_sens_sm1 + f_sens_sm2


# Facets ------------------------------------------------

dfout1_gather1$case <- factor(
  dfout1_gather1$case, 
  levels = c('iwue1', 'iwue2', 'iwue3'), 
  labels = c(expression(paste(italic(a),' = 0.5')), 
             expression(paste(italic(a),' = 1')), 
             expression(paste(italic(a),' = 2')))
)

f2 <- ggplot() + 
  geom_line(data = dfout1_gather1,
            aes(
              x = sm,
              y = iwue,
              linetype = clade,
              color = factor(vpd) # or vpd
            )) +
  labs(
    # y = expression('iWUE (' * mu * 'mol' ~ 'mol' ^ {
    #   -1
    # } * ')'),
    y = expression(paste(italic(W) ~ '(' * mu * 'mol' ~ 'mol' ^ {
      -1
    } * ')')),
    x = expression('SM (m' ^ {
      3
    } ~ 'm' ^ {
      -3
    } * ')'),
    color = expression(italic('D')~'(kPa)'),
    linetype = NULL
    ) + 
  scale_linetype_discrete(labels = c('Angiosperms', 'Gymnosperms')) + 
  # 3 levels for VPD
  scale_color_manual(values = c('#018571', '#80cdc1', '#dfc27d')) +
  theme_bw() + 
  theme(
    panel.grid = element_blank(), 
    text = element_text(size = 7), 
    axis.title = element_text(size = 8), 
    axis.text = element_text(size = 7.5), 
    legend.text = element_text(size = 7), 
    legend.title = element_text(size = 7), 
    strip.background = element_rect(fill = 'white'), 
    strip.text = element_text(size = 8)
  ) + 
  facet_wrap(~case, labeller = label_parsed)
f2

if (F) {
# Ca as x variable
f2_v2 <- 
  ggplot() + 
  geom_line(data = dfout1_gather1,
            aes(
              x = ca,
              y = iwue,
              linetype = clade,
              color = factor(sm) # or vpd
            )) +
  labs(y = expression('iWUE ('*mu*'mol'~'mol'^{-1}*')'), 
       color = expression('SM (m'^{3}~'m'^{-3}*')'), 
       # color = expression(italic('D')~'(kPa)'), 
       x = expression(italic('Ca')~'('*mu*'mol'~'mol'^{-1}*')'), 
       linetype = NULL) + 
  scale_linetype_discrete(labels = c('Angiosperms', 'Gymnosperms')) + 
  # 3 levels for VPD
  # scale_color_manual(values = c('#018571', '#80cdc1', '#dfc27d')) +
  theme_bw() + 
  theme(
    panel.grid = element_blank(), 
    text = element_text(size = 7)
  ) + 
  facet_wrap(~case, labeller = label_parsed)
f2_v2
}

# Save figures ------------------------------------------------------------

# Sensitivity figures

# SM sensitivity
ggsave(
  f_combine1, 
  filename = paste0(
    outfilepath1,
    'Figures/',
    'fig4_iwue_model_smsens.jpeg'
  ),
  units = 'mm',
  dpi = 900,
  width = 135,
  height = 60
)

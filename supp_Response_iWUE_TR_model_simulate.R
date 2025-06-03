# New simulations
# Consider g1 acclimation to ca and climate
# For new figures in the main manuscript

source('/Users/mw/Documents/Proj_Postdoc/Proj_WUETrait/Revision_202212/iWUE_setting.R')


# Set input data ---------------------------------------------------------


# Set data 
thetaw_value <- 0.05 * 100  # from m3 m-3 to kg m-2
thetac_value <- 0.8 * 100
D <- seq(0.5, 3.5, 0.1)

# D using (near) realistic values based on site values
# D <- seq(0.5, 0.7, 0.005)

SM <- seq(50, 10, length.out = length(D)) # Inversely related to D
SMmin = 10
SMmax = 50

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
# Theoretical model -------------------------------------------------------

# Beta soil moisture ---------------------------------------------

beta1_fun <- function(M, thetaw, thetac, c=0.7) {
  
  # Refer to Sabot et al., 2020
  # M - variable indicating moisture (e.g., precipitation or soil moisture)
  # thetaw - wilting point
  # thetac - field capacity
  # c - constant parameter, 1 or 0.5
  
  beta <- ((M-thetaw) / (thetac-thetaw))^c
  return(beta)
  
}

# Test different c parameters
if (F) {
  beta_05 <- beta1_fun(SM, thetaw_value, thetac_value, c = 0.5)
  beta_20 <- beta1_fun(SM, thetaw_value, thetac_value, c = 2.0)
  
  plot(x = SM, y = beta_05)
  plot(x = SM, y = beta_20)
}


# Theoretical sensitivity -----------------------------------------------


# Derived from an optimality-theory based model

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

w_sm_fun <- function(D, ca, g1, M, 
                     thetaw = thetaw_value, thetac = thetac_value) {
  # sensitivity of W to SM
  
  beta <- beta1_fun(M, thetaw, thetac, c=c_value)
  dwdsm <- -D^0.5*ca*c_value*g1*beta / (1.6*(D^0.5 + g1*beta)^2*(M - thetaw))
  return(dwdsm)
}


###
# Consider g1 acclimation to ca
# Function 2
# initial conditions: g1g; ca = 300
# Use the same slope: k_g1_ca*(eca - cainit) = g1ac - g1g

w_ca_fun2 <- function(D, g1, beta) {
  # sensitivity of W to Ca
  
  
  return(dwdca)
}

w_vpd_fun2 <- function(ca, g1, beta, D) {
  # sensitivity of W to VPD
  # cainit: initial value of Ca
  # k_g1_ca: the slope of the relationship between g1 and ca
  
  g1ac <- (ca - cainit)*k_g1_ca + g1
  dwdvpd <- ca*g1ac*beta/(2*D^0.5) / (1.6*(D^0.5 + g1ac*beta)^2)
  return(dwdvpd)
}

w_sm_fun2 <- function(D, ca, g1, M, thetaw, thetac) {
  # sensitivity of W to SM
  # 
  
  g1ac <- (ca - cainit)*k_g1_ca + g1 # Acclimation to Ca
  beta <- beta1_fun(M, thetaw, thetac, c = c_value)
  dwdsm <- -D^0.5*ca*c_value*g1ac*beta / (1.6*(D^0.5 + g1ac*beta)^2*(M - thetaw))
  return(dwdsm)
}

# acclimation to climate change
w_sm_fun3 <- function(D, ca, g1ac, M, thetaw, thetac) {
  
  beta <- beta1_fun(M, thetaw, thetac, c = c_value)
  dwdsm <- -D^0.5*ca*c_value*g1ac*beta / (1.6*(D^0.5 + g1ac*beta)^2*(M - thetaw))
  
  return(dwdsm)
}


# Test variation ----------------------------------------------------------

# g1arr <- seq(1.5, 3.1, 0.1) # Simulate g1 variation influenced by climate change
g1arr <- seq(3.37, 5.53, 0.1) # Simulate g1 variation of angiosperms
# g1arr <- (Ca - cainit)*k_g1_ca + g1g

wsens_g1arr <- w_sm_fun(2, 400, g1arr, 40,
                        thetac = thetac_value, thetaw = thetaw_value)


# Separately compute SW at two co2 levels
vpdarrlong = rep(D, times = length(g1arr))
smarrlong = rep(SM, times = length(g1arr))
g1arrlong = rep(g1arr, each = length(D))
wsens_g1_ca300 <- w_sm_fun(vpdarrlong, 300, g1arrlong, smarrlong)
wsens_g1_ca400 <- w_sm_fun(vpdarrlong, 400, g1arrlong, smarrlong)

df_wsens_g1var <- data.frame(
  D = vpdarrlong, 
  SM = smarrlong, 
  g1 = g1arrlong, 
  Sens_ca300 = wsens_g1_ca300, 
  Sens_ca400 = wsens_g1_ca400
)
df_wsens_g1var <- df_wsens_g1var %>% 
  mutate(Sens_diff = Sens_ca400 - Sens_ca300)

# Consider different ca levels
caarr <- c(300, 400)
caarrlong <- rep(caarr, times = length(g1arr))
g1arrlong <- rep(g1arr, each = length(caarr))
wsens_g1arrlong <- w_sm_fun(2, caarrlong, g1arrlong, 40,
                            thetac = thetac_value, thetaw = thetaw_value)

df_wsens_g1_ca <- data.frame(
  ca = caarrlong, 
  g1 = g1arrlong, 
  sens = wsens_g1arrlong
)
# Wide the data frame
df_wsens_g1_ca_wide <- df_wsens_g1_ca %>% 
  pivot_wider(names_from = ca, 
              values_from = sens)

if (F) {
# Visualization test
ggplot(data = df_wsens_g1var) + 
  geom_tile(
    aes(
      x = SM, 
      y = g1, 
      fill = Sens_diff
    )
  ) 

ggplot(data = df_wsens_g1var) + 
  geom_density(aes(x = Sens_diff)) + 
  geom_vline(xintercept = 0, linetype = 'dashed')
}

# Simulation g1 acclimate to ca --------------------------------------------------------------------

# Prepare data for simulations
# Make long data frames
# Make a long vector
Dlong <- rep(D, times = length(Ca))
Calong <- rep(Ca, each = length(D))
SMlong <- rep(SM, times = length(Ca))


# Use the acclimation to ca in the function
wsens_sm1_acc <- w_sm_fun2(Dlong, Calong, g1g, SMlong, 
                           thetac = thetac_value, thetaw = thetaw_value)

# Acclimation
df_wsens_sm1_acc <- data.frame(
  VPD = Dlong, 
  Ca = Calong, 
  SM = SMlong, 
  Sens = wsens_sm1_acc
)


# Simulation g1 No acclimation----------------
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

# Simulation g1 acclimate to climate --------------------------------------

# acclimation to climate change
# Test simulations

# Spatial
g1_aclim_fun <- function(mat, mi, pft) {
  
  # Set parameters
  if (pft == 'g') {
    a4 <- -0.97
  } else if (pft == 'a') {
    a4 <- -0.37
  }
  
  a0 <- 1.32
  a1 <- 0.03 
  a2 <- 0.02
  a3 <- 0.01
  
  g1acc <- exp(a0 + a1*mi + a2*mat + a3*mi*mat + a4)
  return(g1acc)
}



# Time1
Tt1 <- 20
MIt1 <- 2
# Delta of the values
dT <- 1.5
dMI <- -0.2
dTMI <- -Tt1*MIt1+(Tt1+dT)*(MIt1+dMI)

#ratio_gt2_gt1 <- exp(a1*dT + a2*dMI + a3*dTMI)

x <- seq(0.5, 5, 0.5)
df <- data.frame(x= x, 
                 y1 = x, 
                 y2 = x^2, 
                 y3 = x/(x+0.5^0.5)^2)
# Simple tests
ggplot(data = df) + 
  # geom_point(aes(x = x, y = y1)) + 
  # geom_point(aes(x = x, y = y2), 
  #            shape = 21) + 
  geom_point(aes(x = x, y = y3), 
             shape = 24)

###
# Temporal
# Consistent with simulations in Mastrotheodoros et al., 2017, the 1% decline
years = 1951:2010
nyears = length(years)
g1_0 = 2.35 # kPa^0.5 the first year
g1_change_perc = 0.01  # 
D_0 = 0.5
SM_0 = 50

ca_years = seq(300, 400, length.out = nyears)
D_years = seq(0.5, 3.5, length.out = nyears)
SM_years = seq(50, 10, length.out = nyears)

# Estimate g1
g1_seq_fun <- function(g1_year1, change_perc, nyear) {
  g1_seq <- g1_year1*(1-change_perc)^(nyear-1)
  return(g1_seq)
}

g1_years = g1_seq_fun(g1_0, g1_change_perc, 1:nyears)

Dlong <- rep(D_years, times = length(ca_years))
SMlong <- rep(SM_years, times = length(ca_years))
Calong <- rep(ca_years, each = length(D_years))
g1long <- rep(g1_years, each = length(D_years))

wsens_sm1_acc2 <- w_sm_fun3(Dlong, Calong, g1long, SMlong, 
                           thetac = thetac_value, thetaw = thetaw_value)

# Acclimation to climate
df_wsens_sm1_acc2 <- data.frame(
  VPD = Dlong, 
  Ca = Calong, 
  SM = SMlong, 
  g1 = g1long, 
  Sens = wsens_sm1_acc2
)


# Simulation SM equation --------------------------------------------------

# Only simulate the beta curve
c_value = seq(0.2, 0.7, 0.1)
df_beta_sim = data.frame()
for (cone in c_value) {
  beta_arr = beta1_fun(SM, thetaw_value, thetac_value, c = cone)
  df_beta_temp = data.frame(beta = beta_arr, 
                            a = cone, 
                            SM = SM)
  df_beta_sim <- bind_rows(df_beta_sim, df_beta_temp)
}

# Visual the results
# Effect of different values in parameter c(a)
f_beta_sim = ggplot(data = df_beta_sim, 
       aes(x = SM, y = beta)) + 
  geom_point(aes(shape = factor(a))) + 
  labs(x = expression(SM~'(kg '*m^-2*')'), 
       y = expression(beta~'factor'), 
       shape = expression('Parameter'~italic(a))) + 
  theme_char


# Consider different formulas for SM effect

w_sm_beta_fun <- function(D, ca, g1, M, c, 
                     thetaw = thetaw_value, thetac = thetac_value) {
  # sensitivity of W to SM
  beta <- beta1_fun(M, thetaw, thetac, c) # c, input from outside of the function
  dwdsm <- -D^0.5*ca*c_value*g1*beta / (1.6*(D^0.5 + g1*beta)^2*(M - thetaw))
  return(dwdsm)
}

# Gymnosperms
wsens_sm_c <- w_sm_beta_fun(Dlong, Calong, g1g, SMlong, 2, 
                            thetac = thetac_value, thetaw = thetaw_value)

df_wsens_sm_c <- data.frame(
  VPD = Dlong,
  # VPD = 2, 
  Ca = Calong, 
  SM = SMlong, 
  Sens = wsens_sm_c
)


# Figures -----------------------------------------------------------------

# Tile, SM Sensitivity ------------
# x = SM
barlimt_sim <- c(-5, 0)
barlimt_real <- c(-2, -1)

f_sens_sm <- ggplot(data = df_wsens_sm) + 
  geom_tile(
    aes(
      x = SM, 
      y = Ca, 
      fill = Sens
    )
  ) + 
  scale_fill_distiller(
    direction = -1, 
    # palette = 'OrRd',
    palette = 'YlOrBr', 
    limits = barlimt_real
    # limits = c(-6, 0)
  ) + 
  labs(fill = expression(S[italic(W)]), 
       x = expression(SM~'(kg '*m^-2*')'), 
       y = expression(italic(C[a])~'('*mu*mol~mol^-1*')'), 
       tag = 'a') + 
  theme_char + 
  theme(
    legend.key.width= unit(2.5, 'mm')
  )
f_sens_sm

# x = VPD
if (T) {
  f_sens_sm2 <- ggplot(data = df_wsens_sm) + 
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
      limits = c(-5, 0)
    ) + 
    labs(fill = expression(S[italic(W)]), 
         x = expression(italic(D)~'(kPa)'), 
         y = expression(italic(C[a])~'('*mu*mol~mol^-1*')'), 
         tag = 'a') + 
    theme_char + 
    theme(
      legend.key.width= unit(2.5, 'mm')
    )
}


# Tile, Consider the acclimation -------- 

# Acclimation to Ca
f_sens_sm_acc <- ggplot(data = df_wsens_sm1_acc) + 
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
    limits = c(-5, 0)
  ) + 
  labs(fill = expression(S[italic(W)]), 
       x = expression(SM~'(kg '*m^-2*')'), 
       y = expression(italic(C[a])~'('*mu*mol~mol^-1*')'), 
       tag = 'a') + 
  theme_char + 
  theme(
    legend.key.width= unit(2.5, 'mm')
  )


# Acclimation to climate change
f_sens_sm_acc2 <- ggplot(data = df_wsens_sm1_acc2) + 
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
    limits = c(-5, 0)
  ) + 
  labs(fill = expression(S[italic(W)]), 
       x = expression(SM~'(kg '*m^-2*')'), 
       y = expression(italic(C[a])~'('*mu*mol~mol^-1*')'), 
       tag = 'a') + 
  theme_char + 
  theme(
    legend.key.width= unit(2.5, 'mm')
  )

f_sens_sm_acc2

# Tile, different SM equations --------------------------------------------

f_sens_sm_c_v1 <- ggplot(data = df_wsens_sm_c) + 
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
    # limits = c(-6, 0)
  ) + 
  labs(fill = expression(S[italic(W)]), 
       x = expression(SM~'(kg '*m^-2*')'), 
       y = expression(italic(C[a])~'('*mu*mol~mol^-1*')'), 
       tag = 'a') + 
  theme_char + 
  theme(
    legend.key.width= unit(2.5, 'mm')
  )
f_sens_sm_c_v1


###

# Point, change in SW against SM ------------------------------------------

# Relationship between SW change under different ca and SM variation
# Prepare data
sens_diff_fun <- function(df) {
  # Estimate difference of the sensitivity under two levels of ca
  df_wsens_ca1 <- df %>% 
    filter(Ca == 300) %>% 
    rename(Ca1 = Ca, Sens1 = Sens)
  df_wsens_ca2 <- df %>% 
    filter(Ca == 400) %>% 
    rename(Ca2 = Ca, Sens2 = Sens)
  # Compute the difference
  df_wsens_diff <- df_wsens_ca1 %>% 
    left_join(df_wsens_ca2, by = c('VPD'='VPD', 'SM'='SM')) %>% 
    mutate(Sens_diff = Sens2- Sens1)
  return(df_wsens_diff)
}

df_wsens_sm1_diff <- sens_diff_fun(df_wsens_sm)

# Revise the y limits for different data ranges
# In reality, the SM range is small
f_senschange <- ggplot(data = df_wsens_sm1_diff) + 
  geom_point(aes(x = SM, y = Sens_diff), 
             size = point_size, 
             alpha = 0.7) + 
  geom_hline(yintercept = 0, linetype = 'dashed') +
  # geom_hline(yintercept = -0.33, linetype = 'dashed') +
  labs(
    x = expression(SM~'(kg '*m^-2*')'), 
    y = expression(atop(paste('Change in '*S[italic(W)]),
                            paste('Between two levels of '*italic(C[a])))), 
    tag = 'b'
  ) + 
  scale_y_continuous(breaks = seq(-1.2, 0, 0.3)) +
  # scale_y_continuous(breaks = seq(-0.43, -0.33, 0.02)) +
  theme_char + 
  theme(axis.title.y = element_text(size = axis_text_size))
f_senschange


# For VPD variation
f_senschange_vpdvar <- ggplot(data = df_wsens_sm1_diff) + 
  geom_point(aes(x = VPD, y = Sens_diff), 
             size = point_size, 
             alpha = 0.7) + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  labs(
    x = expression(italic(D)~'(kPa)'), 
    y = expression(atop(paste('Change in ' * S[italic(W)]),
                        paste('Between two levels of ' * italic(C[a])))), 
    tag = 'b'
  ) + 
  scale_y_continuous(breaks = seq(-1.2, 0, 0.3)) + 
  theme_char + 
  theme(axis.title.y = element_text(size = axis_text_size))
f_senschange_vpdvar

if (T) {
  # Consider acclimation
  # Acclimation to Ca
  df_wsens_sm1_acc_diff <- sens_diff_fun(df_wsens_sm1_acc)
  
  f_senschange_acc <- ggplot(data = df_wsens_sm1_acc_diff) +
    geom_point(aes(x = SM, y = Sens_diff),
               size = point_size,
               alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    labs(
      x = expression(SM ~ '(kg ' * m ^ -2 * ')'),
      y = expression(atop(paste('Change in ' * S[italic(W)]),
                          paste('Between two levels of ' * italic(C[a])))), 
      tag = 'b'
    ) +
    #scale_y_continuous(breaks = seq(-1.2, 0, 0.3)) +
    theme_char +
    theme(axis.title.y = element_text(size = axis_text_size))
   f_senschange_acc
   
   # Acclimation to climate change
   df_wsens_sm1_acc_diff2 <- sens_diff_fun(df_wsens_sm1_acc2)
   
   f_senschange_acc2 <- ggplot(data = df_wsens_sm1_acc_diff2) +
     geom_point(aes(x = SM, y = Sens_diff),
                size = point_size,
                alpha = 0.7) +
     geom_hline(yintercept = 0, linetype = 'dashed') +
     labs(
       x = expression(SM ~ '(kg ' * m ^ -2 * ')'),
       y = expression(atop(paste('Change in ' * S[italic(W)]),
                           paste('Between two levels of ' * italic(C[a])))), 
       tag = 'b'
     ) +
     theme_char +
     theme(axis.title.y = element_text(size = axis_text_size))
   f_senschange_acc2
}

# Point, change in SW against Ca -----------------------------------------------------------

# Prepare data
sens_diff_fun2 <- function(df, smmin, smmax) {
  # Estimate difference of the sensitivity under two levels of SM
  df_wsens1 <- df %>% 
    filter(SM == smmax) %>% 
    rename(SM1 = SM, Sens1 = Sens)
  df_wsens2 <- df %>% 
    filter(SM == smmin) %>% 
    rename(SM2 = SM, Sens2 = Sens)
  # Compute the difference
  df_wsens_diff <- df_wsens1 %>% 
    left_join(df_wsens2, by = c('Ca'='Ca')) %>% 
    mutate(Sens_diff = Sens2- Sens1)
  return(df_wsens_diff)
}

df_wsens_sm1_diff_cavar <- sens_diff_fun2(df_wsens_sm, SMmin, SMmax)

# Revise the y limits for different data ranges
# 
f_senschange_cavar <- ggplot(data = df_wsens_sm1_diff_cavar) + 
  geom_point(aes(x = Ca, y = Sens_diff), 
             size = point_size, 
             alpha = 0.7) + 
  # geom_hline(yintercept = -3, linetype = 'dashed') + 
  geom_hline(yintercept = -0.22, linetype = 'dashed') + 
  labs(
    x = expression(italic(C[a])~'('*mu*mol~mol^-1*')'), 
    y = expression(atop(paste('Change in '*S[italic(W)]),
                        paste('Between two levels of SM'))), 
    tag = 'c'
  ) + 
  # scale_y_continuous(breaks = seq(-4.2, -3, 0.3)) +
  scale_y_continuous(breaks = seq(-0.32, -0.22, 0.02)) +
  theme_char + 
  theme(axis.title.y = element_text(size = axis_text_size))
f_senschange_cavar

if (T) {
  # Consider acclimation
  # Acclimateion to Ca
  df_wsens_sm1_acc_diff_cavar <- sens_diff_fun2(df_wsens_sm1_acc, SMmin, SMmax)
  
  f_senschange_acc_cavar <- ggplot(data = df_wsens_sm1_acc_diff_cavar) +
    geom_point(aes(x = Ca, y = Sens_diff),
               size = point_size,
               alpha = 0.7) +
    geom_hline(yintercept = -3, linetype = 'dashed') +
    labs(
      x = expression(italic(C[a]) ~ '(' * mu * mol ~ mol ^ -1 * ')'),
      y = expression(atop(paste('Change in ' * S[italic(W)]),
                          paste('Between two levels of SM'))),
      tag = 'c'
    ) +
    scale_y_continuous(breaks = seq(-4.2, -3, 0.3)) +
    theme_char +
    theme(axis.title.y = element_text(size = axis_text_size))
  f_senschange_acc_cavar
  
  # Acclimation to climate change
  df_wsens_sm1_acc_diff_cavar2 <- sens_diff_fun2(df_wsens_sm1_acc2, SMmin, SMmax)
  f_senschange_acc_cavar2 <- ggplot(data = df_wsens_sm1_acc_diff_cavar2) +
    geom_point(aes(x = Ca, y = Sens_diff),
               size = point_size,
               alpha = 0.7) +
    geom_hline(yintercept = -2.5, linetype = 'dashed') +
    labs(
      x = expression(italic(C[a]) ~ '(' * mu * mol ~ mol ^ -1 * ')'),
      y = expression(atop(paste('Change in ' * S[italic(W)]),
                          paste('Between two levels of SM'))),
      tag = 'c'
    ) +
    # scale_y_continuous(breaks = seq(-4.2, -3, 0.3)) +
    theme_char +
    theme(axis.title.y = element_text(size = axis_text_size))
  f_senschange_acc_cavar2
}


# Point, different SM equations -------------------------------------------

df_wsens_sm_c_cadiff <- sens_diff_fun(df_wsens_sm_c)
df_wsens_sm_c_smdiff <- sens_diff_fun2(df_wsens_sm_c)

# change in SW between two levels of ca

ggplot(data = df_wsens_sm_c_cadiff) + 
  geom_point(aes(x = SM, y = Sens_diff), 
             size = point_size, 
             alpha = 0.7) + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  labs(
    x = expression(SM~'(kg '*m^-2*')'), 
    y = expression(atop(paste('Change in '*S[italic(W)]),
                        paste('Between two levels of '*italic(C[a])))), 
    tag = 'b'
  ) + 
  # scale_y_continuous(breaks = seq(-1.2, 0, 0.3)) + 
  theme_char + 
  theme(axis.title.y = element_text(size = axis_text_size))

# change in SW between two levels of SM

ggplot(data = df_wsens_sm_c_smdiff) + 
  geom_point(aes(x = Ca, y = Sens_diff), 
             size = point_size, 
             alpha = 0.7) + 
  # geom_hline(yintercept = -3, linetype = 'dashed') + 
  labs(
    x = expression(italic(C[a])~'('*mu*mol~mol^-1*')'), 
    y = expression(atop(paste('Change in '*S[italic(W)]),
                        paste('Between two levels of SM'))), 
    tag = 'c'
  ) + 
  # scale_y_continuous(breaks = seq(-4.2, -3, 0.3)) +
  theme_char + 
  theme(axis.title.y = element_text(size = axis_text_size))

# Tile, Change in SW under climate acclimation ----------

f_senschange_g1var <- ggplot(data = df_wsens_g1var) + 
  geom_tile(
    aes(
      x = SM, 
      y = g1, 
      fill = Sens_diff
    )
  ) + 
  scale_fill_distiller(
    direction = -1, 
    palette = 'PuBu',
    limits = c(-1.35, 0), 
    breaks = c(-1.2, -0.8, -0.4, 0)
  ) + 
  labs(fill = expression('Change in'~S[italic(W)]), 
       x = expression(SM~'(kg '*m^-2*')'), 
       y = expression(italic(g[1])~'('*k*P*a^0.5*')'), 
       tag = 'a') + 
  theme_char + 
  theme(
    legend.key.width= unit(3.5, 'mm')
  )

# Density plot
color_sensdiff <- '#5e3c99'
f_senschange_dens <- ggplot(data = df_wsens_g1var) + 
  geom_density(
    aes(x = Sens_diff), 
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
  theme_char



# Combine figures -----------

f_sens_combine1 <- f_sens_sm + f_senschange + f_senschange_cavar
f_sens_combine2 <- f_sens_sm2 + f_senschange_vpdvar + f_senschange_cavar
# Acclimation to Ca
f_sens_combine3 <- f_sens_sm_acc + f_senschange_acc + f_senschange_acc_cavar
# Acclimation to climate change
f_sens_combine4 <- f_sens_sm_acc2 + f_senschange_acc2 + f_senschange_acc_cavar2

if (F) {
  f_sens_combine2 <- f_sens_sm_acc + f_senschange_acc
  f_senschange_combine <- f_senschange_g1var + f_senschange_dens
}



# Save data ---------------------------------------------------------------

data_outpath = "out file paths"
write_csv(df_wsens_sm, 
          file = paste0(data_outpath, 'SuppFig5/', 'iwue_sens_against_climate_noacclim_tile.csv'))
write_csv(df_wsens_sm1_diff, 
          file = paste0(data_outpath, 'SuppFig5/', 'iwue_sens_change_against_sm_noacclim_point.csv'))
write_csv(df_wsens_sm1_diff_cavar, 
          file = paste0(data_outpath, 'SuppFig5/', 'iwue_sens_change_against_ca_noacclim_point.csv'))

# Acclimate to ca
write_csv(df_wsens_sm1_acc, 
          file = paste0(data_outpath, 'SuppFig6/', 'iwue_sens_against_climate_acclimca_tile.csv'))
write_csv(df_wsens_sm1_acc_diff, 
          file = paste0(data_outpath, 'SuppFig6/', 'iwue_sens_change_against_sm_acclimca_point.csv'))
write_csv(df_wsens_sm1_acc_diff_cavar, 
          file = paste0(data_outpath, 'SuppFig6/', 'iwue_sens_change_against_ca_acclimca_point.csv'))

# Acclimate to climate
write_csv(df_wsens_sm1_acc2, 
          file = paste0(data_outpath, 'SuppFig7/', 'iwue_sens_against_climate_acclimclim_tile.csv'))
write_csv(df_wsens_sm1_acc_diff2, 
          file = paste0(data_outpath, 'SuppFig7/', 'iwue_sens_change_against_sm_acclimclim_point.csv'))
write_csv(df_wsens_sm1_acc_diff_cavar2, 
          file = paste0(data_outpath, 'SuppFig7/', 'iwue_sens_change_against_ca_acclimclim_point.csv'))

# Beta factor simulations
write_csv(df_beta_sim, 
          file = paste0(data_outpath, 'SuppFig10/', 'iwue_beta_factor_sim.csv'))

# Save figures ------------------------------------------------------------

outfilepath1 <- '/Users/mw/Documents/Proj_Postdoc/Proj_WUETrait/Revision_202212/Results_revise_part1/'

ggsave(
  f_sens_combine1, 
  filename = paste0(
    outfilepath1,
    'Figures/',
    # 'fig4_iwue_model_smsens_v3_smvar_cavar.jpeg'
    'fig_iwue_model_smsens_v3_real_smvar_cavar.jpeg'
  ),
  units = 'mm',
  dpi = 900,
  width = 155,
  height = 50
)

ggsave(
  f_sens_combine2, 
  filename = paste0(
    outfilepath1,
    'Figures/',
    'fig4_iwue_model_smsens_v3_vpdvar_cavar.jpeg'
  ),
  units = 'mm',
  dpi = 900,
  width = 155,
  height = 50
)

# Acclimation 
# To Ca
ggsave(
  f_sens_combine3, 
  filename = paste0(
    outfilepath1,
    'Figures/',
    'fig4_iwue_model_smsens_v3_acc_smvar_cavar.jpeg'
  ),
  units = 'mm',
  dpi = 900,
  width = 155,
  height = 50
)
# To climate change
ggsave(
  f_sens_combine4, 
  filename = paste0(
    outfilepath1,
    'Figures/',
    'fig4_iwue_model_smsens_v3_acc2_smvar_cavar.jpeg'
  ), 
  units = 'mm', 
  dpi = 900, 
  width = 155, 
  height = 50
)

# Simulate different beta values
ggsave(
  f_beta_sim, 
  filename = paste0(
    outfilepath1,
    'Figures/',
    'supp_beta_curve_different_para.jpeg'
  ),
  units = 'mm',
  dpi = 900,
  width = 75,
  height = 75
)

# Original, v1
ggsave(
  f_senschange_combine, 
  filename = paste0(
    outfilepath1,
    'Figures/',
    'supp_iwue_model_smsens_change_g1var.jpeg'
  ),
  units = 'mm',
  dpi = 900,
  width = 135,
  height = 60
)

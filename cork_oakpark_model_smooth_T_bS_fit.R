library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

###################################
# CORK & OAKPARK TEMPERATURE FITS
###################################

# Station coordinates
tower_locations <- data.frame(
  location  = c("Cork","Oakpark"),
  eastings  = c(182624, 272895),
  northings = c(101150, 179267)
)

# Distances to towers
stations <- meteo %>% distinct(east, north)
stations <- stations %>%
  rowwise() %>%
  mutate(
    dist_cork    = sqrt((east - tower_locations$eastings[1])^2 +
                          (north - tower_locations$northings[1])^2),
    dist_oakpark = sqrt((east - tower_locations$eastings[2])^2 +
                          (north - tower_locations$northings[2])^2)
  ) %>%
  ungroup()

# Closest station to each tower
closest_to_cork    <- stations %>% slice_min(dist_cork, n = 1)
closest_to_oakpark <- stations %>% slice_min(dist_oakpark, n = 1)

# Subset meteo for closest stations
meteo_cork    <- meteo %>% semi_join(closest_to_cork, by = c("east","north"))
meteo_oakpark <- meteo %>% semi_join(closest_to_oakpark, by = c("east","north"))

# Daily mean temperature
cork_ts <- meteo_cork %>% mutate(meanT = (TX + TN)/2) %>% select(Dates, meanT) %>% mutate(location = "Cork")
oakpark_ts <- meteo_oakpark %>% mutate(meanT = (TX + TN)/2) %>% select(Dates, meanT) %>% mutate(location = "Oakpark")

# Plot daily temperatures
temps <- bind_rows(cork_ts, oakpark_ts) %>% mutate(Dates = as.Date(Dates))
ggplot(temps, aes(x = Dates, y = meanT)) +
  geom_line(color = "red", size = 1) +
  facet_wrap(~location, ncol = 1, scales = "free_y") +
  labs(title = "Daily Mean Temperature", x = "Date", y = "Temperature (°C)") +
  theme_minimal()

# Fit temperature model: Cork
cork_ts <- meteo_cork %>% mutate(meanT = (TX + TN)/2, t = 1:n()) %>% select(t, meanT)
fit_cork <- nls(meanT ~ a - b * cos(2 * pi / 365 * t - c),
                data = cork_ts,
                start = list(a = mean(cork_ts$meanT),
                             b = (max(cork_ts$meanT) - min(cork_ts$meanT))/2,
                             c = 0))
cork_ts$fitted <- predict(fit_cork)

# Fit temperature model: Oakpark
oakpark_ts <- meteo_oakpark %>% mutate(meanT = (TX + TN)/2, t = 1:n()) %>% select(t, meanT)
fit_oakpark <- nls(meanT ~ a - b * cos(2 * pi / 365 * t - c),
                   data = oakpark_ts,
                   start = list(a = mean(oakpark_ts$meanT),
                                b = (max(oakpark_ts$meanT) - min(oakpark_ts$meanT))/2,
                                c = 0))
oakpark_ts$fitted <- predict(fit_oakpark)

# Extract fitted parameters
a_cork     <- coef(fit_cork)["a"];     b_cork     <- coef(fit_cork)["b"];     delta_cork <- coef(fit_cork)["c"]
a_oakpark  <- coef(fit_oakpark)["a"];  b_oakpark  <- coef(fit_oakpark)["b"];  delta_oakpark <- coef(fit_oakpark)["c"]

##########################################
# SI MODEL WITH SMOOTH TEMPS AND b_S FIT
##########################################

heaviside <- function(x) ifelse(x > 0, 1, 0)

briere_P <- function(T, a_B, T_min, T_max, b) {
  out <- rep(0, length(T))
  valid <- (T >= T_min) & (T <= T_max)
  out[valid] <- a_B * T[valid] * (T[valid] - T_min) / ((T_max - T[valid])^(1/b))
  out
}

b_S_func <- function(T, b_opt=10, T_opt=25, sigma_T=5) {
  b_opt * exp(- (T - T_opt)^2 / (2 * sigma_T^2))
}

sech_kernel <- function(x, s) {
  2 / (s * (exp(2 * x / s) + 2 + exp(-2 * x / s)))
}

# Temperatures from fitted model
T_cork_func <- function(t) a_cork - b_cork * cos(2 * pi * t / 365 - delta_cork)
T_oakpark_func <- function(t) a_oakpark - b_oakpark * cos(2 * pi * t / 365 - delta_oakpark)

# Choose location
selected_location <- "Oakpark" 

# Random parameters
rand_params <- list(
  t_max    = 365,
  P_type   = "briere",
  B        = 6.09,
  K        = 3.31,
  a_B      = 0.59,
  T_min    = 9.31,
  T_max    = 37.21,
  b_briere = 2,
  b_opt    = 5.464,
  T_opt    = 23.929,
  sigma_T  = 5.451,
  K_choice = "sech",
  t_D      = 7.45,
  sech_s   = 0.62,
  tau_m    = 10,
  tau_p    = 5,
  tau_A    = 8,
  delta_P  = 4,
  delta_A  = 4,
  sigma_SP = 0.61,
  sigma_IP = 0.54,
  sigma_SA = 0.78,
  sigma_IA = 0.66,
  K_carry  = 50,
  mu       = 0.9
)

# History initialization
History <- max(rand_params$tau_m, rand_params$tau_p, rand_params$tau_A)
rand_params$init_A_S <- rpois(History+1, lambda=5)
rand_params$init_A_I <- rpois(History+1, lambda=5)
rand_params$init_P_I <- rep(0, History+1)
list2env(rand_params, envir = .GlobalEnv)

start_day <- 1
t_all <- seq(-History, t_max)
n_all <- length(t_all)

# Temperature series
T_series <- if(selected_location == "Cork") {
  T_cork_func(start_day + t_all)
} else {
  T_oakpark_func(start_day + t_all)
}

# Plant abundance (Brière)
hatP_series <- briere_P(T_series, a_B, T_min, T_max, b_briere)
P_tot <- hatP_series

P_I <- numeric(n_all); P_S <- numeric(n_all)
A_S <- numeric(n_all); A_I <- numeric(n_all)

fill_history <- function(x) {
  if(length(x) == 1) rep(x, History+1) else x[1:(History+1)]
}

init_P_I_vec <- pmin(fill_history(init_P_I), hatP_series[1:(History+1)])
P_I[1:(History+1)] <- init_P_I_vec
P_S[1:(History+1)] <- hatP_series[1:(History+1)] - P_I[1:(History+1)]
A_S[1:(History+1)] <- fill_history(init_A_S)
A_I[1:(History+1)] <- fill_history(init_A_I)

# Developmental kernel
svec <- 0:tau_m
K_of_s <- function(s) {
  if(K_choice=="indicator") ifelse(s >= t_D, 1, 0) else sech_kernel(s, sech_s)
}
rawK <- K_of_s(svec)
# Optional: I don't know if we have to normalise K(x)
# if so, then I think we should use:
# if(sum(rawK) > 0) rawK <- rawK / sum(rawK)

# Simulation loop
for(idx in (History+1):(n_all-1)) {
  P_t  <- P_tot[idx]
  PI_t <- P_I[idx]
  AS_t <- A_S[idx]
  AI_t <- A_I[idx]
  
  s_inds <- idx - svec; s_inds[s_inds < 1] <- 1
  b_vals <- b_S_func(T_series[s_inds], b_opt, T_opt, sigma_T)
  B_val <- sum(rawK * b_vals * A_S[s_inds]) * exp(-(AS_t + AI_t)/K_carry)
  
  idx_tauP <- idx - tau_p
  PI_delay <- if(idx_tauP >= 1) P_I[idx_tauP] else 0
  P_delay  <- if(idx_tauP >= 1) P_tot[idx_tauP] else 0
  Phi_A_t <- if(P_delay > 0) exp(- delta_A * (PI_delay / P_delay)) else 1
  
  AS_next <- B_val + Phi_A_t * sigma_SA * AS_t
  AI_next <- (1 - Phi_A_t) * sigma_SA * AS_t + sigma_IA * AI_t
  A_S[idx+1] <- AS_next
  A_I[idx+1] <- AI_next
  
  idx_tauA <- idx - tau_A
  AI_delay_forP <- if(idx_tauA >= 1) A_I[idx_tauA] else 0
  A_delay_forP  <- if(idx_tauA >= 1) (A_S[idx_tauA] + A_I[idx_tauA]) else 0
  Phi_P_t <- if(A_delay_forP > 0) exp(- delta_P * (AI_delay_forP / A_delay_forP)) else 1
  
  PI_next <- (1 - Phi_P_t) * sigma_SP * (P_t - PI_t) + sigma_IP * PI_t
  P_I[idx+1] <- min(max(PI_next, 0), hatP_series[idx+1])
  P_S[idx+1] <- hatP_series[idx+1] - P_I[idx+1]
}

########
# PLOT 
########

res <- data.frame(
  t     = t_all,
  Temp  = T_series,
  hatP  = hatP_series,
  P_tot = P_tot,
  P_I   = P_I,
  P_S   = P_S,
  AS    = A_S,
  AI    = A_I
) %>% filter(t >= 0)

df_plants_long <- res %>% select(t, Susceptible = P_S, Infected = P_I) %>% pivot_longer(-t, names_to="Compartment", values_to="Count")
df_aphids_long <- res %>% select(t, Susceptible = AS, Infected = AI) %>% pivot_longer(-t, names_to="Compartment", values_to="Count")

p_plants <- ggplot(df_plants_long, aes(x=t, y=Count, color=Compartment)) +
  geom_line(linewidth=0.8) +
  scale_color_manual(values=c("darkblue","gold")) +
  labs(x="Day", y="Plant Count") + theme_minimal(base_size=14)

p_aphids <- ggplot(df_aphids_long, aes(x=t, y=Count, color=Compartment)) +
  geom_line(linewidth=0.8) +
  scale_color_manual(values=c("darkblue","gold")) +
  labs(x="Day", y="Aphid Count") + theme_minimal(base_size=14)

(p_plants + p_aphids) +
  plot_layout(ncol = 1, guides = "collect") & 
  theme(legend.position = "right")


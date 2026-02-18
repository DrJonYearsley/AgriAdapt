rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)

setwd("C:/Users/blake/Desktop/R_Code")
load("C:/Users/blake/Desktop/R_Code/meteo_minmaxtemps_2024.RData")

# -------------------------
# aphid DD parameters
# -------------------------
temp_baseline  <- 4
temp_threshold <- 150

model_year <- unique(year(meteo$Dates))[1]

# -------------------------
# Oakpark temperature data
# -------------------------
tower_locations <- data.frame(
  location  = c("Cork","Oakpark"),
  eastings  = c(182624, 272895),
  northings = c(101150, 179267)
)

stations <- meteo %>% distinct(east, north)

stations <- stations %>%
  rowwise() %>%
  mutate(
    dist_oakpark = sqrt((east - tower_locations$eastings[2])^2 +
                          (north - tower_locations$northings[2])^2)
  ) %>% 
  ungroup()

closest_to_oakpark <- stations %>% slice_min(dist_oakpark, n = 1)

meteo_oakpark <- meteo %>%
  semi_join(closest_to_oakpark, by = c("east","north")) %>%
  filter(year(Dates) == model_year) %>%
  arrange(Dates)

oakpark_ts <- meteo_oakpark %>%
  mutate(meanT = (TX + TN)/2) %>%
  select(Dates, meanT)

N <- nrow(oakpark_ts)
Tseries <- oakpark_ts$meanT
day <- 1:N

# -------------------------
# degree day model
# -------------------------
meteo_oakpark$tavg <- (meteo_oakpark$TX + meteo_oakpark$TN)/2
meteo_oakpark$dd <- pmax(meteo_oakpark$tavg - temp_baseline, 0)
cum_gdd <- cumsum(meteo_oakpark$dd)

# -------------------------
# parameters
# -------------------------
b_opt <- runif(1, 1, 15)
T_opt <- runif(1, 15, 30)
sigma_T <- runif(1, 1, 5)

# sigma_SP + sigma_IP < 1 is sufficient (strong assumption) 
# for P_S and P_I are both nonnegative and P_S + P_I = Ptot
sigma_SP <- runif(1, 0.1, 0.1)
sigma_IP <- runif(1, 0.1, 0.9)
sigma_SA <- runif(1, 0.1, 0.45)
sigma_IA <- runif(1, 0.1, 1)

delta_P <- runif(1, 0.1, 1)
delta_A <- runif(1, 0.1, 1)

tau_P <- 7
tau_A <- 14

Ptot <- runif(1,100,1000)
Gamma <- runif(1,10,100)

# -------------------------
# history length
# -------------------------

tau <- max(ceiling(temp_threshold / max(meteo_oakpark$dd)), 
           tau_A, tau_P)

if(tau >= N){
  stop("tau is larger than temperature time series length.")
}

# -------------------------
# emergence times within year
# -------------------------

emerge_map <- vector("list", N)

for(i in 1:N){
  emerge_idx <- which(cum_gdd >= cum_gdd[i] + temp_threshold)[1]
  if(!is.na(emerge_idx) && emerge_idx <= N){
    emerge_map[[emerge_idx]] <- c(emerge_map[[emerge_idx]], i)
  }
}

# -------------------------
# fecundity
# -------------------------

b_S <- b_opt*exp(-(Tseries - T_opt)^2 / (2 * sigma_T^2))

# -------------------------
# initial conditions
# -------------------------

P_I <- numeric(N)
P_S <- numeric(N)
A_I <- numeric(N)
A_S <- numeric(N)
r   <- rep(0,N)
Phi_A <- numeric(N)
Phi_P <- numeric(N)

# History window = 1:tau
P_I[1:tau] <- 0
P_S[1:tau] <- Ptot
A_S[1:tau] <- rpois(tau,5)
#c(rep(0, tau-1), rpois(1,1))
A_I[1:tau] <- rpois(tau,5)
#c(rep(0, tau-1), rpois(1,1))

# -------------------------
# simulation
# -------------------------

start_t <- tau + 1
end_t   <- N - 1

for(t in start_t:end_t){
  
  r_t <- 0
  
  if(!is.null(emerge_map[[t]]) && length(emerge_map[[t]]) > 0){
    r_t <- sum(b_S[emerge_map[[t]]] * A_S[emerge_map[[t]]])
  }
  
  r[t] <- r_t
  
  # plant escape: density-dependent
  if (Tseries[t] < 12) {
    # no transmission (so 100% escape) below 12 degrees
    Phi_P[t] <- 1
  } else {
      Phi_P[t] <- exp(-delta_P*A_I[t - tau_P])
  }
  
  # plant escape: frequency-dependent
  #if (Tseries[t] < 12) {
  #   no transmission below 12°C
  #  Phi_P[t] <- 1
  #} else {
  #  if(A_S[idxP] + A_I[idxP] == 0){
  #    Phi_P[t] <- 1
  #  } else {
  #    Phi_P[t] <- exp(-delta_P*A_I[idxP] / (A_S[idxP] + A_I[idxP]))
  #  }
  #}
  
  # aphid escape: density-dependent
   Phi_A[t] <- exp(-delta_A*P_I[t - tau_A])
  
  # aphid escape: frequency-dependent
  # Phi_A[t] <- exp(-delta_A*P_I[idxA] / Ptot)
  
  # plant dynamics
  P_I[t+1] <- (1 - Phi_P[t])*sigma_SP*P_S[t] + sigma_IP*P_I[t]
  P_S[t+1] <- Ptot - P_I[t+1]
  
  # aphid dynamics
  A_I[t+1] <- (1 - Phi_A[t])*sigma_SA*A_S[t] + sigma_IA*A_I[t]
  A_S[t+1] <- r[t]*exp(-(A_S[t] + A_I[t]) / Gamma) + Phi_A[t]*sigma_SA*A_S[t]
}

df <- tibble(
  day = day,
  Temperature = Tseries,
  r = r,
  P_S = P_S,
  P_I = P_I,
  A_S = A_S,
  A_I = A_I
)

df_plants <- df |> 
  select(day, P_S, P_I) |> 
  pivot_longer(-day, names_to = "state", values_to = "count")

df_aphids <- df |> 
  select(day, A_S, A_I) |> 
  pivot_longer(-day, names_to = "state", values_to = "count")

# -------------------------
# plot
# -------------------------

p_temp <- ggplot(df, aes(day, Temperature)) +
  geom_line(linewidth = 1) +
  theme_minimal()

p_r <- ggplot(df, aes(day, r)) +
  geom_line(linewidth = 1) + 
  theme_minimal()

p_plants <- ggplot(df_plants, aes(day, count, colour = state)) +
  geom_line(linewidth = 1) + 
  theme_minimal()

p_aphids <- ggplot(df_aphids, aes(day, count, colour = state)) +
  geom_line(linewidth = 1) + 
  theme_minimal()

p_temp; p_r; p_plants; p_aphids

# -------------------------
# total aphids plot
# -------------------------

#df$total_aphids <- df$A_S + df$A_I

#p_total_aphids <- ggplot(df, aes(day, total_aphids)) +
#  geom_line(linewidth = 1) +
#  theme_minimal() +
#  labs(y = "Total Aphids")

#p_total_aphids


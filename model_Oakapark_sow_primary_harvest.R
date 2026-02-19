# proposed model for spring barley and S. avenae
# with start date 1st Jan, end date 31st Dec

# growing season:
# sow day -> harvest day
# plant dynamics only within this interval

# primary invasion/infection day: 
# where A_I > 0 for first time
# primary invasion/infection count could be the first observed alate 
# aphid in suction tower dataset of given year
# aphids can grow and decline past harvest day?


rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(patchwork)

setwd("C:/Users/blake/Desktop/R_Code")
load("C:/Users/blake/Desktop/R_Code/meteo_minmaxtemps_2024.RData")

# 2024 temp data selected for now to see how things behave

# aphid DD parameters gotten from 
# https://ipm.ucanr.edu/PHENOLOGY/ma-english_grain_aphid.html
temp_baseline  <- 4
temp_threshold <- 150

model_year <- unique(year(meteo$Dates))[1]

# Oakpark temperature data
tower_locations <- data.frame(
  location  = c("Cork","Oakpark"),
  eastings  = c(182624, 272895),
  northings = c(101150, 179267)
)

# Euclidean distance used
stations <- meteo %>% distinct(east, north) %>%
  rowwise() %>%
  mutate(dist_oakpark = sqrt((east - tower_locations$eastings[2])^2 +
                               (north - tower_locations$northings[2])^2)) %>% ungroup()

closest_to_oakpark <- stations %>% slice_min(dist_oakpark, n = 1)

meteo_oakpark <- meteo %>%
  semi_join(closest_to_oakpark, by = c("east","north")) %>%
  filter(year(Dates) == model_year) %>%
  arrange(Dates)

# mean daily temp data
oakpark_ts <- meteo_oakpark %>%
  mutate(meanT = (TX + TN)/2) %>%
  select(Dates, meanT)

N <- nrow(oakpark_ts)
Tseries <- oakpark_ts$meanT
day <- 1:N

# degree day model
meteo_oakpark$tavg <- (meteo_oakpark$TX + meteo_oakpark$TN)/2
meteo_oakpark$dd <- pmax(meteo_oakpark$tavg - temp_baseline, 0)
cum_gdd <- cumsum(meteo_oakpark$dd)

# emergence times, given that start times known
# is this good enough for telescoping,
# or is this valid only for egg laying (overwintering/diapause)?

emerge_map <- vector("list", N)
for(i in 1:N){
  emerge_idx <- which(cum_gdd >= cum_gdd[i] + temp_threshold)[1]
  if(!is.na(emerge_idx)){
    emerge_map[[emerge_idx]] <- c(emerge_map[[emerge_idx]], i)
  }
}

# parameters
b_opt   <- runif(1, 1, 8)
T_opt   <- runif(1, 20, 35)
sigma_T <- runif(1, 1, 10)

sigma_SP <- runif(1, 0.1, 0.5)
sigma_IP <- runif(1, 0.1, 0.5)
sigma_SA <- runif(1, 0.1, 1)
sigma_IA <- runif(1, 0.1, 1)

delta_P <- runif(1, 1, 30)
delta_A <- runif(1, 1, 30)

tau_P <- 7
tau_A <- 14

Ptot  <- runif(1,100,1000)
Gamma <- runif(1,10,100)

# fecundity
b_S <- b_opt * exp(-(Tseries - T_opt)^2 / (2 * sigma_T^2))

# state variables
P_I  <- numeric(N)
P_S  <- numeric(N)
A_I  <- numeric(N)
A_S  <- numeric(N)
r    <- numeric(N)
Phi_A <- numeric(N)
Phi_P <- numeric(N)

# growing season
sow_day     <- 32    # start of Feb
harvest_day <- 243   # end of Aug

# primary infection day between 10 and 40 days of sow day
primary_day <- round(runif(1, sow_day + 10, sow_day + 40))

# simulation
for(t in 1:(N-1)){
  
  # days before sowing we know everything is ~0
  if(t < sow_day){
    P_S[t+1] <- 0
    P_I[t+1] <- 0
    A_S[t+1] <- 0
    A_I[t+1] <- 0
    r[t]     <- 0
    Phi_P[t] <- 1
    Phi_A[t] <- 1
    next
  }
  
  # initialise crops on sow day (no infectious plants, only susceptible)
  if(t == sow_day){
    P_S[t] <- Ptot
    P_I[t] <- 0
  }
  
  # plant dynamics within growing season
  if(t >= sow_day && t < harvest_day){
    idxP <- ifelse(t - tau_P > 0, t - tau_P, NA)
    if(is.na(idxP) || (A_S[idxP] + A_I[idxP]) == 0){
      Phi_P[t] <- 1
    } else {
      Phi_P[t] <- exp(-delta_P * A_I[idxP] / (A_S[idxP] + A_I[idxP]))
    }
    
    P_I[t+1] <- (1 - Phi_P[t]) * sigma_SP * P_S[t] + sigma_IP * P_I[t]
    P_S[t+1] <- Ptot - P_I[t+1]
  }
  
  # primary infection
  if(t == primary_day){
    A_obs <- round(runif(1, 2, 15))
    A_I[t] <- sample(1:(A_obs-1), 1)   # at least 1 infectious
    A_S[t] <- A_obs - A_I[t]           # at least 1 susceptible
  }
  
  # aphid dynamics
  if(t >= primary_day){
    
    # r is sum of delayed fecundity x abundance 
    # for aphids that have emerged as adults on day t
    
    r_t <- 0
    if(!is.null(emerge_map[[t]]) && length(emerge_map[[t]]) > 0){
      r_t <- sum(b_S[emerge_map[[t]]] * A_S[emerge_map[[t]]])
    }
    r[t] <- r_t
    
    idxA <- ifelse(t - tau_A > 0, t - tau_A, NA)
    if(is.na(idxA)){
      Phi_A[t] <- 1
    } else {
      Phi_A[t] <- exp(-delta_A * P_I[idxA] / Ptot)
    }
    
    A_I[t+1] <- (1 - Phi_A[t]) * sigma_SA * A_S[t] + sigma_IA * A_I[t]
    A_S[t+1] <- r[t] * exp(-(A_S[t] + A_I[t]) / Gamma) + Phi_A[t] * sigma_SA * A_S[t]
  }
  
  # remove plants on harvest day
  if(t >= harvest_day){
    P_S[t+1] <- 0
    P_I[t+1] <- 0
  }
}

# plot
df <- tibble(
  day = day,
  Temperature = Tseries,
  r = r,
  P_S = P_S,
  P_I = P_I,
  A_S = A_S,
  A_I = A_I
)

df$total_aphids <- df$A_S + df$A_I

df_plants <- df |> 
  select(day, P_S, P_I) |> 
  pivot_longer(-day, names_to = "state", values_to = "count")

df_aphids <- df |> 
  select(day, A_S, A_I) |> 
  pivot_longer(-day, names_to = "state", values_to = "count")

vlines <- data.frame(
  event = c("Sowing", "Primary infection", "Harvest"),
  day   = c(sow_day, primary_day, harvest_day)
)

p_plants <- ggplot(df_plants, aes(day, count, colour = state)) +
  geom_line(linewidth = 1) +
  geom_vline(data = vlines, aes(xintercept = day), linetype = "dashed", colour = "black") +
  theme_minimal() +
  labs(title = "SI Plant Dynamics")

p_aphids <- ggplot(df_aphids, aes(day, count, colour = state)) +
  geom_line(linewidth = 1) +
  geom_vline(data = vlines, aes(xintercept = day), linetype = "dashed", colour = "black") +
  theme_minimal() +
  labs(title = "SI Aphid Dynamics")

p_total_aphids <- ggplot(df, aes(day, total_aphids)) +
  geom_line(linewidth = 1) +
  geom_vline(data = vlines, aes(xintercept = day), linetype = "dashed", colour = "black") +
  theme_minimal() +
  labs(y = "count", title = "Total Aphid Dynamics")

p_r <- ggplot(df, aes(day, r)) +
  geom_line(linewidth = 1) +
  geom_vline(data = vlines, aes(xintercept = day), linetype = "dashed", colour = "black") +
  theme_minimal() +
  labs(title = "r")

p_aphids / p_plants

p_infected_plants <- ggplot(df, aes(day, P_I)) +
  geom_line(linewidth = 1, colour = "darkblue") +
  geom_vline(data = vlines, aes(xintercept = day), linetype = "dashed", colour = "black") +
  theme_minimal() +
  labs(title = "Infected Plant Dynamics", y = "Count", x = "Day")

p_susceptible_aphids <- ggplot(df, aes(day, A_S)) +
  geom_line(linewidth = 1, colour = "gold") +
  geom_vline(data = vlines, aes(xintercept = day), linetype = "dashed", colour = "black") +
  theme_minimal() +
  labs(title = "Susceptible Aphid Dynamics", y = "Count", x = "Day")

p_infected_aphids <- ggplot(df, aes(day, A_I)) +
  geom_line(linewidth = 1, colour = "darkblue") +
  geom_vline(data = vlines, aes(xintercept = day), linetype = "dashed", colour = "black") +
  theme_minimal() +
  labs(title = "Infected Aphid Dynamics", y = "Count", x = "Day")

p_infected_plants / p_susceptible_aphids / p_infected_aphids





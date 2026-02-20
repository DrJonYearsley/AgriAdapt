# =============================================================
# A discrete-time temperature-dependent host-vector model of 
# persistent virus dynamics; with application to BYDV 
# transmission between spring barley and S. avenae.
# =============================================================


############
############
#### Q1 ####
############
############

# growing season approach?
# sow day -> harvest day
# plant dynamics only within this interval
# 0 before and after

############
############
#### Q2 ####
############
############

# primary invasion/infection day: 
# where A_I and/or A_S > 0 for first time (both non-zero for now)
# primary invasion/infection abundance could be the first observed 
# alate aphid count in the suction tower dataset for given year

############
############
#### Q3 ####
############
############

# aphids dynamics assumed to change past harvest day (A_S, A_I
# can be >0 past harvest_day)
# is this a fair assumption??

rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(patchwork)

year_select <- sample(2010:2024, 1)

# reset to your own work directory where temp files are
setwd("C:/Users/blake/Desktop/R_Code")

# load temp file for selected year
load(paste0("C:/Users/blake/Desktop/R_Code/meteo_minmaxtemps_", 
            year_select, ".RData"))
message("year chosen: ", year_select)

# Oakpark temperature data
tower_locations <- data.frame(
  location  = c("Cork","Oakpark"),
  eastings  = c(182624, 272895),
  northings = c(101150, 179267)
)

# Euclidean distance used to select closest station
stations <- meteo %>% 
  distinct(east, north) %>%
  rowwise() %>%
  mutate(dist_oakpark = sqrt((east - tower_locations$eastings[2])^2 +
                               (north - tower_locations$northings[2])^2)) %>% 
  ungroup()

closest_to_oakpark <- stations %>% slice_min(dist_oakpark, n = 1)

meteo_oakpark <- meteo %>%
  semi_join(closest_to_oakpark, by = c("east","north")) %>%
  filter(year(Dates) == year_select) %>%
  arrange(Dates)

# mean daily temperature time series
oakpark_ts <- meteo_oakpark %>%
  mutate(meanT = (TX + TN)/2) %>%
  select(Dates, meanT)

N <- nrow(oakpark_ts)
Tseries <- oakpark_ts$meanT
day <- 1:N

############
############
#### Q4 ####
############
############

# is this DD model okay for telescoping??

# aphid DD parameters from 
# https://ipm.ucanr.edu/PHENOLOGY/ma-english_grain_aphid.html
temp_baseline  <- 4
temp_threshold <- 150

# degree day model
meteo_oakpark$tavg <- (meteo_oakpark$TX + meteo_oakpark$TN)/2
meteo_oakpark$dd   <- pmax(meteo_oakpark$tavg - temp_baseline, 0)
cum_gdd <- cumsum(meteo_oakpark$dd)

# aphid adult emergence
emerge_map <- vector("list", N)
max_delay  <- 0

for(i in 1:N){
  emerge_idx <- which(cum_gdd >= cum_gdd[i] + temp_threshold)[1]
  if(!is.na(emerge_idx)){
    emerge_map[[emerge_idx]] <- c(emerge_map[[emerge_idx]], i)
    # get max start -> emergence time
    max_delay <- max(max_delay, emerge_idx - i)
  }
}

# parameters
b_opt   <- 15.3
T_opt   <- 20
sigma_T <- 3.3

# sigma_SP + sigma_IP \leq 1 for well-defined P dynamics
sigma_SP <- 0.62
sigma_IP <- 0.31
sigma_SA <- 0.56
sigma_IA <- 0.53

delta_P <- 0.9
delta_A <- 0.87

tau_P <- 7
tau_A <- 9

Ptot  <- 1000
Gamma <- 100

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
sow_day     <- 32    # sow at start of Feb
harvest_day <- 240   # harvest end of Aug

############
############
#### Q5 ####
############
############

# is this a good way of doing this?

# primary infection/invasion day between 10 and 40 days after sowing
primary_day <- round(runif(1, sow_day + 10, sow_day + 40))

# get max time delay, so can initialise model
tau <- max(max_delay, tau_P, tau_A)

# we need to know tau initial conditions a priori

# plant initial conditions
if(sow_day > 1){
  # zero before sow_day
  P_S[1:(sow_day-1)] <- 0
  P_I[1:(sow_day-1)] <- 0
}

# Ptot susceptible plants only at sow_day
P_S[sow_day] <- Ptot
P_I[sow_day] <- 0

# keep plants fully susceptible up to primary_day
# as no infectious aphids haveentered system
if(primary_day > sow_day){
  P_S[(sow_day + 1):(primary_day-1)] <- Ptot
  P_I[(sow_day + 1):(primary_day-1)] <- 0
}

# if tau > primary_day we need more known initial conditions
extra <- max(0, tau - primary_day)
if(extra > 0){
  start_day <- primary_day
  end_day   <- min(N, primary_day + extra)
  if(start_day <= end_day){
    seed_days <- start_day:end_day
    P_S[seed_days] <- Ptot  #keep as Ptot???
    P_I[seed_days] <- 0
  }
}

# aphid initial conditions
# zeros before primary_day
A_S[1:(primary_day-1)] <- 0
A_I[1:(primary_day-1)] <- 0

# small immigration into field on primary_day
A_obs <- round(runif(1, 2, 5))
A_I[primary_day] <- sample(1:(A_obs-1), 1)
A_S[primary_day] <- A_obs - A_I[primary_day]

# need extra initial conditions if tau > primary_day
if(extra > 0){
  start_day <- primary_day + 1
  end_day   <- min(N, primary_day + extra)
  if(start_day <= end_day){
    seed_days <- start_day:end_day
    # small random counts???
    A_S[seed_days] <- round(runif(length(seed_days), 1,5))
    A_I[seed_days] <- round(runif(length(seed_days), 1,5))
  }
}

# simulation
for(t in tau:(N-1)){
  
  # plant dynamics
  if(t < harvest_day){
    
      denom <- A_S[t - tau_P] + A_I[t - tau_P]
      # if denom==0 then infection pressure is zero => Phi_P = 1
      if(denom > 0){
        Phi_P[t] <- exp(-delta_P * A_I[t - tau_P] / denom)
      } else {
        Phi_P[t] <- 1
      }
    
    P_I[t+1] <- (1 - Phi_P[t]) * sigma_SP * P_S[t] + sigma_IP * P_I[t]
    P_S[t+1] <- Ptot - P_I[t+1]
  }
  
  # aphid dynamics
  if(t > primary_day){
    
    # r is sum of delayed fecundity x abundance 
    # for aphids that have emerged as adults on day t
    if(length(emerge_map[[t]]) > 0){
      # only sum over emergence indices that are >=1 and <=N, i.e. within year
      idxs <- emerge_map[[t]]
      idxs <- idxs[idxs >= 1 & idxs <= N]
      if(length(idxs) > 0){
        r[t] <- sum(b_S[idxs] * A_S[idxs])
      } else {
        r[t] <- 0
      }
    } else {
      r[t] <- 0
    }
    
    Phi_A[t] <- exp(-delta_A * P_I[t - tau_A] / Ptot)
    
    A_I[t+1] <- (1 - Phi_A[t]) * sigma_SA * A_S[t] + sigma_IA * A_I[t]
    A_S[t+1] <- r[t] * exp(-(A_S[t] + A_I[t]) / Gamma) + Phi_A[t] * sigma_SA * A_S[t]
  }
  
  # harvest
  if(t >= harvest_day){
    P_S[t+1] <- 0
    P_I[t+1] <- 0
    # do we choose to make aphids go to 0 here???
    # A_S[t+1] <- 0 
    # A_I[t+1] <- 0
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
  geom_vline(data = vlines, aes(xintercept = day),
             linetype = "dashed", colour = "black") +
  theme_minimal() +
  labs(title = "SI Plant Dynamics")

p_aphids <- ggplot(df_aphids, aes(day, count, colour = state)) +
  geom_line(linewidth = 1) +
  geom_vline(data = vlines, aes(xintercept = day),
             linetype = "dashed", colour = "black") +
  theme_minimal() +
  labs(title = "SI Aphid Dynamics")

p_total_aphids <- ggplot(df, aes(day, total_aphids)) +
  geom_line(linewidth = 1) +
  geom_vline(data = vlines, aes(xintercept = day),
             linetype = "dashed", colour = "black") +
  theme_minimal() +
  labs(y = "count", title = "Total Aphid Dynamics")

(p_aphids / p_plants) | p_total_aphids
message("year chosen: ", year_select)

p_infected_aphids <- ggplot(df, aes(day, A_I)) +
  geom_line(colour = "darkblue", linewidth = 1) +
  geom_vline(data = vlines, aes(xintercept = day),
             linetype = "dashed", colour = "black") +
  theme_minimal() +
  labs(y = "count", title = "Infected Aphids")

p_infected_plants <- ggplot(df, aes(day, P_I)) +
  geom_line(colour = "darkblue", linewidth = 1) +
  geom_vline(data = vlines, aes(xintercept = day),
             linetype = "dashed", colour = "black") +
  theme_minimal() +
  labs(y = "count", title = "Infected Plants")

p_total_aphids_sep <- ggplot(df, aes(day, total_aphids)) +
  geom_line(colour = "black", linewidth = 1) +
  geom_vline(data = vlines, aes(xintercept = day),
             linetype = "dashed", colour = "black") +
  theme_minimal() +
  labs(y = "count", title = "Total Aphids")

p_susceptible_aphids <- ggplot(df, aes(day, A_S)) +
  geom_line(colour = "gold", linewidth = 1) +
  geom_vline(data = vlines, aes(xintercept = day),
             linetype = "dashed", colour = "black") +
  theme_minimal() +
  labs(y = "count", title = "Susceptible Aphids")

(print(p_infected_aphids) + print(p_infected_plants))/ 
  (print(p_total_aphids_sep) + print(p_susceptible_aphids))

# double check that Ptot is conserved during growing season
all(P_I[sow_day:harvest_day] + P_S[sow_day:harvest_day] == Ptot)

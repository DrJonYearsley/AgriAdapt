# ==============================================================================
# Temperature phenology model
# ==============================================================================
#
# Project: Agriadapt
#
# Description: 
# Use daily weather data to calculate degree day model
# 
# Name of author: Jon Yearsley & Ultan O'Donnell
# Contact email: jon.yearsley@ucd.ie
# Date created: 10/02/2026
# ==============================================================================

# Clear Workspace
rm(list = ls())

# Load required libraries
library(data.table) # For importing data
library(dplyr) # For Pipe Operations
library(lubridate) # For floor_tabe
library(ggplot2) # For Plotting
library(Matrix)

# ============================================
# 1. SETUP & PARAMETERS (Define these first!)  -------
# ============================================

setwd("~/git_repos/AgriAdapt/")

# Parameters
start_year <- 2010
nYears = 10   # Number of years for model
temp_baseline = 5    # Degree day baseline (deg C)
temp_threshold = 170  # Degree day threshold  (day deg C)
meteofile =  "Data/dly575.csv"   # Data set for Oakpark

# +++++++++++++++++ No need to edit below this line ++++++++++++++++++++++++++++


year_list = start_year + c(0:(nYears-1))







# ============================================
# 2. Load met eireann data -----
# ============================================

meteo <- fread(meteofile)

# Wrangle
meteo$date = as.IDate(meteo$date, format="%d-%b-%Y")
meteo$year = year(meteo$date)
meteo$doy = yday(meteo$date)

# Subset data
meteo_sub = subset(meteo, year %in% (start_year-1 + c(1:nYears)))

# Inspect structure and plot map
str(meteo_sub)
meteo_agg = aggregate(cbind(maxtp,mintp)~doy, 
                      FUN=mean, 
                      data=meteo_sub, 
                      na.rm=TRUE)

ggplot(data=meteo_agg,
       aes(x=doy,
           y=maxtp,
          colour=mintp)) + 
  geom_point() + 
  scale_colour_fermenter(palette="RdBu") + 
  theme_bw()



# ============================================
# 3. Preprocess: daily avg temp, degree days-------
# ============================================


# Calculate degree days and order rows  date
meteo_dd = meteo_sub %>%
  select(date,maxtp,mintp) %>%
  arrange(date) %>% 
  mutate(
    tavg = (maxtp + mintp) / 2,
    dd = pmax(tavg - temp_baseline, 0),
    day_index = floor(as.numeric(difftime(date,min(date, na.rm=TRUE),units='days')))
  )


# ============================================
# 4. Calculate emergence dates   -------
# ============================================

# Accumulated degree days
cum_gdd_sub <- c(0,cumsum(meteo_dd$dd)) 

results <- data.frame(start=meteo_dd$date,
                      emerge=as.IDate(NA)) 

for (i in 1:nrow(meteo_dd)) { 
  # Find day when development is complete
  emerge_idx  <- which(cum_gdd_sub >= cum_gdd_sub[i] + temp_threshold)[1] 
  
  if (!is.na(emerge_idx)) { 
    results$emerge[i] <- meteo_dd$date[emerge_idx]
  }
}


# Plot results
ggplot(data=results,
       aes(x=yday(start),
           y=yday(emerge),
           colour=factor(year(start)))
       ) + 
  geom_point() +
  coord_equal() +
  scale_colour_brewer("Year", palette="Dark2") +
  theme_bw()



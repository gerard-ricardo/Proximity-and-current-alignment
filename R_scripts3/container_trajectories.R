library(dplyr)
# 1 Import data -----------------------------------------------------------
load("./Rdata/end_waypoints.RData")
start_df = data2
load("./Rdata/start_waypoints.RData")
end_df = data1
# wrangling ---------------------------------------------------------------
str(start_df)
str(end_df)
start_df$Colony  
end_df = end_df %>% rename(Colony = SampleNumber)
missing_colonies <- setdiff(end_df$Colony, start_df$Colony)
missing_colonies
end_df <- end_df %>% 
  mutate(Colony = case_when(
    Colony == "C.4" ~ "C.5",
    Colony == "C.2" ~ "C.9",
    Colony == "C.1" ~ "C.6",
    TRUE ~ Colony
  ))
joined_df = left_join(start_df, end_df, by = 'Colony')
joined_df <- joined_df[complete.cases(joined_df), ]
nrow(joined_df)
# -------------------------------------------------------------------------
library(ggmap)
library(ggplot2)
library(dplyr)
df <- joined_df %>% 
  rename(StartLat = N.x, StartLon = E.x, EndLat = N.y, EndLon = E.y)
(centre_lat <- mean(c(df$StartLat, df$EndLat)) )
(centre_lon <- mean(c(df$StartLon, df$EndLon)))
basemap <- get_googlemap(center = c(lon = centre_lon, lat = centre_lat), zoom = 16, maptype = "satellite")
ggmap(basemap) +
  geom_segment(data = df,
               aes(x = StartLon, y = StartLat, xend = EndLon, yend = EndLat),
               arrow = arrow(length = unit(0.2,"cm")), colour = "blue", size = 0.7) +
  geom_point(data = df, aes(x = StartLon, y = StartLat), colour = "red", size = 2) +
  geom_point(data = df, aes(x = EndLon, y = EndLat), colour = "darkgreen", size = 2) +
  theme_minimal()
# analysis -------------------------------------------------------------------
library(geosphere)
df$Bearing <- bearingRhumb(cbind(df$StartLon, df$StartLat), cbind(df$EndLon, df$EndLat)) %% 360
df$Bearing2 = df$Bearing +180
bear_rad <- df$Bearing * pi / 180
C <- mean(cos(bear_rad))
S <- mean(sin(bear_rad))
mean_rad <- atan2(S, C) %% (2*pi)
R <- sqrt(C^2 + S^2)
n <- length(bear_rad)
se <- sqrt((2 * (1 - R)) / (n * R))
ci_lower <- (mean_rad - 1.96 * se) %% (2*pi)
ci_upper <- (mean_rad + 1.96 * se) %% (2*pi)
mean_deg <- mean_rad * 180 / pi
ci_deg <- c(ci_lower, ci_upper) * 180 / pi
mean_deg +180
ci_deg +180
circular_sd <- sqrt(-2 * log(R)) * 180 / pi
circular_sd
# advection ---------------------------------------------------------------
df$Distance_m <- distHaversine(cbind(df$StartLon, df$StartLat), cbind(df$EndLon, df$EndLat))
df$DeltaTime_hr <- as.numeric(difftime(df$EndTime, df$StartTime, units = "hours"))
df$Velocity_ms <- df$Distance_m / (df$DeltaTime_hr * 3600)

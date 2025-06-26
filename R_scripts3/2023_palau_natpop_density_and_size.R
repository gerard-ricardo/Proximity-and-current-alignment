# load libraries ----------------------------------------------------------
library(dplyr)
library(ggplot2)
library(truncnorm)
library(tidybayes)
library(spatstat)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2")
# import data -------------------------------------------------------------
load("./Rdata/2023palau.natural.den.RData")
# labeling and wrangling -----------------------------------------------
str(data1)
data1$side <- as.factor(as.character(data1$side))
data2 <- data1 %>% split(., data1$side)
data2$left$neg <- -data2$left$dist_m_x
data1 <- data.frame(y = data1$dist_m_y, x = c(data2$right$dist_m_x, data2$left$neg))
# Data exploration ------------------------------------------------------
p0 <- ggplot() +
  geom_point(data1, mapping = aes(x = x, y = y), position = position_jitter(width = .02, height = .02), alpha = 0.50, size = 3)
p0
# spatial clustering ------------------------------------------------------
rslt2 <- as.ppp(data1, W = owin(c(0, 50), c(-2.5, 2.5)))
plot(as.solist(rslt2))
df <- data.frame(summary1 = c(summary(rslt2$`1`)[2], summary(rslt2$`2`)[2], summary(rslt2$`3`)[2]))
summary(rslt2)[2]
donnelly <- unname(clarkevans(rslt2)[2])
donnelly
clarkevans(rslt2)
plot(Kest(rslt2, correction = c("best")))
fitted_models <- list(
  Poisson = ppm(rslt2, ~1),
  Thomas = kppm(rslt2 ~1, clusters = "Thomas", method = "clik2"),
  Matern = kppm(rslt2 ~1, clusters = "MatClust", method = "clik2"),
  HardCore = kppm(rslt2 ~1, clusters = "Cauchy", method = "clik2")
)
(model_comparison <- sapply(fitted_models, AIC))
fitted_thomas <- kppm(rslt2 ~1, clusters = "Thomas", method = "clik2")
fitted_matern <- kppm(rslt2 ~1, clusters = "MatClust", method = "clik2")
scale_thomas <- sqrt(fitted_thomas$par[2])
dist <- nndist(rslt2)
quantile(dist)
median(dist)
n_n.med <- lapply(dist, median)
n_n.mean_all <- mean(unlist(n_n.med))
n_n.sd_all <- sd(unlist(n_n.med))
print(paste0('spacing slope: ', round(median(dist),2), '+/-', round(n_n.sd_all,2) ))
sd(dist)
str(data1)
transect_width = 5
transect_length = 50
den1 = nrow(data1) / (transect_width * transect_length)
print(paste0('density slope: ',den1))
sd_count = sd(data1$y)
sd_density = sd_count / (transect_width * transect_length)
print(paste0('SD of density slope: ', sd_density))
plot(density(unlist(dist)))
par(mfrow = c(3, 1), mar = c(2, 2, 2, 2) + 0.1)
plot(envelope(rslt2, fun = Kest, nsim = 780, nrank = 20))
plot(envelope(rslt2, Kest, correction = "Ripley", verbose = F))
plot(envelope(rslt2, Lest, correction = "Ripley", verbose = F))
plot(density(rslt2), main = "Transect 1")
# plot -------------------------------------------------------------
p1 <- ggplot() +
  geom_density(aes(dist), alpha = 0.3, color = "steelblue", fill = "steelblue") +
  tidybayes::stat_pointinterval(aes(y = 0.00, x = dist), .width = c(.66, .95))
p1 <- p1 + coord_cartesian(ylim = c(0.0, 0.7))
p1 <- p1 + scale_x_continuous(name = "Nearest neighbour distance (m)")
p1 <- p1 + scale_y_continuous(name = "Frequency")
p1
# colony size --------------------------------------------------------------------
# 1 Import data -----------------------------------------------------------
load("./Rdata/2023_palau_transect_all.RData")
# 2 Labelling and wrangling -----------------------------------------------
str(data1)
data1$side <- as.factor(as.character(data1$side))
data1$id <- as.factor(as.character(data1$id))
data1$photo_ID <- as.factor(as.character(data1$photo_ID))
## Wrangling
quantile(data1$mean, c(0.025, 0.5, 0.975), na.rm = T)
sd(data1$mean, na.rm = T)
quantile(data1$mean / 2, c(0.025, 0.5, 0.975), na.rm = T)
sd(data1$mean / 2, na.rm = T)
p1 <- ggplot(data1) +
  geom_density(aes(mean), alpha = 0.3, color = "steelblue", fill = "steelblue") +
  tidybayes::stat_pointinterval(aes(y = 0.00, x = mean), .width = c(.66, .95))
p1 <- p1 + coord_cartesian(ylim = c(0.0, 0.03))
p1 <- p1 + scale_x_continuous(name = "Nearest neighbour distance (m)")
p1 <- p1 + scale_y_continuous(name = "Frequency")
p1
# 2 Labelling and wrangling -----------------------------------------------
str(data1)
data1$side <- as.factor(as.character(data1$side))
data2 <- data1 %>% split(., data1$side)
data2$left$neg <- -data2$left$dist_m_x
df <- data.frame(y = data1$dist_m_y, x = c(data2$right$dist_m_x, data2$left$neg))
grid_size <- 0.01
df$buffer_radius <- rtruncnorm(nrow(df), a = 0.05, mean = 0.12975, sd = 0.06547)
grid_points <- expand.grid(x = seq(min(df$x), max(df$x), by = grid_size), y = seq(min(df$y), max(df$y), by = grid_size))
points_in_buffer <- sapply(1:nrow(grid_points), function(i) {
  any(sqrt((df$x - grid_points$x[i])^2 + (df$y - grid_points$y[i])^2) < df$buffer_radius)
})
plot(0, 0, type = "n", xlim = c(min(df$x), max(df$x)), ylim = c(min(df$y), max(df$y)), xlab = "X", ylab = "Y")
points(grid_points$x, grid_points$y, pch = 15, col = ifelse(points_in_buffer, "red", "black"))
points(df$x, df$y, pch = 19)
## Draw circles around grid points
calculate_nearest_edge_distances <- function(df) {
  distances <- numeric(nrow(df))
  for (i in 1:nrow(df)) {
    other_points <- df[-i, ]
    dist_to_others <- sqrt((df$x[i] - other_points$x)^2 + (df$y[i] - other_points$y)^2)
    nearest_distance <- min(dist_to_others)
    edge_to_edge_distance <- nearest_distance - (df$buffer_radius + df$buffer_radius)
    distances[i] <- max(edge_to_edge_distance, 0)
  }
  return(distances)
}
(nn <- calculate_nearest_edge_distances(df) %>% quantile(., c(0.025, 0.5, 0.975)))
points_in_buffer <- sapply(1:nrow(grid_points), function(i) {
  any(sqrt((df$x - grid_points$x[i])^2 + (df$y - grid_points$y[i])^2) < df$buffer_radius)
})
total_points_in_buffer <- sum(points_in_buffer)
grid_area <- length(seq(0, 10, by = grid_size))^2
(percent_cover <- (total_points_in_buffer / grid_area) * 100)
# compare 2022 with 2023 crest slope --------------------------------------
load("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects/allee_experiments/Rdata/2022_adult_nat_intercol.RData")
df3 <- data.frame(dist = dist, habitat = "slope")
df4 <- data.frame(dist = ad, habitat = "crest")
df5 <- rbind(df3, df4)
intercol_dist_natural = df5
save(intercol_dist_natural, file = file.path("./Rdata", "intercol_dist_natural.RData"))
load("./Rdata/intercol_dist_natural.RData")
p1 <- ggplot(df5, aes(x = dist)) +
  geom_density(aes(group = habitat, color = habitat, fill = habitat), alpha = 0.3)
p1 <- p1 + scale_fill_manual(values = c("steelblue4", "orchid4", "red"), name = "Habitat", labels = c("Reef slope", "Reef crest"), guide = guide_legend(override.aes = list(color = NA))) +
  scale_color_manual(values = c("steelblue4", "orchid4", "steelblue1", "steelblue4", "grey", "grey", "grey", "grey"))
p1 <- p1 + geom_vline(xintercept = 0.707, linetype = "dashed", color = "steelblue4", size = 1)
p1 <- p1 + geom_vline(xintercept = 0.632, linetype = "dashed", color = "orchid4", size = 1)
p1 <- p1 + scale_y_continuous(name = "Density")
p1 <- p1 + scale_x_continuous(name = "Intercolonial distance (m)")
p1 <- p1 + theme_sleek2()
p1 <- p1 + theme(
  legend.position = c(0.80, 0.80),
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 12),
  legend.key.size = unit(1.5, "lines")
)
p1 <- p1 + guides(color = "none")
p1

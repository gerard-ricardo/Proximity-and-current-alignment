# load libraries ----------------------------------------------------------
library(spatstat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(foreach)
library(doParallel)
library(proxy)
library(scales)
library(nls2)
library(dplyr)
library(purrr)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2")
source("./R_scripts/simulate_process.R")
load("./Rdata/palau_unweight_dist.RData")
load("./Rdata/glm_coefficients_for_sim.RData")
den_grid_size <- 100
buffer_radius <- 0.1
grid_size1 <- 0.1
density_fixed <- 0.01
types <- c("random","crest","slope")
set.seed(123)
results <- lapply(types, function(t) {
  simulate_process(density_fixed, den_grid_size, buffer_radius, grid_size1, max_distance, unweight_dist, type = t)
})
par(mfrow = c(1,3))
par(pty = "s")
par(mar = c(4,4,3,1))
for(i in seq_along(results)) {
  res <- results[[i]]
  pts <- res$coral_points
  segs <- res$pairwise_df
  plot(pts$x, pts$y, xlab = "x (m)", ylab = "y (m)", pch = 16, col = "blue",
       xlim = c(0, den_grid_size), ylim = c(0, den_grid_size),
       xaxs = "i", yaxs = "i",
       main = paste0(types[i], "\nN = ", nrow(pts)))
  if(nrow(segs) > 0) {
    segments(x0 = pts$x[match(segs$id1, pts$id)],
             y0 = pts$y[match(segs$id1, pts$id)],
             x1 = pts$x[match(segs$id2, pts$id)],
             y1 = pts$y[match(segs$id2, pts$id)],
             col = scales::alpha("red", 0.4))
  }
  box()
}

## Todo
# load libraries ----------------------------------------------------------
library(sf)
library(lattice)
library(plotKML)
library(RColorBrewer)
library(dartR)
# inputs ------------------------------------------------------------------
resolution <- 0.01
# import gmx --------------------------------------------------------------
# add to grid -------------------------------------------------------------
load("./Rdata/data_gl.RData")
meta = data_gl@other$ind.metrics
meta = meta  %>% rename(id2 = id) %>% mutate(id = paste0("X", id2)) 
meta2 <- meta %>% dplyr::select(c(lat, lon, genotype, total_mean_dia ))
coords <- meta2 %>% rename(x = lon, y = lat, ID = genotype, colony_diam = total_mean_dia) %>% filter(complete.cases(.)) %>% 
  distinct(ID, .keep_all = TRUE)
head(coords)
coords_sf <- st_as_sf(coords, coords = c("x", "y"), crs = 4326)
utm_zone <- 52
crs_utm <- st_crs(sprintf("+proj=utm +zone=%d +datum=WGS84 +units=m +no_defs", utm_zone))
coords_utm <- st_transform(coords_sf, crs = crs_utm)
coords_converted <- as.data.frame(st_coordinates(coords_utm))
grid_dims <- ceiling(as.numeric((sapply(coords_converted, max) - sapply(coords_converted, min))) / resolution) + 1
grid <- array(0, dim = c(10^ceiling(log10(grid_dims[1])), 10^ceiling(log10(grid_dims[2]))))
coords_idx <- as.data.frame(apply(coords_converted, 2, function(col) ceiling((col - min(col)) / resolution) + 1))
min(coords_idx$X)
min(coords_idx$Y)
coords_idx$ID <- coords$ID
coords_idx$colony_diam <- coords$colony_diam
for (i in 1:nrow(coords_idx)) {
  grid[coords_idx[i, "X"], coords_idx[i, "Y"]] <- grid[coords_idx[i, "X"], coords_idx[i, "Y"]] + 1
}
plot(coords_idx$X, coords_idx$Y, main = " Points Distribution", xlab = "X", ylab = "Y", pch = 19)
text(coords_idx$X, coords_idx$Y, labels = coords_idx$ID, pos = 3, cex = 0.7)
## uses large memory
# rotate around centre colony 0101, and extract centre ----------------------------------------
rotation_center <- coords_idx[coords_idx$ID == "c17", c("X", "Y")]
theta <- 130 * (pi / 180)
rotate_coords <- function(coords_df, center, angle_rad) {
  coords_shifted <- coords_df[, c("X", "Y")] - matrix(unlist(center), nrow = nrow(coords_df), ncol = 2, byrow = TRUE)
  rotation_matrix <- matrix(c(cos(angle_rad), -sin(angle_rad), sin(angle_rad), cos(angle_rad)), nrow = 2, byrow = TRUE)
  coords_rotated <- as.matrix(coords_shifted) %*% rotation_matrix
  coords_rotated <- coords_rotated + matrix(unlist(center), nrow = nrow(coords_df), ncol = 2, byrow = TRUE)
  coords_rotated <- round(coords_rotated, 0)
  coords_rotated_df <- data.frame(X = coords_rotated[, 1], Y = coords_rotated[, 2], coords_df[, !(names(coords_df) %in% c("X", "Y"))])
  coords_rotated_df$X = coords_rotated_df$X + abs(min(coords_rotated_df$X))
  coords_rotated_df$Y = coords_rotated_df$Y + abs(min(coords_rotated_df$Y))
  return(coords_rotated_df)
}
coords_rotated <- rotate_coords(coords_df = coords_idx, center = rotation_center, angle_rad = theta)
str(coords_rotated)
coords_rotated <- coords_rotated %>% data.frame()
coords_rotated$ID <- coords_idx$ID
plot(coords_rotated$X, coords_rotated$Y, main = "Rotated Points Distribution", xlab = "X", ylab = "Y", pch = 19)
text(coords_rotated$X, coords_rotated$Y, labels = coords_rotated$ID, pos = 3, cex = 0.7)
unique(coords_rotated)
# extract spe and egg maps ------------------------------------------------
coords_rotated$ID
spoke_1 <- c("1_05", "1_3", "1_10", "1_30")
spoke_2 <- c("2_05", "2_3", "2_10", "2_30")
spoke_3 <- c("3_05", "3_3", "3_10", "3_30")
spoke_4 <- c("4_05", "4_3", "4_10", "4_30")
spoke_5 <- c("5_05", "5_3", "5_10", "5_30")
spoke_6 <- c("6_05", "6_3", "6_10", "6_30")
spoke_7 <- c("7_05", "7_3", "7_10", "7_30")
spoke_8 <- c("8_05", "8_3", "8_10", "8_30")
all_spokes <- c(spoke_1, spoke_2, spoke_3, spoke_4, spoke_5, spoke_6, spoke_7, spoke_8)
spokes <- coords_rotated %>% dplyr::filter(., ID %in% all_spokes)
centre <- coords_rotated %>% dplyr::filter(!(ID %in% all_spokes))
all_colonies = coords_rotated
## to map a single colony
ID2 = "c19"
target <- coords_rotated %>% dplyr::filter(., ID == ID2)
## to map all other or multiple colonies
target$fecund_area = pi * (target$colony_diam / 2)^2 * 0.7
sum(target$fecund_area)
myfun1 <- function(full_patch, subset_patch, resolution) {
  coords1 <- full_patch
  coords2 <- subset_patch
  
  
  offset_m <- 2
  offset <- ceiling(offset_m / resolution)
  grid_dims_x <- 10^ceiling(log10(max(coords1$X))) + offset
  grid_dims_y <- (10^ceiling(log10(max(coords1$Y))) + offset) * 0.8
  grid <- array(0, dim = c(grid_dims_x, grid_dims_y))
  grid_egg <- array(0, dim = c(grid_dims_x, grid_dims_y))
  
  dim(grid)
  coords1$Y <- round((dim(grid)[2] / 2 - max(coords1$Y) / 2) + coords1$Y, 0)
  coords1$X <- coords1$X + (dim(grid)[1] - max(coords1$X) - 2) - offset
  
  
  
  
  spemap <- coords1[coords1$ID %in% coords2$ID, ]
  
  
  for (i in 1:nrow(spemap)) {
    x_idx <- spemap[i, "X"]
    y_idx <- spemap[i, "Y"]
    grid[x_idx, y_idx] <- grid[x_idx, y_idx] + 1
  }
  
  
  n_finecells_colony <- numeric(nrow(spemap))
  
  
  grid_points <- expand.grid(x = seq(1, nrow(grid), by = 1), y = seq(1, ncol(grid), by = 1))
  
  for (i in 1:nrow(spemap)) {
    fecun_zone = 0.7
    colony_diam <- spemap$colony_diam[i]
    if (is.na(colony_diam)) next
    colony_diam_meters <- (colony_diam / 100) * sqrt(fecun_zone)
    buffer_radius <- (colony_diam_meters / 2) / resolution
    x_idx <- spemap[i, "X"]
    y_idx <- spemap[i, "Y"]
    distances <- sqrt((grid_points$x - x_idx)^2 + (grid_points$y - y_idx)^2)
    points_in_buffer <- distances < buffer_radius
    n_finecells_colony[i] <- sum(points_in_buffer)
    grid[grid_points$x[points_in_buffer], grid_points$y[points_in_buffer]] <-
      grid[grid_points$x[points_in_buffer], grid_points$y[points_in_buffer]] + 1
  }
  
  
  message("Fine-grid cells for each colony:")
  print(data.frame(
    ID = spemap$ID,
    Diam_cm = spemap$colony_diam,
    FineCells = n_finecells_colony,
    Actual = pi * (spemap$colony_diam / 2)^2  * 0.7
  ))
  
  res_back <- 1 / resolution
  
  coarse_dim_x <- floor(nrow(grid) / res_back)
  coarse_dim_y <- floor(ncol(grid) / res_back)
  
  grid_coarse_sperm <- matrix(0, nrow = coarse_dim_x, ncol = coarse_dim_y)
  
  
  for (i in 1:coarse_dim_x) {
    for (j in 1:coarse_dim_y) {
      x_indices <- ((i - 1) * res_back + 1):(i * res_back)
      y_indices <- ((j - 1) * res_back + 1):(j * res_back)
      
      block_sum <- sum(grid[x_indices, y_indices])  
      block_sum2 = block_sum/(1/resolution^2) * 100 * 100
      
      grid_coarse_sperm[i, j] <- block_sum2
    }
  }
  
  return(grid_coarse_sperm)
}
out <- myfun1(full_patch = coords_rotated, subset_patch = target, resolution = resolution)
grid_coarse <-  out
coul <- colorRampPalette(brewer.pal(8, "Purples"))(25)
levelplot(t(apply(grid_coarse, 2, rev)),
  col.regions = coul, xlab = "Transverse",
  ylab = "Longitudinal", main = "Conc. (cells/m^3)"
)
sum(grid_coarse)
sum(grid_coarse) - sum(target$fecund_area)
target
(which_cells <- which(grid_coarse > 1, arr.ind = TRUE))
grid_coarse[which_cells]
dim(grid_coarse)
# saved maps --------------------------------------------------------------
ID2
save(grid_coarse, file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects_project/coral_fert_model/Rdata",
                                   "2023palau_coarse_grid_xc19.RData"))
# bind eggmaps to list ----------------------------------------------------
load("./Rdata/2023palau_coarse_grid.RData")
load("./Rdata/2023palau_coarse_grid_centre.RData")
load("./Rdata/2023palau_coarse_grid_spokes.RData")
load("./Rdata/2023palau_coarse_grid_rand_centres.RData")
load("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects_project/coral_fert_model/Rdata/2023palau_coarse_grid_xc13.RData")
xc13 = grid_coarse 
dim(xc13)
load("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects_project/coral_fert_model/Rdata/2023palau_coarse_grid_xc5.RData")
xc5 = grid_coarse 
dim(xc5)
load("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects_project/coral_fert_model/Rdata/2023palau_coarse_grid_xc10.RData")
xc10 = grid_coarse 
dim(xc10)
load("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects_project/coral_fert_model/Rdata/2023palau_coarse_grid_xc19.RData")
xc19 = grid_coarse 
dim(xc19)
rand_center = list('xc13' = xc13, 'xc5' = xc5, 'xc10' = xc10,'xc19' = xc19)
load("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects_project/coral_fert_model/Rdata/2023palau_rand_center_list.RData")
load("./Rdata/2023palau_coarse_grid_0106.RData")
x0106 = grid_coarse 
load("./Rdata/2023palau_coarse_grid_0107.RData")
x0107 = grid_coarse 
load("./Rdata/2023palau_coarse_grid_0108.RData")
x0108 = grid_coarse 
load("./Rdata/2023palau_coarse_grid_0109.RData")
x0109 = grid_coarse 
spoke2 = list('x0106' = x0106, 'x0107' = x0107, 'x0108' = x0108,'x0109' = x0109)
load("./Rdata/2023palau_spoke2_list.RData")
load("./Rdata/2023palau_coarse_grid_0110.RData")
x0110 = grid_coarse 
load("./Rdata/2023palau_coarse_grid_0111.RData")
x0111 = grid_coarse 
load("./Rdata/2023palau_coarse_grid_0112.RData")
x0112 = grid_coarse 
load("./Rdata/2023palau_coarse_grid_0113.RData")
x0113 = grid_coarse 
spoke3 = list('0110' = x0110, '0111' = x0111, '0112' = x0112, '0113' = x0113)
save(spoke3, file = file.path("./Rdata", "2023palau_spoke3_list.RData"))
load("./Rdata/2023palau_spoke3_list.RData")
load("./Rdata/2023palau_coarse_grid_0114.RData")
x0114 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0115.RData")
x0115 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0116.RData")
x0116 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0117.RData")
x0117 = grid_coarse
spoke4 = list('0114' = x0114, '0115' = x0115, '0116' = x0116, '0117' = x0117)
save(spoke4, file = file.path("./Rdata", "2023palau_spoke4_list.RData"))
load("./Rdata/2023palau_spoke4_list.RData")
load("./Rdata/2023palau_coarse_grid_0118.RData")
x0118 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0119.RData")
x0119 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0120.RData")
x0120 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0121.RData")
x0121 = grid_coarse
spoke5 = list('0118' = x0118, '0119' = x0119, '0120' = x0120, '0121' = x0121)
save(spoke5, file = file.path("./Rdata", "2023palau_spoke5_list.RData"))
load("./Rdata/2023palau_spoke5_list.RData")
load("./Rdata/2023palau_coarse_grid_0122.RData")
x0122 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0123.RData")
x0123 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0124.RData")
x0124 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0125.RData")
x0125 = grid_coarse
spoke6 = list('0122' = x0122, '0123' = x0123, '0124' = x0124, '0125' = x0125)
save(spoke6, file = file.path("./Rdata", "2023palau_spoke6_list.RData"))
load("./Rdata/2023palau_spoke6_list.RData")
load("./Rdata/2023palau_coarse_grid_0126.RData")
x0126 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0127.RData")
x0127 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0128.RData")
x0128 = grid_coarse
load("./Rdata/2023palau_coarse_grid_0129.RData")
x0129 = grid_coarse
spoke7 = list('0126' = x0126, '0127' = x0127, '0128' = x0128, '0129' = x0129)
save(spoke7, file = file.path("./Rdata", "2023palau_spoke7_list.RData"))

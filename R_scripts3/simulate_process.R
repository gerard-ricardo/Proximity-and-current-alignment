##notes, changed kappa 
simulate_process <- function(density1 = 0.01, den_grid_size = 100, buffer_radius = 0.1, grid_size1 = 0.1, max_distance = 3, 
                             unweight_dist, type = 'slope') {
  grid_area <- den_grid_size^2
  
  print(paste("density:", density1))
  
  if (type == "random") {
    pp <- rpoispp(density1, win = owin(c(0, den_grid_size), c(0, den_grid_size)))
  
    } else if (type == "crest") {
    mu_crest <- 0.005117504
    scale_crest <- 2189.77211
    kappa_adjusted <- density1/mu_crest
    
    pp <- rMatClust(kappa = kappa_adjusted, r = scale_crest, mu = mu_crest, 
                    win = owin(c(0,den_grid_size),c(0,den_grid_size)))
    
    
  } else if (type == "slope") {
    mu_slope <- 0.08374283
    scale_slope <- 0.1427898
    kappa_adjusted <- density1/mu_slope
    
    pp <- rMatClust(kappa = kappa_adjusted, r = scale_slope, mu = mu_slope, 
                    win = owin(c(0,den_grid_size),c(0,den_grid_size)))
    
    
  } else {
    stop("Invalid type. Choose either 'random', 'crest', or 'slope'.")
  }
  
  
  
  if (is.null(pp) || pp$n == 0) {
    message("Warning: Generated process returned no colonies.")
    return(NULL)
  }
  
  df <- data.frame(x = pp$x, y = pp$y)
  df$id <- seq_len(nrow(df))
  
  print(paste("colonies:", pp$n))
  
  ##Coral cover (not really needed)
  calculate_nearest_edge_distances <- function(df, buffer_radius) {
    distances <- numeric(nrow(df))
    for (i in 1:nrow(df)) {
      other_points <- df[-i, ]
      dist_to_others <- sqrt((df$x[i] - other_points$x)^2 + (df$y[i] - other_points$y)^2)
      nearest_distance <- min(dist_to_others)
      edge_to_edge_distance <- nearest_distance - 2 * buffer_radius
      distances[i] <- max(edge_to_edge_distance, 0)
    }
    return(distances)
  }
  distances <- calculate_nearest_edge_distances(df, buffer_radius)
  
  grid_points <- expand.grid(
    x = seq(0, den_grid_size, by = grid_size1),
    y = seq(0, den_grid_size, by = grid_size1)
  )
  points_in_buffer <- sapply(1:nrow(grid_points), function(i) {
    any(sqrt((df$x - grid_points$x[i])^2 + (df$y - grid_points$y[i])^2) < buffer_radius)
  })
  total_points_in_buffer <- sum(points_in_buffer)
  percent_cover <- (total_points_in_buffer / grid_area) * 100
  
# pdfs --------------------------------------------------------------------
   max_orig_dist <- max(unweight_dist)
  ## Don't normalize to max=1, keep original density shape
  all_pairs <- t(combn(nrow(df), 2))
  dists <- sqrt((df$x[all_pairs[,1]] - df$x[all_pairs[,2]])^2 +
                  (df$y[all_pairs[,1]] - df$y[all_pairs[,2]])^2)
  valid_pairs <- dists <= max_orig_dist
  all_pairs <- all_pairs[valid_pairs,]
  dists <- dists[valid_pairs]
  dx <- df$x[all_pairs[,2]] - df$x[all_pairs[,1]]
  dy <- df$y[all_pairs[,2]] - df$y[all_pairs[,1]]
  angles <- atan2(dy, dx)
  current_dir = 0
  cos_ang <- cos(angles - current_dir)
  
  
#  dist vs angle using glm coefs (might be better than joint pdfs because includes zeros------------------------------------------------------
  beta0_c <- glm_export$beta0_c
  beta_dist <- glm_export$beta_dist
  beta_angle <- glm_export$beta_angle
  mean_dist <- glm_export$mean_dist
  mean_cos_ang <- glm_export$mean_cos_ang
  
  dists_c <- dists - mean_dist
  cos_ang_c <- cos_ang - mean_cos_ang
  
  lambda <- exp(beta0_c + beta_dist * dists_c + beta_angle * cos_ang_c)
  lambda[is.na(lambda)] <- 0
  emb_seq = 7
  total_eggs <- 1e5
  mean_fert_rate_per_dam = 0.328
  sampling_fraction <- emb_seq / (total_eggs * mean_fert_rate_per_dam)  
  lambda_total <- lambda / sampling_fraction
  
  prob_10_5_probs <- 1 - exp(-lambda_total)
  
  successful <- runif(length(dists)) < prob_10_5_probs
  
  pairwise_df_all <- data.frame(
    id1 = df$id[all_pairs[,1]],
    id2 = df$id[all_pairs[,2]],
    distance = dists,
    angle = angles,
    cos_angle = cos_ang,
    effective_strength = lambda_total,
    prob_10_5_probs = prob_10_5_probs,
    success = successful
  )
  pairwise_df_all = pairwise_df_all %>% arrange(., distance)
  
  
  pairwise_df  = pairwise_df_all %>% dplyr::filter(success == "TRUE")
  
  
  pairwise_df = pairwise_df %>% arrange(., distance)
  
  max_simulated_distance <- if(length(dists) == 0) NA_real_ else max(dists)
  
  pairwise_merged <- merge(pairwise_df, df, by.x = "id1", by.y = "id")
  pairwise_merged <- merge(pairwise_merged, df, by.x = "id2", by.y = "id", suffixes = c("1", "2"))
  
  
  return(list(
    distances = distances,
    percent_cover = percent_cover,
    pairwise_df_all = pairwise_df_all,
    pairwise_df = pairwise_df,
    actual_number = pp$n,
    coral_points = df
  ))
}

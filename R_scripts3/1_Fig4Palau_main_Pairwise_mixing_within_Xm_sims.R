## This script is the main test script for the reef type simulations. It loads the simulate_process function. 
## The density of the distance from the 2023 Palau experiment are loaded in to create probabilities around each pair. 
##Adding on extra until 0.03
##maybe faulty logic because of this discretisation....need to think about this more
##Workflow is to:
## check correct function run at bottom
filename = 'slope_26'
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
num_cores <- max(1, min(detectCores() - 1, as.integer(Sys.getenv("SLURM_CPUS_ON_NODE", detectCores() - 1)), 8))
cl <- if(.Platform$OS.type == "unix") makeForkCluster(num_cores) else makeCluster(num_cores)
registerDoParallel(cl)
set.seed(Sys.time())
den_grid_size <- 100
buffer_radius <- 0.1
grid_size1 <- 0.1
max_distance = 3
non_participating_list <- vector()
sims = 201
length = 16
type = 'slope'
den_min = 0.0003
den_max = 0.3
# single run --------------------------------------------------------------
# # SINGLE RUN - ONLY RUN TO HERE --------------------------------------------------------
process_data3 <- function(type = type, den_grid_size = 100, sims = 10, length = 3, den_min = 0.0003, den_max = 0.3) {
  
  
  # inputs ------------------------------------------------------------------
  set.seed(Sys.time())
  den_grid_size <- 100
  
  sims = sims
  length = length
  density_seq = 10^(seq(log10(den_min), log10(den_max), length.out = length))
  buffer_radius <- 0.1
  grid_size1 <- 0.1
  max_distance = 3
  
  non_participating_list <- vector()
  # bootstrap runs ----------------------------------------------------------
  density_results <- vector("list", length(density_seq))
  result_list <- vector("list", sims)
  pb <- txtProgressBar(min = 0, max = length(density_seq), style = 3)
  for (d in seq_along(density_seq)) {
    density <- density_seq[d]
    result_list <- vector("list", sims)
    for (i in 1:sims) {
      result_list[[i]] <- tryCatch({
        simulate_process(density, den_grid_size, buffer_radius, grid_size1, max_distance, unweight_dist, type = type)
      }, error = function(e) {
        message(paste("Skipping failed simulation:", i, "for density:", density))
        NULL
      })
    }
    
    result_list <- result_list[!sapply(result_list, is.null)]
    message(sprintf("Density: %f | Successful runs: %d", density, length(result_list)))
    
    
    if (length(result_list) == 0) {
      density_results[[d]] <- data.frame(
        density = density,
        type = type,
        mean_num_pairs = NA_real_,
        lower_95_num_pairs = NA_real_,
        upper_95_num_pairs = NA_real_,
        mean_non_participating = NA_real_,
        lower_95_non_participating = NA_real_,
        upper_95_non_participating = NA_real_,
        mean_fertilized = NA_real_,
        lower_95_fertilized = NA_real_,
        upper_95_fertilized = NA_real_,
        mean_cv = NA_real_,
        lower_95_cv = NA_real_,
        upper_95_cv = NA_real_
      )
    } else {
      num_pairs_vector <- sapply(result_list, function(res) nrow(res$pairwise_df))
      num_non_participating_vector <- sapply(result_list, function(res) {
        all_ids <- res$coral_points$id
        participating_ids <- unique(c(res$pairwise_df$id1, res$pairwise_df$id2))
        length(setdiff(all_ids, participating_ids))
      })
      fertilised_counts_vector <- sapply(result_list, function(res) {
        res$pairwise_df %>%
          mutate(partner_id1 = id2, partner_id2 = id1) %>%
          pivot_longer(cols = c(id1, id2), names_to = "source", values_to = "id") %>%
          pivot_longer(cols = c(partner_id1, partner_id2), names_to = "partner_source", values_to = "partner_id") %>%
          filter(source != partner_source) %>%
          group_by(id) %>%
          summarise(num_partners = n_distinct(partner_id)) %>%
          pull(num_partners)
      })
      fertilised_counts_vector <- unlist(fertilised_counts_vector)
      cv_eff_strength_vector <- sapply(result_list, function(res) {
        if(is.null(res$pairwise_df_all) || nrow(res$pairwise_df_all) == 0) return(NA_real_)
        s <- res$pairwise_df_all$effective_strength
        if(length(s) == 0 || all(is.na(s)) || mean(s, na.rm=TRUE) == 0) return(NA_real_)
        sd(s, na.rm = TRUE) / mean(s, na.rm = TRUE)
      })
      
      
      
      
    
      density_results[[d]] <- data.frame(
        density = density,
        type = type,
        mean_num_pairs = mean(num_pairs_vector, na.rm = TRUE),
        lower_95_num_pairs = quantile(num_pairs_vector, 0.025, na.rm = TRUE),
        upper_95_num_pairs = quantile(num_pairs_vector, 0.975, na.rm = TRUE),
        mean_non_participating = mean(num_non_participating_vector, na.rm = TRUE),
        lower_95_non_participating = quantile(num_non_participating_vector, 0.025, na.rm = TRUE),
        upper_95_non_participating = quantile(num_non_participating_vector, 0.975, na.rm = TRUE),
        mean_fertilized = mean(fertilised_counts_vector, na.rm = TRUE),
        lower_95_fertilized = quantile(fertilised_counts_vector, 0.025, na.rm = TRUE),
        upper_95_fertilized = quantile(fertilised_counts_vector, 0.975, na.rm = TRUE),
        mean_cv = mean(cv_eff_strength_vector, na.rm = TRUE),
        lower_95_cv = quantile(cv_eff_strength_vector, 0.025, na.rm = TRUE),
        upper_95_cv = quantile(cv_eff_strength_vector, 0.975, na.rm = TRUE)
      )
      
      
      
      if(length(num_pairs_vector) != length(num_non_participating_vector) ||
         length(num_pairs_vector) != length(fertilised_counts_vector)) {
        warning(paste("Skipping per-simulation detailed table for density", density, "due to length mismatch"))
      } else {
        detailed_results_list[[d]] <- data.frame(
          density = rep(density, length(num_pairs_vector)),
          type = rep(type, length(num_pairs_vector)),
          num_pairs = num_pairs_vector,
          num_non_participating = num_non_participating_vector,
          fertilised_counts = fertilised_counts_vector,
          cv_eff_strength = cv_eff_strength
        )
        if(!is.null(detailed_results_list[[d]])) {
          cat("Columns in detailed_results_list[[", d, "]]:", 
              paste(names(detailed_results_list[[d]]), collapse = ", "), "\n")
          cat("Number of rows:", nrow(detailed_results_list[[d]]), "\n")
        } else {
          cat("detailed_results_list[[", d, "]] is NULL\n")
        }
        expected_cols <- c("density","type","num_pairs","num_non_participating","fertilised_counts","cv_eff_strength")
        if(!is.null(detailed_results_list[[d]]) && !all(expected_cols %in% names(detailed_results_list[[d]]))) {
          cat("Mismatch detected at density", density, "\n")
          cat("Missing columns:", paste(setdiff(expected_cols, names(detailed_results_list[[d]])), collapse = ", "), "\n")
        }
        
      }
      
      
    }
    
    setTxtProgressBar(pb, d)
}
  close(pb)
  
  
  # end boot ----------------------------------------------------------------
  
  
  
  
  
  density_summary <- do.call(rbind, density_results)
  rownames(density_summary) <- NULL
  nrow(density_summary)
  
  return(density_summary)
}
  
slope = process_data3(type = type, sims = sims, length = length, den_min = den_min, den_max = den_max)
save(slope, file = file.path("./Rdata", paste0("mixing_", filename, '_', sims, '_', length, ".RData")))
stopCluster(cl)  
# sim end -----------------------------------------------------------------
print(paste0('sims = ', sims))
print(paste0('length = ', length))
print(paste0('type = ', type))

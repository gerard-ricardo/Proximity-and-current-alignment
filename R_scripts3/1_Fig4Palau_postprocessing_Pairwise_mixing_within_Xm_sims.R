library(ggraph)
library(igraph)
library(tidygraph)
library(tidyverse)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2")
load("./Bunya_mixing/mixing_random_24_random_201_16.RData")
den_sum_pois = random
load("./Bunya_mixing/mixing_crest_25_201_16.RData")
den_sum_crest = crest
load("./Bunya_mixing/mixing_slope_26_201_16.RData")
den_sum_slope = slope
## recompile
data3 = rbind(den_sum_pois, den_sum_crest, den_sum_slope)
data3 = rbind(random, crest, slope)
data3$type <- factor(data3$type, levels = c("random", "crest", "slope"))
data3 %>% group_by(type) %>% summarise(runs = nrow(.))
custom_colors <- c("random" = "#7F7F7F",  # Grey
                   "crest" = "#1F77B4",  # Nice blue
                   "slope" = "#D62728")  # Nice red
data_long <- data3 %>% pivot_longer(cols = starts_with("mean_"), names_to = "metric", values_to = "mean") %>%
  mutate(metric = recode(metric,
                         mean_num_pairs = "Pairwise crosses",
                         mean_fertilized = "Colonies/group",
                         mean_non_participating = "Non-part. colonies"))
ci_lower <- data3 %>%
  pivot_longer(cols = starts_with("lower_95_"), names_to = "metric_lower", values_to = "lower") %>%
  mutate(metric = gsub("lower_95_", "", metric_lower)) %>%
  mutate(metric = recode(metric,
                         num_pairs = "Pairwise crosses",
                         fertilized = "Colonies/group",
                         non_participating = "Non-part. colonies"))
ci_upper <- data3 %>%
  pivot_longer(cols = starts_with("upper_95_"), names_to = "metric_upper", values_to = "upper") %>%
  mutate(metric = gsub("upper_95_", "", metric_upper)) %>%
  mutate(metric = recode(metric,
                         num_pairs = "Pairwise crosses",
                         fertilized = "Colonies/group",
                         non_participating = "Non-part. colonies"))
ci_long <- left_join(ci_lower %>% select(density, type, metric, lower),
                     ci_upper %>% select(density, type, metric, upper),
                     by = c("density", "type", "metric"))
data_plot <- left_join(data_long, ci_long, by = c("density", "type", "metric"))
str(data_plot)
facet_labels <- c(
  "Pairwise crosses" = "Pairwise crosses (no.)",
  "Colonies/group" = "Colonies/group (no.)",
  "Non-part. colonies" = "Isolated colonies (no.)"
)
data_plot$metric <- factor(data_plot$metric, levels = c("Pairwise crosses", "Colonies/group", "Non-part. colonies"))
data_plot <- data_plot %>% filter(!is.na(metric)) %>% data.frame()
data_plot$no_col = data_plot$density * 100^2
data_plot$stand <- ifelse(data_plot$metric == "Colonies/group",
                          data_plot$mean,
                          data_plot$mean / data_plot$no_col)
data_plot$lower_stand <- ifelse(data_plot$metric == "Colonies/group",
                                data_plot$lower,
                                data_plot$lower / data_plot$no_col)
data_plot$upper_stand <- ifelse(data_plot$metric == "Colonies/group",
                                data_plot$upper,
                                data_plot$upper / data_plot$no_col)
sim_mix_plot <- ggplot(data_plot, aes(x = density, y = mean, color = type, fill = type)) +
  geom_line(size = 1) +
  geom_point() +
  geom_ribbon(data = data_plot %>% filter(!is.na(lower) & !is.na(upper)),
              aes(ymin = lower, ymax = upper), alpha = 0.3)+
  scale_x_log10() +
  scale_y_continuous(labels = scales::label_number()) +
  scale_color_manual(values = custom_colors, name = "Reef Type") +
  scale_fill_manual(values = custom_colors, name = "Reef Type") +
  labs(x = "Colony density (per m²)", y = NULL) +
  facet_wrap(~ metric,
             scales = "free_y",
             ncol = 1,
             labeller = labeller(metric = facet_labels),
             strip.position = "left",
             switch = "y") +
  theme_sleek2() +
  theme(strip.text = element_text(size = 12),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 90),
        legend.position = c(0.9, 0.2),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA))
sim_mix_plot
head(data_plot)
str(data_plot)
sim_mix_plot2 <- ggplot(data_plot, aes(x = density, y = stand, color = type, fill = type)) +
  geom_line(size = 1) +
  geom_point() +
  geom_ribbon(data = data_plot %>% filter(!is.na(lower_stand) & !is.na(upper_stand)),
              aes(ymin = lower_stand, ymax = upper_stand), alpha = 0.3) +
  scale_x_log10() +
  scale_y_continuous(labels = scales::label_number()) +
  scale_color_manual(values = custom_colors, name = "Reef Type") +
  scale_fill_manual(values = custom_colors, name = "Reef Type") +
  labs(x = "Colony density (per m²)", y = NULL) +
  facet_wrap(~ metric,
             scales = "free_y",
             ncol = 1,
             labeller = labeller(metric = facet_labels),
             strip.position = "left",
             switch = "y") +
  theme_sleek2() +
  theme(strip.text = element_text(size = 12),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 90),
        legend.position = c(0.9, 0.2),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA))
sim_mix_plot2
data_plot <- data_plot %>% filter(metric %in% c("Colonies/group","Non-part. colonies")) %>% data.frame()
facet_labels <- c(
  "Colonies/group" = "Colonies/group (no.)",
  "Non-part. colonies" = "Isolated colonies (prop.)"
)
sim_mix_plot3 <- ggplot(data_plot, aes(x = density, y = stand, color = type, fill = type)) +
  geom_line(size = 1) +
  geom_point() +
  geom_ribbon(data = data_plot %>% filter(!is.na(lower_stand) & !is.na(upper_stand)),
              aes(ymin = lower_stand, ymax = upper_stand), alpha = 0.3) +
  scale_x_log10() +
  scale_y_continuous(labels = scales::label_number()) +
  scale_color_manual(values = custom_colors, name = "Reef Type") +
  scale_fill_manual(values = custom_colors, name = "Reef Type") +
  labs(x = "Colony density (per m²)", y = NULL) +
  facet_wrap(~ metric,
             scales = "free_y",
             ncol = 1,
             labeller = labeller(metric = facet_labels),
             strip.position = "left",
             switch = "y") +
  theme_sleek2() +
  theme(strip.text = element_text(size = 12),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 90),
        legend.position = c(0.9, 0.2),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA))
sim_mix_plot3
# # model fits --------------------------------------------------------------
data_colony <- data3 %>% filter(!is.na(mean_fertilized))
start_grid <- expand.grid(a = seq(min(data_colony$mean_fertilized, na.rm = TRUE),
                                  max(data_colony$mean_fertilized, na.rm = TRUE), length.out = 5),
                          b = seq(0.0001, 0.01, length.out = 5))
fit_nls <- function(start_vals) {
  tryCatch(
    nls(mean_fertilized ~ a * exp(b * density),
        data = data_colony,
        start = list(a = start_vals$a, b = start_vals$b)),
    error = function(e) NULL
  )
}
fits <- start_grid %>% split(1:nrow(.)) %>% map(fit_nls) %>% compact()
if(length(fits) > 0){
  best_model <- fits[[which.min(sapply(fits, AIC))]]
  summary(best_model)
  params <- coef(best_model)
  a_hat <- params["a"]
  b_hat <- params["b"]
  inflection_point <- 1 / b_hat
  cat("Inflection point occurs at density =", inflection_point, "\n")
} else {
  cat("No valid model fits found.\n")
}
density_threshold <- data_colony %>%
  filter(mean_fertilized < 10) %>%
  arrange(density) %>%
  slice(1) %>%
  pull(density)
if(!is.na(density_threshold)){
  cat("Density at which Colonies/group falls below 10:", density_threshold, "\n")
} else {
  cat("No density found where Colonies/group < 10.\n")
}
################### 
# Summary plots for review paper? -----------------------------------------
num_pairs_vector <- sapply(analysis_results, `[[`, "num_pairs")
num_non_participating_vector <- sapply(analysis_results, `[[`, "num_non_participating")
act_num_vector <- sapply(analysis_results, `[[`, "actual_number")
fertilised_counts_vector <- sapply(analysis_results, `[[`, "fertilised_counts")
num_pairs_vector
num_non_participating_vector
med_non_participating_perc = median(num_non_participating_vector / act_num_vector * 100)
p2 = ggplot() +
  geom_density(aes(x = num_pairs_vector), fill = "steelblue", alpha = 0.5) +
  labs(title = "Unique pairwise crosses (no.)", x = paste0("Number of pairs (", sims, " sims)"), y = "Frequency") +
  geom_vline(aes(xintercept = median(num_pairs_vector)), colour = "red", linetype = "dashed", size = 1) +
  annotate("text", x = max(num_pairs_vector), y = 0.9 * max(density(num_pairs_vector)$y),
           label = paste("Median:", round(median(num_pairs_vector), 2)), colour = "black", hjust = 1) +
  theme_minimal()
p3 = ggplot() +
  geom_histogram(aes(x = num_non_participating_vector), bins = length(unique(num_non_participating_vector)), fill = "pink3", alpha = 0.5) +
  labs(title = "Non-participating (indiv.)", x = paste0("Indiv. (", sims, " sims)"), y = "Frequency") +
  geom_vline(aes(xintercept = median(num_non_participating_vector)), colour = "red", linetype = "dashed", size = 1) +
  annotate("text", x = max(num_non_participating_vector)+0.5, y =  max(hist(num_non_participating_vector, plot = FALSE)$counts) * 0.9,
           label = paste("Median:", round(median(num_non_participating_vector), 2)), colour = "black", hjust = 1)+
  annotate("text", x = max(num_non_participating_vector)+0.5, y = max(hist(num_non_participating_vector, plot = FALSE)$counts) * 0.8,
           label = paste("Percent:", round(med_non_participating_perc, 2)), colour = "black", hjust = 1) +
  theme_minimal()
p4 = ggplot() +
  geom_density(aes(x = fertilised_counts$n), fill = "purple2", alpha = 0.5) +
  labs(title = "Colonies per group (no.)", x = paste0("Colonies per group (", sims, " sims)"), y = "Frequency") +
  geom_vline(aes(xintercept = median(fertilised_counts$n)), colour = "red", linetype = "dashed", size = 1) +
  annotate("text", x = max(fertilised_counts$n), y = 0.9 * max(density(fertilised_counts$n)$y),
           label = paste("Median:", round(median(fertilised_counts$n), 2)), colour = "black", hjust = 1) +
  theme_minimal()
den_panel = grid.arrange(
  p1,
  p2,
  p3,
  p4,
  layout_matrix = rbind(
    c(1, 1, 1, 1, 2, 2),
    c(1, 1, 1, 1, 2, 2),
    c(1, 1, 1, 1, 4, 4),
    c(1, 1, 1, 1, 4, 4),
    c(1, 1, 1, 1, 3, 3),
    c(1, 1, 1, 1, 3, 3)
  )
)
den_panel
# # model fits --------------------------------------------------------------

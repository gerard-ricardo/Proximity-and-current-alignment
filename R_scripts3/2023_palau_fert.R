# notes -------------------------------------------------------------------
# 1. Load Libraries ------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(tidyr)
library(glmmTMB)
library(lubridate)
library(gamm4)
library(mgcv)
library(nlme)
library(dplyr)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2")
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek3")
# 1 Import data -----------------------------------------------------------
load("./Rdata/2023_palau_fert.RData")
# 2. Labelling and wrangling -----------------------------------------------
str(data1)
data1$spoke <- as.factor(as.character(data1$spoke))
data1$id <- as.factor(as.character(data1$id))
data1$obs <- factor(formatC(1:nrow(data1), flag = "0", width = 3))
data1$prop <- data1$suc / data1$tot
data1$time <- sprintf("%s PM", data1$time)
data1$time <- as.POSIXct(paste("2000-01-01", data1$time), format="%Y-%m-%d %I:%M %p", tz="UTC")
base_time <- min(data1$time, na.rm = T)
data1$time_from_base <- as.numeric(difftime(data1$time, base_time, units = "secs"))
data1$quality_score = 1
data1 <- data1[!is.na(data1$prop), ]
data1$suc = ifelse(data1$suc == 0, 1, data1$suc)
data1$quality_score[data1$id == "5_05" & data1$spoke == 5] <- 0.5
data1$quality_score[data1$id == "7_05" & data1$spoke == 7] <- 0.5
# 3. Data Exploration ----------------------------------------------------
hist(data1$suc / data1$tot)
centre_indiv <- data1[which(data1$spoke == 'c'),]
radial_indiv <- data1[which(data1$spoke != 'c'),]
radial_indiv <- radial_indiv %>% mutate(deg_from_north = case_when(spoke == 2 ~ 77, spoke == 3 ~ 101, spoke == 4 ~ 124, 
                                                                   spoke == 5 ~ 148, spoke == 6 ~ 173, spoke == 7 ~ 201, TRUE ~ NA_real_))
diff(unique(radial_indiv$deg_from_north))
radial_indiv$deg = radial_indiv$deg_from_north - 146.5
range(radial_indiv$deg )
radial_indiv$deg_rad <- radial_indiv$deg * pi / 180
radial_indiv$sin_deg <- sin(radial_indiv$deg_rad)
radial_indiv$cos_deg <- cos(radial_indiv$deg_rad)
rad_lines = unique(radial_indiv$deg)
## overall mean and weighted means
(with(data1, sum(suc, na.rm = T) / sum(tot, na.rm = T)))
nrow(data1)
with(data1, mean(prop, na.rm = T))
with(data1, sd(prop, na.rm = T))
(wei_mean_prop_cent <- sum(data1$prop * data1$tot, na.rm = TRUE) / sum(data1$tot, na.rm = TRUE))
(wei_mean_prop_cent <- sum(centre_indiv$prop * centre_indiv$tot, na.rm = TRUE) / sum(centre_indiv$tot, na.rm = TRUE))
print(paste0('unweighted centre fert: ' ,mean(centre_indiv$prop)))
set.seed(123)
target_n <- min(centre_indiv$tot)
resampled_data <- centre_indiv %>% group_by(id) %>% mutate(prob_suc = suc / tot) %>%  filter(tot >= target_n) %>%
  mutate(suc1 = rbinom(1, target_n, prob_suc), fail1 = target_n - suc1, tot1 = target_n, prop1 = suc1 / tot1) %>%
  ungroup()%>% data.frame()
mean(resampled_data$prop1)
sd(resampled_data$prop1)
print(paste0('unweighted downsampled centre fert: ' , 100*round(mean(resampled_data$prop1),4)))
set.seed(123)
target_n_a <- 50
resampled_data_all <- data1 %>% group_by(id) %>% mutate(prob_suc = suc / tot) %>%  filter(tot >= target_n_a) %>%
  mutate(suc1 = rbinom(1, target_n, prob_suc), fail1 = target_n - suc1, tot1 = target_n, prop1 = suc1 / tot1) %>%
  ungroup() %>% data.frame()
mean(resampled_data_all$prop1)
print(paste0('unweighted downsampled centre fert: ' , 100 * round(mean(resampled_data_all$prop1), 4)))
median_prop <- median(resampled_data_all$prop1)
iqr_prop <- IQR(resampled_data_all$prop1)
print(paste0('Median fertilisation: ', 100*round(median_prop, 4), '%, IQR: ', 100*round(iqr_prop, 4), '%'))
range(resampled_data_all$prop1)
set.seed(123)
target_n_a <- 50
resampled_data_rad <- radial_indiv %>% group_by(id) %>% mutate(prob_suc = suc / tot) %>%  filter(tot >= target_n_a) %>%
  mutate(suc1 = rbinom(1, target_n, prob_suc), fail1 = target_n - suc1, tot1 = target_n, prop1 = suc1 / tot1) %>%
  ungroup() %>% data.frame()
mean(resampled_data_rad$prop1)
range(resampled_data_rad$prop1)
print(paste0('unweighted downsampled centre fert: ' , 100 * round(mean(resampled_data_rad$prop1), 4))  )
(wei_mean_prop_rad <- sum(radial_indiv$prop * radial_indiv$tot, na.rm = TRUE) / sum(radial_indiv$tot, na.rm = TRUE))
plot(radial_indiv$prop ~ radial_indiv$dist)
plot(radial_indiv$prop ~ radial_indiv$deg)
plot(radial_indiv$prop ~ radial_indiv$time_from_base)
text(radial_indiv$time_from_base, radial_indiv$prop, labels = radial_indiv$id, pos = 3, cex = 0.8, col = "red")
# Modelling ------------------------------------------------------------
## glmm
# GAMM --------------------------------------------------------------------
md1 <- gamm(cbind(suc, tot - suc) ~ s(deg, k = 3) + dist,   random = list(obs = ~1), family = binomial,  method = "REML", verbosePQL = F, 
            data = radial_indiv)
md1$gam
coef(md1)
summary(md1)
summary(md1$gam)
AIC(md1)
coef_dist <- summary(md1$gam)$p.table["dist", "Estimate"]
pval_dist <- summary(md1$gam)$p.table["dist", "Pr(>|t|)"]
odds_ratio <- exp(coef_dist)
print(paste0("For each 1 m increase in distance, the odds of fertilisation decrease by ", 
             round((1 - odds_ratio) * 100, 1), "% (p = ", signif(pval_dist, 3), ")"))
fixef <- coef(md1$gam)
beta0 <- fixef["(Intercept)"]
beta1 <- fixef["dist"]
dist_vals <- seq(0, 10, by = 1)
eta <- beta0 + beta1 * dist_vals
prob <- exp(eta) / (1 + exp(eta))
data.frame(Distance = dist_vals, Probability = prob)
                   
############################
radial_indiv$clonemate_in_sperm <- ifelse(radial_indiv$egg_clone_in_sperm, 1, 0)
radial_indiv$effective_sperm_sources <- ifelse(radial_indiv$id %in% c("5_30", "3_30", "4_30"), 14, 15)
md2 <- gamm(cbind(suc, tot - suc) ~ s(deg, k = 3) + dist + effective_sperm_sources, 
              random = list(obs = ~1), family = binomial, method = "REML", data = radial_indiv)
summary(md2$gam)
AIC(md2)
################
k <- 5
set.seed(123)
radial_indiv$fold <- sample(rep(1:k, length.out = nrow(radial_indiv)))
cv_results <- lapply(1:k, function(i) {
  train_data <- radial_indiv %>% filter(fold != i)
  test_data  <- radial_indiv %>% filter(fold == i)
  model <- gamm(cbind(suc, tot - suc) ~ s(deg, k = 3) + dist, 
                random = list(obs = ~1), family = binomial, 
                method = "REML", data = train_data)
  pred_probs <- predict(model$gam, newdata = test_data, type = "response")
  actual <- test_data$suc / test_data$tot
  log_lik <- sum(actual * log(pred_probs) + (1 - actual) * log(1 - pred_probs), na.rm = TRUE)
  return(log_lik)
})
(mean_log_lik <- mean(unlist(cv_results)))
################
plot(fitted(md1$gam), residuals(md1$gam))
abline(h = 0)
new_data <- expand.grid(
  dist = seq(min(radial_indiv$dist), max(radial_indiv$dist), length.out = 100),
  deg = seq(min(radial_indiv$deg), max(radial_indiv$deg), length.out = 100),
  quality_score = 1  
)
quantile(new_data$predicted)
new_data$obs <- factor("001")
new_data$predicted <- predict(md1$gam, newdata = new_data, type = "response")
radial_indiv$residuals <- residuals(md1$gam, type = "deviance")
p0 <- ggplot(new_data, aes(x = dist, y = deg)) +
  geom_raster(aes(fill = predicted)) +
  scale_fill_gradient(low = "#3B9AB2", high = "#F21A00", breaks = seq(0.1, 0.9, by = 0.1)) +
  labs(x = "Distance from center patch (m)", y = "Downstream of center patch (°)", fill = "Fertilisation \n success") +
  geom_contour(aes(z = predicted), breaks = seq(0.1, 0.9, by = 0.1), color = "white") +
  theme_minimal() +
  scale_y_reverse() + 
  geom_hline(yintercept = rad_lines, color = "#E1AF00", lty = 2) +
  annotate("text", x = max(new_data$dist)+1, y = rad_lines[1]-5, label = "Radial line 2", 
           color = "#E1AF00", size = 5, hjust = 1.2) +  # Adjust position
  annotate("text", x = max(new_data$dist)+1, y = rad_lines[length(rad_lines)]-5, label = "Radial line 7", 
           color = "#E1AF00", size = 5, hjust = 1.2)+
  theme(
    legend.key.height = unit(1.5, "cm"),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12)
  )
p0 <- p0 +
  geom_point(data = radial_indiv, aes(x = dist, y = deg, color = prop), size = 2) +
  geom_text(data = radial_indiv, aes(x = dist, y = deg, label = round(prop, 2)), 
            color = "black", size = 3, vjust = -1)
p0

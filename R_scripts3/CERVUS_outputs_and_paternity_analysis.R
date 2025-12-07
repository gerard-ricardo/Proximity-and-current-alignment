## CERVUS outputs
## PRIORITES (for review)
## Notes
##THOUGHTS ON ASYNCRONY
# 1. Load Libraries ------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(tidyr)
library(ggmap)
library(units)
library(sf)
library(ggmap)
library(Hmisc)
library(purrr)
library(gridExtra)
library(magick)
library(geosphere)
library(plotrix)
library(networkD3)
library(circular)
library(CircStats)
library(dplyr)
library(reshape2)
library(spatstat)
library(glmmTMB)
library(RVAideMemoire)
library(performance)
library(DHARMa)
library(MuMIn)
library(effects)
library(countreg)
library(GGally)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2")
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek1")
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek3")
# Import cervus out data and prep-----------------------------------------------------------
(data1 <- read.csv(file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Palau genetics mixing/Cervus", 
                                    "summary_out_docent1.csv"), check.names = T))
print(paste0('sequenced larvae: ', nrow(data1)))
data1 <- data1 %>% mutate(across(c(Mother.ID , Candidate.father.ID), ~ gsub("05", "0.7", .)))
str(data1)
data1$offsp_id <- as.factor(as.character(data1$Offspring.ID))
data1$moth_id <- as.factor(as.character(data1$Mother.ID))
data1$fath_id <- as.factor(as.character(data1$Candidate.father.ID))
nrow(data1)
(data1 <- data1[grep("\\*", data1$Trio.confidence), ])
print(paste0('sequenced confidence larvae: ', nrow(data1)))
data2 = data1 %>% dplyr::select(., c(offsp_id, moth_id, fath_id))
data2 = data2 %>% mutate(id = moth_id)
data2 %>% group_by(id) %>% summarise(count = n())
# Import meta data -----------------------------------------------------------
meta <- read.csv(file = file.path("./data", "2023_palau_meta_2.csv"))
meta <- meta %>% mutate(across(c(genotype, id), ~ gsub("05", "0.7", .)))
meta = meta  %>% rename(id2 = id) %>% mutate(id = paste0("X", id2)) 
meta2 <- meta %>% dplyr::select(c(id, lat, lon, genotype, total_mean_dia ))
meta2_2 = meta2 %>% na.omit()
str(meta2_2)
unique_meta <- meta2_2 %>% distinct(genotype, .keep_all = TRUE) %>% mutate(genotype = as.character(genotype))
pairs_df <- expand.grid(genotype.x = as.character(unique_meta$genotype),
          genotype.y = as.character(unique_meta$genotype), stringsAsFactors=FALSE) %>% 
          filter(genotype.x < genotype.y)
str(pairs_df)
pairs_df <- pairs_df %>% left_join(unique_meta %>% dplyr::select(genotype,lat,lon, total_mean_dia), by=c("genotype.x"="genotype")) %>% 
          rename(lat.x=lat,lon.x=lon)
pairs_df <- pairs_df %>% left_join(unique_meta %>% dplyr::select(genotype,lat,lon,total_mean_dia), by=c("genotype.y"="genotype")) %>% 
          rename(lat.y=lat,lon.y=lon)
points_x <- st_as_sf(pairs_df,coords=c("lon.x","lat.x"),crs=4326)
points_y <- st_as_sf(pairs_df,coords=c("lon.y","lat.y"),crs=4326)
points_x_proj <- st_transform(points_x,crs=32653)
points_y_proj <- st_transform(points_y,crs=32653)
pairs_df$dist <- st_distance(points_x_proj,points_y_proj,by_element=TRUE)
pairs_df$dist_m <- as.numeric(pairs_df$dist)
ggplot(pairs_df, aes(x=dist_m)) + geom_density(fill="blue",alpha=0.3,colour="blue") + labs(x="Distance (m)",y="Probability Density",title="Pairwise Distance Distribution") + theme_minimal()
quantile(pairs_df$dist_m)
meta2_2 <- meta2 %>% na.omit() %>% distinct(genotype, .keep_all=TRUE) 
meta2_2_sf <- st_as_sf(meta2_2, coords = c("lon", "lat"), crs = 4326)
meta2_2_proj <- st_transform(meta2_2_sf, crs = 32653)
coords <- st_coordinates(meta2_2_proj)
win <- owin(xrange = range(coords[,1]), yrange = range(coords[,2])) 
points_ppp <- ppp(coords[,1], coords[,2], window=win) 
density_map <- density(points_ppp, sigma = bw.diggle(points_ppp))
plot(density_map, main="Spatial Density of Population") 
points(points_ppp, pch=20, col="red")
mean_density <- mean(density_map$v)
join_df1 <- left_join(data2, meta2, by = 'id')
join_df1 <- join_df1 %>% dplyr::select(-id) %>% mutate(id = fath_id)
join_df2 <- left_join(join_df1, meta2, by = 'id')
nrow(join_df2)
join_df2 <- join_df2 %>% mutate(across(c(moth_id, fath_id, id), ~ gsub("05", "0.7", .)))
nrow(join_df2)
join_df2 <- left_join(join_df2, real_geno_df.x, by = "moth_id")
join_df2 <- left_join(join_df2, real_geno_df.y, by = "fath_id")
selfs1 = join_df2[which(join_df2$real_geno.x == join_df2$real_geno.y),]
nrow(selfs1)
nrow(join_df2)
8/102
adult_colonies1 = adult_colonies_sort
adult_colonies1 <- adult_colonies1 %>% mutate(across(c(Individual1, Individual2), ~ gsub("05", "0.7", .)))
adult_colonies1 = adult_colonies1 %>% rename(genotype.x = Individual1, genotype.y = Individual2) %>% 
  mutate(genotype.x = as.factor(genotype.x), genotype.y = as.factor(genotype.y)) 
str(adult_colonies1)
adult_colonies1 <- adult_colonies1 %>%
  mutate(genotype_x_trimmed = sub("_[^_]*$", "", genotype.x))
adult_colonies1 <- adult_colonies1 %>%
  mutate(genotype_y_trimmed = sub("_[^_]*$", "", genotype.y))
genetic_dist_ave <- adult_colonies1 %>%
  group_by(genotype_x_trimmed, genotype_y_trimmed ) %>% summarise(mean_Distance = mean(Distance, na.rm = TRUE)) %>%
  ungroup() %>% rename(genotype.x = genotype_x_trimmed , genotype.y = genotype_y_trimmed, gen_dist = mean_Distance ) %>% 
  data.frame()
## selfing (no selfing from cervus )
## add father genotype and remove multiple options (need group_list2 in 2a) to run
check_cert <- function(position, group_list2) {
  matched_group <- sapply(group_list2, function(x) position %in% x)
  if (any(matched_group)) {
    group_index <- which(matched_group)[1]
    if (length(group_list2[[group_index]]) == 1) {
      return("cert")
    } else {
      return("uncert")
    }
  } else {
    return("uncert")
  }
}
join_df2$cert_status <- mapply(check_cert, join_df2$genotype.y, MoreArgs = list(group_list2 = group_list2))
head(join_df2)
nrow(join_df2)
join_df2 <- join_df2 %>% dplyr::filter(cert_status == "cert")
nrow(join_df2)
print(paste0('non-ambiguos larvae: ', nrow(join_df2)))
larvae_count <- table(join_df2$moth_id)
join_df2$normalised_weight <- 1 / larvae_count[join_df2$moth_id]
load("./Rdata/2023_palau_fert.RData")
data1 <- data1 %>% mutate(across(c(id), ~ gsub("05", "0.7", .)))
join_df2 = join_df2 %>% mutate(id = genotype.x)
join_df2  = left_join(join_df2, data1, by = 'id')
# calculate metrics -----------------------------------------------------
points_x <- st_as_sf(join_df2, coords = c("lon.x", "lat.x"), crs = 4326)
points_y <- st_as_sf(join_df2, coords = c("lon.y", "lat.y"), crs = 4326)
points_x_proj <- st_transform(points_x, crs = 32653)
points_y_proj <- st_transform(points_y, crs = 32653)
join_df2$dist <- st_distance(points_x_proj, points_y_proj, by_element = TRUE)
join_df2$dist_m = as.numeric(join_df2$dist )
##add asycnrony (mum only)
join_df2$minutes <- sub(".*:(\\d+)", "\\1", join_df2$time)  %>% as.numeric()
mean(join_df2$minutes, na.rm = T)
sd(join_df2$minutes, na.rm = T)
min(join_df2$minutes, na.rm = T)
max(join_df2$minutes, na.rm = T)
hist(join_df2$minutes)
range(join_df2$minutes[which(join_df2$minutes > 44)])
join_df2$syncr <- ifelse(join_df2$minutes < 44, "async", "sync")
bin_width <- 20
join_df2 <- join_df2 %>% mutate(angle = (bearing(cbind(lon.x, lat.x), cbind(lon.y, lat.y)) + 360) %% 360)
ds_angle = 323
join_df2 <- join_df2 %>% mutate(ang_rel_ds = (bearing(cbind(lon.x, lat.x), cbind(lon.y, lat.y)) - ds_angle) %% 360, 
                                ang_rel_ds = ifelse(ang_rel_ds > 180, 360 - ang_rel_ds, ang_rel_ds))
join_df2 <- join_df2 %>% mutate(angle_binned = cut(angle, breaks = seq(0, 360, by = bin_width), include.lowest = TRUE, labels = seq(bin_width/2, 360 - bin_width/2, by = bin_width)))
angle_counts_binned <- as.data.frame(table(join_df2$angle_binned))
colnames(angle_counts_binned) <- c("angle_binned", "count")
angle_counts_binned$angle_binned <- as.numeric(as.character(angle_counts_binned$angle_binned))
(polar_plot = polar.plot(lengths = angle_counts_binned$count, polar.pos = angle_counts_binned$angle_binned, 
                         radial.lim = c(0, max(angle_counts_binned$count)), start = 90, lwd = 5, line.col = 4, 
                         clockwise = T, cex.axis = 2, cex = 2,cex.lab = 2))
fig3c_inset_source = angle_counts_binned
save(fig3c_inset_source, file = file.path("./Rdata/fig3c_inset_source.RData"))
create_polar_plot <- function() {
  par(mar = c(0, 0, 0, 0))
  polar.plot(
    lengths = angle_counts_binned$count,
    polar.pos = angle_counts_binned$angle_binned,
    radial.lim = c(0, max(angle_counts_binned$count)),
    start = 90,
    lwd = 4,
    line.col = 4,
    cex.axis = 2, cex = 2,cex.lab = 2,
    clockwise = TRUE
  )
}
center_angle <- 323
angle_range <- 45
lower_bound <- (center_angle - angle_range) %% 360
upper_bound <- (center_angle + angle_range) %% 360
if (lower_bound < upper_bound) {
  filtered_counts <- angle_counts_binned$count[angle_counts_binned$angle_binned >= lower_bound & 
                                                 angle_counts_binned$angle_binned <= upper_bound]
} else {
  filtered_counts <- angle_counts_binned$count[angle_counts_binned$angle_binned >= lower_bound | 
                                                 angle_counts_binned$angle_binned <= upper_bound]
}
percentage <- sum(filtered_counts) / sum(angle_counts_binned$count) * 100
cat("Percentage of counts within ±45° of 323°:", round(percentage, 2), "%\n")
join_df2$angle_rad <- join_df2$angle * pi / 180
angles_circular <- circular(join_df2$angle_rad, units = "radians")
rayleigh.test(angles_circular)
print(paste0('Larvae (no.) after angle and sync calculations: ', nrow(join_df2)))
# Analyses ----------------------------------------------------------------
## Pairwise distances
(quan <- quantile(join_df2$dist_m, probs=c(0, .25, .5, .83, 1)))
unname(quan[3])
print(paste0('unweighted medium: ', unname(quan[3])))
hist(join_df2$dist_m)
table(join_df2$dist_m)
quan_66 <- wtd.quantile(join_df2$dist_m, probs = c(0.17, 0.83))
quan_95 <- wtd.quantile(join_df2$dist_m, probs = c(0.025, 0.975))
lower_66 <- quan_66[1]; upper_66 <- quan_66[2]
lower_95 <- quan_95[1]; upper_95 <- quan_95[2]
unweight_dist = join_df2$dist_m
(quan_w <- wtd.quantile(join_df2$dist_m, weights = join_df2$normalised_weight, probs=c(0, .25, .5, .83, 1)))
print(paste0('weighted medium: ', unname(quan_w[3])))
quan_66_w <- wtd.quantile(join_df2$dist_m, weights = join_df2$normalised_weight, probs = c(0.17, 0.83))
quan_95_w <- wtd.quantile(join_df2$dist_m, weights = join_df2$normalised_weight, probs = c(0.025, 0.975))
lower_66_w <- quan_66_w[1]; upper_66_w <- quan_66_w[2]
lower_95_w <- quan_95_w[1]; upper_95_w <- quan_95_w[2]
cols <- c("Unweighted" = "#A2CFE3",  # pastel blue
          "Weighted" = "#E3A2A8")    # pastel red
fig3c_source = data.frame(dist = join_df2$dist_m, weights = as.vector(join_df2$normalised_weight))
save(fig3c_source, file = file.path("./Rdata/fig3c_source.RData"))
pairwise_dist_plot <- ggplot(join_df2, aes(x = dist_m)) +
  geom_density(aes(fill = "Unweighted", color = "Unweighted"), alpha = 0.3, bw = 3) +
  geom_density(aes(weight = normalised_weight, fill = "Weighted", color = "Weighted"), alpha = 0.3,  bw = 3) +
  geom_errorbarh(aes(y = 0, xmin = lower_66, xmax = upper_66, color = "Unweighted"), 
                 height = 0.0, linetype = "solid", size = 1.5, alpha = 0.5) +
  geom_errorbarh(aes(y = 0, xmin = lower_95, xmax = upper_95, color = "Unweighted"), 
                 height = 0.0, linetype = "solid", size = .5, alpha = 0.5) +
  geom_errorbarh(aes(y = 0, xmin = lower_66_w, xmax = upper_66_w, color = "Weighted"), 
                 height = 0.0, linetype = "solid", size = 1.5, alpha = 0.5) +
  geom_errorbarh(aes(y = 0, xmin = lower_95_w, xmax = upper_95_w, color = "Weighted"), 
                 height = 0.0, linetype = "solid", size = .5, alpha = 0.5) +
  annotate("point", x = unname(quan[3]), y = 0, color = cols["Unweighted"], size = 3, shape = 21, fill = cols["Unweighted"]) +
  annotate("point", x = unname(quan_w[3]), y = 0, color = cols["Weighted"], size = 3, shape = 21, fill = cols["Weighted"]) +
  geom_vline(xintercept = unname(quan[3]), color = '#A2CFE3', lty = 2) + # Unweighted median line
  geom_vline(xintercept = unname(quan_w[3]), color = "#E3A2A8", lty = 2) + # Weighted median line
  scale_fill_manual(values = cols, name = "Density Type") +
  scale_color_manual(values = cols, name = "Density Type") +
  coord_cartesian(ylim = c(0.0, 0.11)) +
  scale_x_continuous(name = "Distance (m)") +
  scale_y_continuous(name = "Probability density") +
  theme_sleek2() +
  theme(
    legend.position = c(0.8, 0.9),
    legend.text = element_text(size = rel(1), colour = "grey20"),
    strip.text.x = element_text(colour = "grey30", size = 8, vjust = -7),
    panel.spacing.y = unit(-1.5, "lines")
  )
pairwise_dist_plot
## weighting for possible pairwise combinations (overrepresetns under values)
dens_est <- density(pairs_df$dist_m,from=min(pairs_df$dist_m),to=max(pairs_df$dist_m),n=512)
dens_func <- approxfun(dens_est$x,dens_est$y,rule=2)
join_df2$weight_distance <- 1/dens_func(join_df2$dist_m)
weighted_quantiles <- wtd.quantile(join_df2$dist_m, weights = join_df2$weight_distance, probs=c(0,0.25,0.5,0.83,1))
weighted_quantiles
join_df2_2 = join_df2 %>%  dplyr::select(c(dist_m, weight_distance))
ggplot(join_df2, aes(x=dist_m)) + geom_density(aes(fill="blue", alpha=0.3, colour="blue", weight = weight_distance))+
  labs(x="Distance (m)",y="Probability Density",title="Pairwise Distance Distribution") + theme_minimal()
join_df2$dist_m
ks.test(join_df2$dist_m, pairs_df$dist_m)
ks.test(join_df2$normalised_weight , pairs_df$dist_m)
## participants
meta3 <- meta2[complete.cases(meta2), ]
part1 = data.frame(id = c(join_df2$moth_id , join_df2$fath_id))
part2 = left_join(part1, meta, by  = 'id') %>% dplyr::select(genotype ) %>%   distinct(genotype)
(nrow(part2))
print(paste0('total positions participants: ', nrow(part2)))
anti_join(meta3, part2,  by = 'genotype' ) %>%   distinct(genotype) %>% nrow()
print(paste0('total non -participants: ', nrow(part2)))
(distinct_sire = join_df2 %>% distinct(genotype.y) %>% nrow())
print(paste0('no sires: ', distinct_sire))
no_dams = join_df2 %>% distinct(genotype.x) %>% nrow()
print(paste0('no dams: ', no_dams))
no_both = join_df2 %>% distinct(genotype.x , genotype.y) %>% nrow() 
print(paste0('no all: ', no_both))
md_sire =join_df2 %>% group_by(genotype.y) %>%
  summarise(count = n_distinct(genotype.x)) %>% summarise(mean_count = mean(count))  %>% round(., 2)
print(paste0('mean dams sired per male: ', md_sire))
mm_d = join_df2 %>% group_by(genotype.x) %>%
  summarise(count = n_distinct(genotype.y)) %>% summarise(mean_count = mean(count))  %>% round(., 2)
print(paste0('mean males per dam: ', mm_d))
data3 = join_df2 %>% dplyr::select(genotype.y, genotype.x) %>% distinct() 
order = data3 %>% count(genotype.y) %>% arrange(desc(n))  %>%  dplyr::select(-n)     
df11 = left_join(order, data3)
top_sires <- df11 %>% count(genotype.y, sort = TRUE) %>% top_n(3, n) %>% pull(genotype.y)
distinct_x_sired_by_top <- df11 %>% filter(genotype.y %in% top_sires) %>% distinct(genotype.x) %>% nrow()
total_distinct_x <- df11 %>% distinct(genotype.x) %>% nrow()
percentage_sired_by_top <- (distinct_x_sired_by_top / total_distinct_x) * 100
percentage_sired_by_top
nrow(data3)
# stats -------------------------------------------------------------------
## 3 types of zero. 1) structural (impossible), 2) missing data (only subsampled), 3) true zeros if sampled whole container
## prepare data for analysis of every pairwise combination as counts. 
unique_dams <- unique(data1$id, na.rm = T)
unique_sires <- unique(na.omit(meta2$genotype))
all_pairs <- expand.grid(genotype.x = unique_dams, genotype.y = unique_sires)
all_pairs$cert_status <- mapply(check_cert, all_pairs$genotype.y, MoreArgs = list(group_list2 = group_list2))
head(all_pairs)
all_pairs <- all_pairs %>% dplyr::filter(cert_status == "cert")
nrow(all_pairs)
observed_crosses_aggregated <- join_df2 %>% group_by(genotype.x, genotype.y, syncr) %>% summarise(count = n(), .groups = 'drop')
print(paste0('After observed_crosses_aggregated: ', nrow(observed_crosses_aggregated)))
complete_data <- left_join(all_pairs, observed_crosses_aggregated, by = c("genotype.x", "genotype.y")) %>% mutate(count = replace_na(count, 0), syncr = replace_na(syncr, "sync"))
meta2_selected <- meta2 %>% mutate(genotype.x = genotype) %>% dplyr::select(genotype.x, lat, lon, total_mean_dia) %>% distinct() %>% na.omit()
meta2_selected_y <- meta2 %>% mutate(genotype.y = genotype) %>% dplyr::select(genotype.y, lat, lon, total_mean_dia) %>% distinct() %>% na.omit()
join_df3 = left_join(complete_data, meta2_selected, by = c("genotype.x"), relationship = "many-to-one") %>% rename(lat.x = lat, lon.x = lon) %>% mutate(syncr = as.factor(syncr))
join_df3 = left_join(join_df3, meta2_selected_y, by = c("genotype.y"), relationship = "many-to-one") %>% rename(lat.y = lat, lon.y = lon)
join_df3 <- join_df3[complete.cases(join_df3), ]
join_df3 = left_join(join_df3, genetic_dist_ave, by = c("genotype.x", "genotype.y"))
head(join_df3)
selfs2 = join_df3 %>% filter(count > 0 & gen_dist <200)
nrow(selfs2) / nrow(join_df2)  * 100
points_x_proj <- st_as_sf(join_df3, coords = c("lon.x", "lat.x"), crs = 4326) %>% st_transform(crs = 32653)
points_y_proj <- st_as_sf(join_df3, coords = c("lon.y", "lat.y"), crs = 4326) %>% st_transform(crs = 32653)
join_df3$dist_m <- as.numeric(st_distance(points_x_proj, points_y_proj, by_element = TRUE))
join_df3 <- join_df3 %>% mutate(ang_rel_ds = (bearing(cbind(lon.x, lat.x), cbind(lon.y, lat.y)) - ds_angle) %% 360, 
                                ang_rel_ds = ifelse(ang_rel_ds > 180, 360 - ang_rel_ds, ang_rel_ds),
                                angle = (bearing(cbind(lon.x, lat.x), cbind(lon.y, lat.y)) + 360) %% 360)
join_df3 <- join_df3 %>% mutate(cos_ang = abs(cos(ang_rel_ds * pi / 180)))
nrow(join_df3)
data4 = data1 %>% dplyr::select(id, suc, prop) %>% na.omit() %>% rename(genotype.x = id)
join_df4 <- join_df3 %>% left_join(data4, by = 'genotype.x') %>% na.omit() %>% dplyr::select(genotype.x, genotype.y, count, suc, syncr, everything())
join_df4  %>% group_by(genotype.x) %>% summarise(count = sum(count), suc = first(suc)) %>% mutate(diff = suc - count) %>% data.frame()
## data exploration
hist(join_df4$dist_m)
plot(join_df4$count~ (join_df4$dist_m))  
plot(join_df4$count~ (join_df4$ang_rel_ds))  
plot(join_df4$count ~ (join_df4$cos_ang ))
join_df4$count
join_df4 = join_df4[which(as.character(join_df4$genotype.x) != as.character(join_df4$genotype.y)),]
str(join_df4)
## truncated count (modelling only the magnitude of success once success is possible/observed)
nrow(join_df4)
df_pos_only <- subset(join_df4, count > 0)
df_pos_only <- df_pos_only %>% mutate(align_cat = case_when(
    (angle >= 281 | ang_rel_ds <= 11) ~ "Aligned",
    (angle >= 101 & ang_rel_ds <= 191) ~ "Aligned",
    TRUE ~ "Not Aligned"
  ))
df_pos_only$align_cat <- factor(df_pos_only$align_cat, levels = c("Not Aligned", "Aligned"))
table(df_pos_only$align_cat)
plot(df_pos_only$count  ~ df_pos_only$dist_m)
plot(df_pos_only$count  ~ df_pos_only$ang_rel_ds)
plot(df_pos_only$count  ~ df_pos_only$cos_ang)
summary(glm(df_pos_only$count  ~ df_pos_only$cos_ang))
range(df_pos_only$ang_rel_ds)
trunc_pois <- zerotrunc(count ~ scale(dist_m) + align_cat, data = df_pos_only,  dist = "poisson") 
summary(trunc_pois)
trunc_nb <- zerotrunc(count ~ scale(dist_m) + cos_ang , data = df_pos_only, dist = "negbin")
summary(trunc_nb) 
AIC(trunc_pois, trunc_nb)
trunc_nb_wei  <- zerotrunc(count ~ scale(dist_m) + scale(ang_rel_ds), data = df_pos_only,
                      dist = "negbin", weights = suc)
summary(trunc_nb_wei)
trunc_pois_offset <- zerotrunc(count ~ scale(dist_m) + cos_ang + offset(log(suc)), 
                               data = df_pos_only, dist = "poisson")
summary(trunc_pois_offset)
trunc_nb_offset <- zerotrunc(count ~ scale(dist_m) + cos_ang + offset(log(suc)), 
                             data = df_pos_only, dist = "negbin")
summary(trunc_nb_offset)
AIC(trunc_nb_wei,trunc_pois_offset, trunc_nb_offset )
## try on colonies sire from central patch
df_pos_only_centre <- df_pos_only[grep('c',df_pos_only$genotype.y),]
trunc_pois <- zerotrunc(count ~ scale(dist_m) + cos_ang, data = df_pos_only_centre,  dist = "poisson") 
summary(trunc_pois)
trunc_nb <- zerotrunc(count ~ scale(dist_m) + cos_ang , data = df_pos_only_centre, dist = "negbin")
summary(trunc_nb) 
AIC(trunc_pois, trunc_nb)
plot(df_pos_only_centre$count  ~ df_pos_only_centre$dist_m)
plot(df_pos_only_centre$count  ~ df_pos_only_centre$cos_ang)
## include zeros and imputations
join_df5 = join_df4 %>% group_by(genotype.x) %>% mutate(
  sum_count = sum(count, na.rm = TRUE),
  final_count = case_when(
    suc == 0 | suc < 5 ~ 0,
    (sum_count == suc & count == 0) ~ 0,
    (sum_count < suc & count == 0) ~ 0,
    TRUE ~ as.numeric(count)
  ),
  zero_type = case_when(
    suc == 0 | suc < 5 ~ "structural_zero", 
    (sum_count == suc & count == 0) ~ "true_zero", 
    (sum_count < suc & count == 0) ~ "missing_zero", 
    count > 0 ~ "positive_count", 
    TRUE ~ NA_character_
  )
) %>% ungroup() %>% data.frame()
table(join_df5$zero_type)
str(join_df5)
join_df5 %>% filter(is.na(zero_type)) %>% group_by(genotype.x) %>% summarise(count_na = n())
plot(final_count ~ dist_m, data = join_df5)
plot(final_count ~ cos_ang, data = join_df5)
# zero-inflated models ----------------------------------------------------
str(join_df5)
join_df5 = join_df5 %>% filter(gen_dist > 200)
join_df5 %>%  dplyr::select(dist_m, cos_ang, total_mean_dia.y , gen_dist) %>% ggpairs()
cor(dplyr::select(join_df5, dist_m, cos_ang, total_mean_dia.y, gen_dist), use = "complete.obs")
join_df5$dist_m_c <- as.vector(scale(join_df5$dist_m, center = TRUE, scale = FALSE))
join_df5$cos_ang_c <- as.vector(scale(join_df5$cos_ang, center = TRUE, scale = FALSE))
cor(dplyr::select(join_df5, dist_m_c, cos_ang_c, total_mean_dia.y, gen_dist), use = "complete.obs")
ggpairs(join_df5 %>% dplyr::select(dist_m_c, cos_ang_c, total_mean_dia.y, gen_dist))
join_df5$weight <- join_df5$sum_count / join_df5$suc
join_df5$weight <- ifelse(join_df5$sum_count == join_df5$suc, 
                          1,
                          join_df5$sum_count / join_df5$suc)
plot(join_df5$final_count ~ join_df5$dist_m_c)
plot(join_df5$final_count ~ join_df5$total_mean_dia.y)
plot(join_df5$final_count ~ join_df5$cos_ang_c)
plot(join_df5$final_count ~ join_df5$gen_dist)
plot(join_df5$final_count ~ join_df5$syncr)
global_mod_nb_zsuc <- glmmTMB(final_count ~ dist_m_c  * cos_ang_c + total_mean_dia.y + poly(gen_dist, 2) + offset(log(suc + 1e-6)) ,  
                        ziformula = ~ suc,
                        family = nbinom2,
                        data = join_df5,
                        na.action = na.fail)
global_mod_nb_z1 <- glmmTMB(final_count ~ dist_m_c  * cos_ang_c + total_mean_dia.y + poly(gen_dist, 2) + offset(log(suc + 1e-6)) ,  
                        ziformula = ~ 1,
                        family = nbinom2,
                        data = join_df5,
                        na.action = na.fail)
global_mod_nb_nozi <- glmmTMB(final_count ~ dist_m_c  * cos_ang_c + total_mean_dia.y + poly(gen_dist, 2) + offset(log(suc + 1e-6)),  
                        family = nbinom2,
                        data = join_df5,
                        na.action = na.fail)
global_mod_poi_zsuc <- glmmTMB(final_count ~ dist_m_c * cos_ang_c + total_mean_dia.y + poly(gen_dist, 2) + offset(log(suc + 1e-6)),
                                       ziformula = ~ suc,
                                       family = poisson,
                                       data = join_df5,
                                       na.action = na.fail)
global_mod_poi_z1 <- glmmTMB(final_count ~ dist_m_c * cos_ang_c + total_mean_dia.y + poly(gen_dist, 2) + offset(log(suc + 1e-6)),
                                             ziformula = ~ 1,
                                             family = poisson,
                                             data = join_df5,
                                             na.action = na.fail)
global_mod_poi_nozi <- glmmTMB(final_count ~ dist_m_c * cos_ang_c + total_mean_dia.y + poly(gen_dist, 2) + offset(log(suc + 1e-6)),
                                      family = poisson,
                                      data = join_df5,
                                      na.action = na.fail)
aic_vals = AIC(global_mod_nb_zsuc, global_mod_nb_z1, global_mod_nb_nozi, global_mod_poi_zsuc, global_mod_poi_z1, global_mod_poi_nozi)
aic_vals[order(aic_vals$AIC), ]
summary(global_mod_nb_nozi)
best_int = global_mod_nb_nozi
performance::check_collinearity(best_int)
sum(resid(best_int, type = "pearson")^2) / (nrow(join_df5) - length(coef(best_int)))
plot(fitted(best_int), resid(best_int))
abline(h = 0)
lines(lowess(fitted(best_int), resid(best_int, type = "pearson")), col = "red")
check_zeroinflation(best_int)
check_singularity(best_int)
sim_res <- simulateResiduals(best_int)
testOutliers(sim_res, type = "bootstrap")
plot(sim_res)
check_model(best_int)
best_add <- glmmTMB(final_count ~ dist_m_c  + cos_ang_c   + total_mean_dia.y + poly(gen_dist,2) + offset(log(suc + 1e-6)) ,  
                    ziformula = ~ 0,
                    family = nbinom1,
                    data = join_df5,
                    na.action = na.fail)
performance::check_collinearity(best_add)
sum(resid(best_add, type = "pearson")^2) / (nrow(join_df5) - length(coef(best_add)))
plot(fitted(best_add), resid(best_add))
abline(h = 0)
lines(lowess(fitted(best_add), resid(best_add, type = "pearson")), col = "red")
check_zeroinflation(best_add)
check_singularity(best_add)
sim_res <- simulateResiduals(best_add)
testOutliers(sim_res, type = "bootstrap")
plot(sim_res)
check_model(best_add)
AIC(global_mod_nb_nozi, best_add)
anova(best_add, global_mod_nb_nozi)
dredged_models <- dredge(global_mod_nb_nozi, rank = "AIC", fixed = ~cond(offset(log(suc + 1e-6))))
(best_models <- subset(dredged_models, delta < 2))
extract_models <- get.models(dredged_models, subset = delta < 2)
len1 = length(extract_models)
best = extract_models[1:len1]
model_avg <- model.avg(best, revised.var = TRUE)
summary(model_avg)
# Extract dist and angle for simulation-------------------------------------------------------------------------
## prepare coefs for export to simulation data
(beta0_c <- summary(model_avg)$coefficients[1,1])
(beta_dist <- summary(model_avg)$coefficients[1,3])
(beta_angle <- summary(model_avg)$coefficients[1,2])
glm_export <- list(
  beta0_c = beta0_c,
  beta_dist = beta_dist,
  beta_angle = beta_angle,
  mean_dist = mean(join_df5$dist_m, na.rm = TRUE),
  mean_cos_ang = mean(join_df5$cos_ang, na.rm = TRUE)
)
# -------------------------------------------------------------------------
## plot effects
##Model averaged plot
## plot effects for model averaged model
pred_data_avg <- with(join_df5, {
  dist_seq <- seq(min(dist_m_c), max(dist_m_c), length.out = 100)
  ang_seq <- seq(min(cos_ang_c), max(cos_ang_c), length.out = 100)
  gen_seq <- seq(min(gen_dist), max(gen_dist), length.out = 100)
  size_seq <- seq(min(total_mean_dia.y), max(total_mean_dia.y), length.out = 100)
  dist_grid <- expand.grid(dist_m_c = dist_seq, cos_ang_c = mean(cos_ang_c), gen_dist = mean(gen_dist), total_mean_dia.y = mean(total_mean_dia.y), suc = mean(suc))
  ang_grid <- expand.grid(dist_m_c = mean(dist_m_c), cos_ang_c = ang_seq, gen_dist = mean(gen_dist), total_mean_dia.y = mean(total_mean_dia.y), suc = mean(suc))
  gen_grid <- expand.grid(dist_m_c = mean(dist_m_c), cos_ang_c = mean(cos_ang_c), gen_dist = gen_seq, total_mean_dia.y = mean(total_mean_dia.y), suc = mean(suc))
  size_grid <- expand.grid(dist_m_c = mean(dist_m_c), cos_ang_c = mean(cos_ang_c), gen_dist = mean(gen_dist), total_mean_dia.y = size_seq, suc = mean(suc))
  list(dist = dist_grid, ang = ang_grid, gen = gen_grid, size = size_grid)
})
get_predictions_avg <- function(newdata, model_avg) {
  pred <- predict(model_avg, newdata = newdata, type = "link", se.fit = TRUE)
  fit <- exp(pred$fit)
  upr <- exp(pred$fit + 1.96 * pred$se.fit)
  lwr <- exp(pred$fit - 1.96 * pred$se.fit)
  data.frame(fit = fit, lwr = lwr, upr = upr)
}
dist_pred_avg <- cbind(pred_data_avg$dist, get_predictions_avg(newdata = pred_data_avg$dist, model_avg = model_avg))
ang_pred_avg <- cbind(pred_data_avg$ang, get_predictions_avg(newdata = pred_data_avg$ang, model_avg = model_avg))
gen_pred_avg <- cbind(pred_data_avg$gen, get_predictions_avg(newdata = pred_data_avg$gen, model_avg = model_avg))
size_pred_avg <- cbind(pred_data_avg$size, get_predictions_avg(newdata = pred_data_avg$size, model_avg = model_avg))
p1_avg <- ggplot(dist_pred_avg, aes(x = dist_m_c)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line(aes(y = fit)) +
  geom_point(data = join_df5, aes(y = final_count), alpha = 0.2) +
  theme_sleek2() +
  labs(x = "Distance (centered)",
       y = "Predicted count",
       title = "A") +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = -0.1))
p2_avg <- ggplot(ang_pred_avg, aes(x = cos_ang_c)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line(aes(y = fit)) +
  geom_point(data = join_df5, aes(y = final_count), alpha = 0.2) +
  theme_sleek2() +
  labs(x = "Cosine angle (centered)",
       y = "Predicted count",
       title = "B") +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = -0.1))
p3_avg <- ggplot(gen_pred_avg, aes(x = gen_dist)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line(aes(y = fit)) +
  geom_point(data = join_df5, aes(y = final_count), alpha = 0.2) +
  theme_sleek2() +
  labs(x = "Genetic distance",
       y = "Predicted count",
       title = "C") +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = -0.1))
p4_avg <- ggplot(size_pred_avg, aes(x = total_mean_dia.y)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line(aes(y = fit)) +
  geom_point(data = join_df5, aes(y = final_count), alpha = 0.2) +
  theme_sleek2() +
  labs(x = "Mean diameter (y)",
       y = "Predicted count",
       title = "D") +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = -0.1))
grid.arrange(
  p1_avg, p2_avg,
  layout_matrix = rbind(c(1, 2)),
  widths = c(1, 1),
  heights = c(1)
)
# # imputation --------------------------------------------------------------
# create map --------------------------------------------------------------
lines_sf <- st_sfc(lapply(1:nrow(join_df2), function(i) {
  st_linestring(rbind(c(join_df2$lon.x[i], join_df2$lat.x[i]),
                      c(join_df2$lon.y[i], join_df2$lat.y[i])))}), crs = 4326)
lines_sf <- st_sf(geometry = lines_sf, dist_m = join_df2$dist_m)
fig3a_source = join_df2 %>% dplyr::select(genotype.x , lat.x, lon.x, genotype.y,lat.y, lon.y) %>%
  rename(position_x = genotype.x, position_y = genotype.y)
save(fig3a_source, file = file.path("./Rdata/fig3a_source.RData"))
distance_plot = ggplot() +
  geom_sf(data = points_x, color = 'blue', size = 2) +
  geom_sf(data = points_y, color = 'red', size = 2) +
  geom_sf(data = lines_sf, aes(color = dist_m), size = 1, alpha = 0.5) +
  scale_color_gradient(low="green", high="red") +
  labs(title="Pairwise Crosses Between Sires and Dams",
       color="Distance (m)") +
  theme_minimal()
## add to google map
lon_range <- range(c(join_df2$lon.x, join_df2$lon.y))
lat_range <- range(c(join_df2$lat.x, join_df2$lat.y))
buffer <- 0.0001
map <- get_googlemap(
  center = c(
    lon = mean(lon_range),
    lat = mean(lat_range)
  ),
  zoom = 21,
  color = "color",
  maptype = "satellite",
  bounds = c(
    left = min(lon_range) - buffer,
    bottom = min(lat_range) - buffer,
    right = max(lon_range) + buffer,
    top = max(lat_range) + buffer
  )
)
                 
base_map <- ggmap(map)
points_y_jittered <- st_jitter(points_y, amount = 0.000002)
points_x_jittered <- st_jitter(points_x, amount = 0.000002)
distance_map_plot <- base_map +
  geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +
  geom_sf(data = points_y_jittered, 
          color = 'red', 
          size = 2.5,           
          alpha = 0.7,          
          inherit.aes = FALSE) + 
  geom_sf(data = points_x, 
          color = 'blue', 
          size = 2.5,           
          alpha = 0.7,          
          inherit.aes = FALSE) + 
  geom_sf(data = lines_sf, 
          aes(color = dist_m), 
          linewidth = 1,
          inherit.aes = FALSE,
          alpha = 0.8) + 
  scale_color_gradient(
    low = "#E6E6FA",    
    high = "#4B0082"    
  ) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  labs(color = "Distance (m)",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(
    legend.position = c(0.8, 0.75),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.margin = margin(5, 5, 5, 5)
  )
  
distance_map_plot
map <- get_googlemap(center = c(lon = join_df2$lon.x[8], lat = join_df2$lat.x[8]), zoom = 20, color = "color", maptype = "satellite")
ggmap(map) +
  geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +
  labs(x = "Longitude", y = "Latitude")
fig1a_source = meta2 %>% dplyr::select(lat, lon, genotype ) %>% rename(position = genotype) %>% tidyr::drop_na() 
save(fig1a_source, file = file.path("./Rdata/palau2023_adult_positions.RData"))
map <- get_googlemap(center = c(lon = join_df2$lon.x[8], lat = join_df2$lat.x[8]), zoom = 20, color = "bw", maptype = "satellite")
p2 = ggmap(map) +
  geom_segment(data = join_df2, aes(x = lon.x, y = lat.x, xend = lon.y, yend = lat.y), color = "grey", size = 1) +
  geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +
  geom_point(data = join_df2, aes(x = lon.x, y = lat.x, color = "Mother"), size = 3) +
  geom_point(data = join_df2, aes(x = lon.y, y = lat.y, color = "Father"), size = 3) +
  scale_color_manual(name = "Parent", values = c("Mother" = "blue", "Father" = "red")) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal()
p2
map <- get_googlemap(center = c(lon = join_df2$lon.x[8], lat = join_df2$lat.x[8]), zoom = 20, color = "bw", maptype = "satellite")
p2 = ggmap(map) +
  geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +
  geom_segment(data = join_df2, mapping = aes(x = lon.x, y = lat.x, xend = lon.y, yend = lat.y), color = "red", size = 1) +
  geom_jitter(data = join_df2, aes(x = lon.x, y = lat.x), color = "blue", size = 3, shape = 16, width = 0.000005, height = 0.000005) +
  geom_jitter(data = join_df2, aes(x = lon.y, y = lat.y), alpha = 0.5, color = "green", size = 3, shape = 17, width = 0.000005, height = 0.000005) +
  labs(x = "Longitude", y = "Latitude") +
  theme_sleek1()
map <- get_googlemap(center = c(lon = join_df2$lon.x[4], lat = join_df2$lat.x[4]), zoom = 22, color = "bw", maptype = "satellite")
p3 = ggmap(map) +
  geom_segment(join_df2, mapping = aes(x = lon.x, y = lat.x, xend = lon.y, yend = lat.y), color = "red", size = 1) +
  geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3)+
  geom_point(data = join_df2, aes(x = lon.x, y = lat.x), color = "blue", size = 3)+
  geom_point(data = join_df2, aes(x = lon.y, y = lat.y), color = "blue", size = 3) +
  labs(x = "Longitude", y = "Latitude") +
  theme_sleek1()
map <- get_googlemap(center = c(lon = join_df2$lon.x[8], lat = join_df2$lat.x[8]), zoom = 20, color = "bw", maptype = "satellite")
(p4_1 = ggmap(map) +
  geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +
  labs(title = "All positions", x = "Longitude", y = "Latitude") +
  theme_sleek1())
map <- get_googlemap(center = c(lon = join_df2$lon.x[8], lat = join_df2$lat.x[8]), zoom = 20, color = "bw", maptype = "satellite")
(p4_2 = ggmap(map) +
    geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +
    geom_point(data = join_df2, aes(x = lon.x, y = lat.x), color = "blue", size = 3)+
    labs(title = "Mother", x = "Longitude", y = "Latitude") +
    theme_sleek1())
map <- get_googlemap(center = c(lon = join_df2$lon.x[8], lat = join_df2$lat.x[8]), zoom = 20, color = "bw", maptype = "satellite")
join_df2_sire <- join_df2 %>% distinct(lon.y, lat.y, .keep_all = TRUE)
(p4_3 = ggmap(map) +
    geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +
    geom_point(data = join_df2_sire, aes(x = lon.y, y = lat.y), color = "red", size = 3, alpha = 0.5)+
    labs(title = "Father", x = "Longitude", y = "Latitude") +
    theme_sleek1())
map <- get_googlemap(center = c(lon = join_df2$lon.x[8], lat = join_df2$lat.x[8]), zoom = 20, color = "bw", maptype = "satellite")
(p4_4 = ggmap(map) +
    geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +
    geom_point(data = join_df2[which(join_df2$dist_m == 0),], aes(x = lon.y, y = lat.y), color = "green", size = 3)+
    labs(title = "Selfing", x = "Longitude", y = "Latitude") +
    theme_sleek1())
grid.arrange(arrangeGrob(p4_1, p4_2, ncol = 2), arrangeGrob(p4_3, p4_4, ncol = 2))
# parental pairwise crosses (based on mother genotypes) -------------------------------------------------------
plot_genotype <- function(genotype) {
  subset_df <- join_df2[join_df2$genotype.x == genotype, ]
  map <- get_googlemap(center = c(lon = join_df2$lon.x[8], lat = join_df2$lat.x[8]), zoom = 20, color = "color", maptype = "satellite")
  p = ggmap(map) +
    geom_point(data = meta2, aes(x = lon, y = lat), color = "white", size = 3) +
    geom_segment(data = subset_df, aes(x = lon.x, y = lat.x, xend = lon.y, yend = lat.y), color = "grey", size = 1) +
    geom_point(data = subset_df, aes(x = lon.y, y = lat.y, color = "Father"), size = 3) +
    geom_point(data = subset_df, aes(x = lon.x, y = lat.x, color = "Mother"), size = 3) +
    scale_color_manual(name = "Parent", values = c("Mother" = "blue", "Father" = "red")) +
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal()+
    theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14), axis.title = element_text(size = 14)) + 
  annotate("text", x = 134.49543, y = 7.3136, label = paste(genotype), hjust = 1.1, vjust = -0.1, size = 8, color = "blue")
  
  return(p)
}
plots <- map(unique(join_df2$genotype.x), plot_genotype)
plots[[1]]
walk2(plots, unique(join_df2$genotype.x), ~ggsave(path = "./plots/pairwise_crosses", 
                      filename = paste0("genotype_", .y, ".jpg"), plot = .x, width = 6, height = 4))
## Add to animation
img_dir <- "./plots/pairwise_crosses"
image_files <- list.files(path = img_dir, pattern = "*.jpg", full.names = TRUE)
images <- image_read(image_files)
gif <- image_animate(images, fps = 0.5)
image_write(gif, path = "./plots/compiled_maps.gif")
image_write(gif, path = "./plots/compiled_maps_docent.gif")
cross_counts <- join_df3 %>%
  group_by(genotype.y) %>%
  summarise(distinct_crosses = n_distinct(genotype.x)) %>%
  arrange(desc(distinct_crosses))  %>%  data.frame()
join_df3 <- join_df3 %>%
  left_join(cross_counts, by = "genotype.y") %>%
  arrange(desc(distinct_crosses), genotype.y)
plot(density(cross_counts$distinct_crosses))
quan = quantile(cross_counts$distinct_crosses)
p3 <- ggplot(cross_counts, aes(x = distinct_crosses)) +
  geom_density(aes(fill = 'steelblue4'), alpha = 0.3) + 
  coord_cartesian(ylim = c(0.0, .75)) +
  scale_fill_manual(values = c("steelblue4", "white", "steelblue1", "white", "grey", "grey")) +
  scale_color_manual(values = c("steelblue4", "grey", "steelblue1", "steelblue4", "grey", "grey", "grey", "grey"))+
  tidybayes::stat_pointinterval(aes(y = 0, x = distinct_crosses), .width = c(.66, .95)) +
  scale_y_continuous(name = "Probability density")+
  scale_x_continuous(name = "Distinct pairwise crosses (no.)") +
  theme_sleek3()+
  geom_point(aes(x = unname(quan[3]), y = 0), color = "black", size = 3, shape = 21, fill = "black")
p3
# breeding units ----------------------------------------------------------
library(tidyverse)
library(tidygraph)
library(ggraph)
library(gridExtra)
mother_list <- unique(join_df2$genotype.x)
mother_list <- sort(mother_list, method = "radix")
mothers_with_c <- grep("^c", mother_list, value = TRUE)
other_mothers <- setdiff(mother_list, mothers_with_c)
sorted_mother_list <- c(mothers_with_c, other_mothers)
spatial_pairs <- map_df(sorted_mother_list, function(mother) {
  mother_data <- join_df2 %>%
    filter(genotype.x == mother) %>%
    select(genotype.x, genotype.y, lat.x, lon.x, lat.y, lon.y, dist_m, angle_rad) %>%
    distinct()
  
  mother_pos <- tibble(
    genotype.x = mother,
    genotype.y = mother,
    lat.x = mother_data$lat.x[1],
    lon.x = mother_data$lon.x[1],
    lat.y = mother_data$lat.x[1],
    lon.y = mother_data$lon.x[1],
    dist_m = 0,
    angle_rad = 0
  )
  
  bind_rows(mother_data, mother_pos)
})
breeding_plots <- map(sorted_mother_list, function(mother) {
  current_data <- spatial_pairs %>%
    filter(genotype.x == mother)
  
  edges <- current_data %>%
    filter(genotype.y != mother) %>%
    select(from = genotype.x, to = genotype.y, lon.x, lat.x, lon.y, lat.y)
  
  nodes <- current_data %>%
    select(name = genotype.y, lon = lon.y, lat = lat.y) %>%
    distinct()
  
  mother_node <- tibble(
    name = mother,
    lon = current_data$lon.x[1],
    lat = current_data$lat.x[1]
  )
  nodes <- bind_rows(nodes, mother_node)
  
  graph <- tbl_graph(
    nodes = nodes,
    edges = edges,
    directed = FALSE
  )
  
  sizee = 0.000015
  p <- ggraph(graph, layout = 'manual', x = nodes$lon, y = nodes$lat) +
    geom_edge_link(alpha = 0.5, color = "grey50") +
    geom_node_point(aes(color = name == mother), size = 2) +
    geom_node_text(aes(label = name), 
                   vjust = -0.3, 
                   hjust = 0.5, 
                   size = 3) +
    scale_color_manual(values = c("red", "blue")) +
    coord_fixed(xlim = c(min(nodes$lon) - sizee, max(nodes$lon) + sizee), 
                ylim = c(min(nodes$lat) - sizee, max(nodes$lat) + sizee)) +
    theme_void() +
    ggtitle(paste("Position:", mother)) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 9),
      plot.margin = margin(0, 0, 0, 0, "mm")
    )
  
  return(p)
})
breeding_plots[[4]]
total_plots <- length(breeding_plots)
n_cols <- 3
n_rows <- 4
layout <- t(matrix(seq_len(n_rows * n_cols), nrow = n_cols, ncol = n_rows))
plot_list <- c(
  lapply(breeding_plots, function(p) {
    if (is.null(p)) {
      grid::nullGrob()
    } else {
      ggplotGrob(p)
    }
  }),
  replicate(n_rows * n_cols - total_plots, grid::nullGrob(), simplify = FALSE)
)
arranged_plots <- grid.arrange(
  grobs = plot_list,
  layout_matrix = layout,
  heights = unit(rep(0.22, n_rows), "npc"),
  widths = unit(rep(0.22, n_cols), "npc"),
  padding = unit(5, "mm")
)
 
# network sankey ----------------------------------------------------------
nrow(df_pos_only)
length(unique(c(df_pos_only$genotype.x, df_pos_only$genotype.y)))
cross_counts1 <- df_pos_only %>%
  group_by(genotype.y) %>%
  summarise(distinct_crosses = n_distinct(genotype.x)) %>%
  arrange(desc(distinct_crosses))  %>%  data.frame()
df_pos_only1 <- df_pos_only %>%
  left_join(cross_counts1, by = "genotype.y") %>%
  arrange(desc(distinct_crosses), genotype.y)
links <- df_pos_only1 %>%
  dplyr::select(genotype.x, genotype.y) %>%
  count(genotype.x, genotype.y) %>%
  mutate(source = paste0(genotype.y, "_s"),
         target = paste0(genotype.x, "_e")) %>%
  dplyr::select(source, target, n) %>%
  arrange(desc(n), target)
sources_ordered <- cross_counts1 %>%
  mutate(source = paste0(genotype.y, "_s")) %>%
  arrange(desc(distinct_crosses)) %>%
  pull(source)
targets_ordered <- links %>%
  arrange(target) %>%
  distinct(target) %>%
  pull(target)
nodes <- data.frame(
  name = c(sources_ordered, targets_ordered) %>%
    unique()
)
links <- links %>%
  mutate(IDsource = match(source, nodes$name) - 1,
         IDtarget = match(target, nodes$name) - 1) %>%
  arrange(IDsource, IDtarget)
sire_ids <- unique(links$IDsource)
nodes$group <- ifelse(1:nrow(nodes) %in% sire_ids, paste0("sire", match(1:nrow(nodes), sire_ids)), "other")
links$group <- nodes$group[links$IDsource + 1]
links$group <- as.factor(links$group)
node_colour_scale <- 'd3.scaleOrdinal()
  .domain(["sire1", "sire2", "sire3", "sire4", "sire5", "sire6", "sire7", "sire8", "sire9", "sire10", "other"])
  .range(["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
          "#ffff33", "#a65628", "#f781bf", "#999999", "#66c2a5", "#cccccc"])' # last colour for 'other'
sankey_plot <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "IDsource",
  Target = "IDtarget",
  Value = "n",
  NodeID = "name",
  NodeGroup = "group",
  LinkGroup = "group",
  fontSize = 17,
  nodeWidth = 10,
  nodePadding = 10,
  width = 800,
  height = 400,
  iterations = 0,
  sinksRight = TRUE,
  colourScale = node_colour_scale
  
)
sankey_plot
# chord diagram -----------------------------------------------------------
library(tidyverse)
library(circlize)
cross_matrix <- join_df2 %>%
  group_by(genotype.x, genotype.y) %>%
  summarise(weight = sum(normalised_weight), .groups = "drop") %>%
  pivot_wider(names_from = genotype.y, values_from = weight, values_fill = 0) %>%
  column_to_rownames(var = "genotype.x") %>%
  as.matrix()
all_sectors <- c(rownames(cross_matrix), colnames(cross_matrix))
sector_colours <- setNames(rainbow(length(all_sectors)), all_sectors)
chordDiagram(
  cross_matrix, 
  transparency = 0.5,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.1), 
  grid.col = sector_colours,
  annotationTrackHeight = c(0.1, 0.02),
  link.sort = TRUE,  
  link.decreasing = FALSE
)
circos.trackPlotRegion(
  track.index = 1, 
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter, 
      CELL_META$ylim[1] + 1, 
      CELL_META$sector.index, 
      facing = "clockwise", 
      niceFacing = TRUE, 
      cex = 0.6
    )
  },
  bg.border = NA
)
# visnetwork --------------------------------------------------------------
library(visNetwork)
library(htmlwidgets)
nodes <- data.frame(id = unique(c(join_df2$genotype.x, join_df2$genotype.y)))
edges <- data.frame(from = join_df2$genotype.y, to = join_df2$genotype.x, value = join_df2$dist_m)
visNetwork = visNetwork(nodes, edges) %>%
  visEdges(arrows = "to") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
# igraph ------------------------------------------------------------------
library(igraph)
library(dplyr)
edge_counts <- join_df3 %>%
  group_by(genotype.x, genotype.y) %>%
  summarise(
    count = n(),
    dist_m = mean(dist_m)
  ) %>%
  ungroup() %>%
  mutate(
    genotype.x = as.character(genotype.x),
    genotype.y = as.character(genotype.y)
  ) %>% data.frame()
nodes <- unique(c(edge_counts$genotype.x, edge_counts$genotype.y))
nodes_df <- data.frame(name = nodes)
g <- graph_from_data_frame(
  d = edge_counts,
  vertices = nodes_df,
  directed = FALSE
)
E(g)$weight <- edge_counts$count
E(g)$dist_m <- edge_counts$dist_m
max_width <- 5
E(g)$width <- (E(g)$weight / max(E(g)$weight)) * max_width
genotypes <- nodes_df$name
n <- length(genotypes)
dist_mat <- matrix(Inf, nrow = n, ncol = n, dimnames = list(genotypes, genotypes))
for (i in 1:nrow(edge_counts)) {
  from <- edge_counts$genotype.x[i]
  to <- edge_counts$genotype.y[i]
  dist <- edge_counts$dist_m[i]
  
  dist_mat[from, to] <- dist
  dist_mat[to, from] <- dist
}
max_dist <- max(dist_mat[is.finite(dist_mat)])
dist_mat[is.infinite(dist_mat)] <- max_dist * 2
mds_result <- cmdscale(dist_mat, k = 2, eig = TRUE)
layout_coords <- mds_result$points
plot(
  g,
  layout = layout_coords,
  vertex.size = 5,
  vertex.label = NA,
  vertex.color = "skyblue",
  edge.width = E(g)$width,
  edge.color = 'grey',
  main = "Network"
)
text(
  x = layout_coords[, 1],
  y = layout_coords[, 2] - 0.0,
  labels = V(g)$name,
  cex = 0.8,
  col = "black"
)
# forecenet - not really showing much -------------------------------------
Nodes <- links %>%
  dplyr::select(name = source, group) %>%
  bind_rows(links %>% dplyr::select(name = target, group)) %>%
  distinct(name, .keep_all = TRUE) %>%
  arrange(name)
Nodes$group <- as.factor(Nodes$group)
Nodes <- Nodes %>%
  mutate(ID = 0:(nrow(Nodes) - 1))
forceNetwork(
  Links = links,
  Nodes = Nodes,
  Source = "IDsource",
  Target = "IDtarget",
  NodeID = "name",
  Group = "group",
  opacity = 0.8,
  zoom = TRUE,
  linkDistance = 100,
  charge = -300,
  fontSize = 12,
  linkWidth = 1.5,
  colourScale = JS("d3.scaleOrdinal()
                    .domain([1,2,3,4,5,6,7,'c'])
                    .range(['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f'])")
)

## notes
library(dbscan)
library(vegan)
library(cluster)
# adult only ---------------------------------------------------------------------
## quick plot
pca_data <- tab(data_gl_adult_unique, freq = TRUE, NA.method = "mean") %>% na.omit()
pca <- dudi.pca(pca_data, center = TRUE, scale = FALSE, nf = 18, scannf = FALSE)
pca_complete <- data.frame(pca$li, pop = data_gl_adult_unique$pop)
set.seed(123)
kmeans_result <- kmeans(pca_data, centers = 2, nstart = 25)
metadata <- data.frame(sample_id = rownames(pca_data), group = as.factor(kmeans_result$cluster))
genetic_dist_matrix1 <- vegdist(pca_data, method = "euclidean")
permanova_result <- adonis2(genetic_dist_matrix1 ~ group, data = metadata, permutations = 9999)
print(permanova_result)
(explained_variance <- pca$eig / sum(pca$eig) * 100)
scree_plot <- data.frame(PC = 1:length(explained_variance), Variance = explained_variance)
ggplot(scree_plot, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_line(aes(y = cumsum(Variance)), group = 1, color = "red") +
  geom_point(aes(y = cumsum(Variance)), color = "red") +
  labs(title = "Scree Plot", x = "Principal Component", y = "Percentage of Variance Explained") +
  theme_sleek2()
set.seed(123)
set.seed(123)
wss <- sapply(1:10, function(k) kmeans(pca_data, centers = k, nstart = 25)$tot.withinss)
plot(1:10, wss, type = "b", pch = 19, frame = FALSE, xlab = "Number of Clusters", ylab = "Total Within-Cluster Sum of Squares")
sil_scores <- sapply(2:10, function(k) mean(silhouette(kmeans(pca_data, centers = k, nstart = 25)$cluster, dist(pca_data))[, 3]))
plot(2:10, sil_scores, type = "b", pch = 19, frame = FALSE, xlab = "Number of Clusters", ylab = "Average Silhouette Score")
kmeans_result <- kmeans(pca_data, centers = 1, nstart = 25)
individuals_in_cluster3 <- which(kmeans_result$cluster == 3)
silhouette_score <- silhouette(kmeans_result$cluster, dist(pca_data))
summary(silhouette_score)
pca_complete$kmeans_cluster <- as.factor(kmeans_result$cluster)
kNNdistplot(pca_data, k = 5)
elbow = 13.5
abline(h = elbow, col = "red", lty = 2)
dbscan_result <- dbscan(pca_data, eps = elbow, minPts = 3)
pca_complete$Cluster_dbscan <- as.factor(dbscan_result$cluster)
t2_dbscan <- ggplot(pca_complete, aes(x = Axis1, y = Axis2, color = Cluster_dbscan)) +
  geom_point(alpha = 0.6) +
  labs(title = paste("PCA Plot with DBSCAN Clusters (eps =", elbow, ")"),
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal()
silhouette_score <- silhouette(dbscan_result$cluster, dist(pca_data))
(dbscan_result)
(summary(silhouette_score))
(t2_dbscan)
# prepare for plotting ----------------------------------------------------
pca_complete2 <- pca_complete %>%
  mutate(
    MumID = str_extract(row.names(pca_complete), "^[^_]+"),
    NewID = paste0('Adu', "_", MumID)
  )
data1 <- dplyr::arrange(pca_complete2, Axis1)
pca_complete2 <- pca_complete2 %>% mutate(across(c(MumID, NewID), as.factor))
str(pca_complete2)
my_palette <- c(
  "dodgerblue", "firebrick", "mediumseagreen", "orchid", "darkorange", "gold",
  "skyblue", "sandybrown", "palevioletred", "mediumturquoise", "khaki",
  "darkslategray", "plum", "lightslategray", "limegreen", "cornflowerblue",
  'tomato', 'pink', 'red'
)
unique_pops <- length(unique(pca_complete2$pop))
if (unique_pops > length(my_palette)) {
  my_palette <- scales::hue_pal()(unique_pops)
}
t2 <- ggplot(pca_complete2, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(fill = factor(pop), color = factor(pop )), shape = 22, size = 3, stroke = 1, alpha = 0.7, position = position_jitter(width = 0.1, height = 0.1)) +
  geom_text_repel(aes(label = pop), color = "grey50", size = 3, max.overlaps = 105, point.padding = 0.5, box.padding = 0.5) +
  scale_fill_manual(values = scales::alpha(my_palette, 0.1)) +
  scale_color_manual(values = my_palette) +
  theme_sleek2() +
  labs(
    x = paste0("PCA1 (", round(explained_variance[1], 2), "%)"),
    y = paste0("PCA2 (", round(explained_variance[2], 2), "%)"),
    color = "Population", fill = "Population", shape = "Stage"
  )
t2
## kmean clustering plot
my_palette <- c("red", "blue")
t2 <- ggplot(pca_complete2, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(fill = factor(kmeans_cluster), color = factor(kmeans_cluster)), shape = 22, size = 3, stroke = 1, alpha = 0.7, position = position_jitter(width = 0.1, height = 0.1)) +
  geom_text_repel(aes(label = pop), color = "grey50", size = 3, max.overlaps = 105, point.padding = 0.5, box.padding = 0.5) +
  scale_fill_manual(values = scales::alpha(my_palette, 0.1)) +
  scale_color_manual(values = my_palette) +
  theme_sleek2() +
  labs(
    x = paste0("PC 1 (", round(explained_variance[1], 2), "%)"),
    y = paste0("PC 2 (", round(explained_variance[2], 2), "%)"),
    color = "K-means Cluster", fill = "K-means Cluster", shape = "Stage"
  )
t2
##dbscan cluster
dbscan_cluster_col <- "Cluster_dbscan"
my_palette <- c("red", "blue")
t2_dbscan <- ggplot(pca_complete2, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(fill = factor(!!sym(dbscan_cluster_col)), color = factor(!!sym(dbscan_cluster_col))), 
             shape = 22, size = 3, stroke = 1, alpha = 0.7, position = position_jitter(width = 0.1, height = 0.1)) +
  geom_text_repel(aes(label = pop), color = "grey50", size = 3, max.overlaps = 105, point.padding = 0.5, box.padding = 0.5) +
  scale_fill_manual(values = scales::alpha(my_palette, 0.1)) +
  scale_color_manual(values = my_palette) +
  stat_ellipse(aes(x = Axis1, y = Axis2, group = !!sym(dbscan_cluster_col), color = factor(!!sym(dbscan_cluster_col))), 
               level = 0.95, linetype = 2, size = 1) +
  theme_sleek2() +
  labs(
    x = paste0("PCA1 (", round(explained_variance[1], 2), "%)"),
    y = paste0("PCA2 (", round(explained_variance[2], 2), "%)"),
    color = "DBSCAN Cluster", fill = "DBSCAN Cluster", shape = "Stage"
  )
t2_dbscan
###################################################################################################################
# adult and larvae --------------------------------------------------------
pca_data <- tab(data_gl_filtered, freq = TRUE, NA.method = "mean") %>% na.omit()
pca <- dudi.pca(pca_data, center = TRUE, scale = FALSE, nf = 2, scannf = FALSE)
pca_complete2 <- data.frame(pca$li, pop = data_gl_filtered$pop, id = data_gl_filtered$ind.names, 
                            genotype = data_gl_filtered$other$ind.metrics$genotype,
                            stage = data_gl_filtered$other$ind.metrics$stage)
pca_complete2 %>% dplyr::filter(stage == 'larvae') %>% group_by(genotype) %>% summarise(n = n()) %>% 
  summarise(mean = mean(n))
(explained_variance <- pca$eig / sum(pca$eig) * 100)
scree_plot <- data.frame(PC = 1:length(explained_variance), Variance = explained_variance)
ggplot(scree_plot, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_line(aes(y = cumsum(Variance)), group = 1, color = "red") +
  geom_point(aes(y = cumsum(Variance)), color = "red") +
  labs(title = "Scree Plot", x = "Principal Component", y = "Percentage of Variance Explained") +
  theme_sleek2()
set.seed(123)
set.seed(123)
kmeans_result <- kmeans(pca_data, centers = 3, nstart = 25)
individuals_in_cluster3 <- which(kmeans_result$cluster == 3)
silhouette_score <- silhouette(kmeans_result$cluster, dist(pca_data))
summary(silhouette_score)
plot(silhouette_score)
pca_complete2$Cluster <- as.factor(kmeans_result$cluster)
pca_complete2 <- pca_complete2 %>%
  mutate(
    MumID = str_sub(stage, 1, 3),
    NewID = paste0(MumID, genotype),
    color = ifelse(grepl("larvae", stage), "2", "1")
  )
my_palette <- c(
  "dodgerblue", "firebrick", "mediumseagreen", "orchid", "darkorange", "gold",
  "skyblue", "sandybrown", "palevioletred", "mediumturquoise", "khaki",
  "darkslategray", "plum", "lightslategray", "limegreen", "cornflowerblue",
  "tomato"
)
unique_pops <- length(unique(pca_complete2$genotype))
if (unique_pops > length(my_palette)) {
  my_palette <- scales::hue_pal()(unique_pops)
}
t2 <- ggplot(pca_complete2, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(fill = genotype, shape = stage, color = color),
             size = 3, stroke = 1, alpha = 0.7, position = position_jitter(width = 0.1, height = 0.1)) +
  geom_text_repel(aes(label = NewID, color = color),
                  size = 3, max.overlaps = 80, point.padding = 0.7, box.padding = 0.6) +
  scale_fill_manual(values = my_palette) +
  scale_color_manual(values = c("1" = "grey", "2" = "lightcoral", "3" = "mediumseagreen", "red" = "red", "black" = "black", "lightcoral", 'grey')) +
  scale_shape_manual(values = c("adults" = 22, "larvae" = 21)) +
  theme_sleek2() +
  labs(
    x = paste0("PCA1 (", round(explained_variance[1], 2), "%)"),
    y = paste0("PCA2 (", round(explained_variance[2], 2), "%)"),
    color = "Cluster", fill = "Population", shape = "Stage")
t2

# clones and genetic relatedness------------------------------------------------------------------
# load libraries ----------------------------------------------------------
library(igraph)
library(ggraph)
library(dplyr)
# from colony -------------------------------------------------------------
text <- '1     1.000      c1_1,c1_2
2     1.000      7_30_1,7_3_2,7_05_2,7_3_1,7_10_2,7_30_2,7_05_1,7_10_1
3     1.000      c7_1,c7_2
4     1.000      c11_1
5     1.000      c16_1,c16_2
6     1.000      c20_1,c20_2
7     1.000      5_3_1,5_3_2,5_05_2,5_05_1,5_10_2
8     1.000      4_3_2,4_10_2,4_05_1,4_10_1,4_3_1,4_05_2
9     1.000      6_30_1,6_30_2,6_3_1,6_05_2,6_3_2,6_10_2,6_10_1,6_05_1
10     1.000      3_3_2,3_3_1,3_10_2,3_05_2,3_10_1,3_05_1
11     1.000      c11_2,c17_1,5_30_1,c17_2,5_30_2
12     1.000      2_30_1,2_10_1,2_30_2
13     1.000      c3_1,c3_2,c13_1,c13_2
14     1.000      c8_1,c8_2
15     1.000      c12_1,c12_2
16     1.000      1_10_1
17     1.000      3_30_2,4_30_,3_30_1,4_30_1,c10_1,c10_2
18     1.000      c9_1,c9_2
19     1.000      c18_1,c18_2
20     1.000      1_05_,c14_1,c14_2,1_30_1
21     1.000      c5_2
22     1.000      c6_1,c6_2
23     1.000      c19_1,c19_2
24     1.000      1_3_
25     1.000      2_05_1,2_3_'
df1 <- read.table(text = text, header = FALSE, sep = "", fill = TRUE, col.names = c("CloneID", "Prob", "group"))
df1$group1 <- sapply(strsplit(df1$group, ","), function(x) x)
df1$group <- NULL
genotype_data_long <- do.call(rbind, lapply(1:nrow(df1), function(i) data.frame(genotype2 = df1$CloneID[i], Prob = df1$Prob[i], id = df1$group[[i]])))
genotype_data_long$genotype2 = paste0('X', genotype_data_long$genotype2)
genotype_data_long = genotype_data_long %>% dplyr::select(-Prob )
str(genotype_data_long)
real_geno_df <- stack(clone_groups)
real_geno_df <- real_geno_df %>% mutate(across(c(values ), ~ gsub("05", "0.7", .)))
real_geno_df <- real_geno_df %>% mutate(values = paste0("X", values))
real_geno_df.x <- real_geno_df %>% rename(moth_id = values, real_geno.x = ind)
real_geno_df.y = real_geno_df %>% rename(fath_id = values, real_geno.y = ind)
data_genind_adult@pop <- factor(rep("population1", nrow(data_genind_adult@tab)))
genetic_dist_matrix <- gd.smouse(data_genind_adult, verbose = TRUE)
genetic_dist_df <- as.data.frame(as.matrix(genetic_dist_matrix))
genetic_dist_df <- tibble::rownames_to_column(genetic_dist_df, "Individual1")
adult_colonies <- pivot_longer(genetic_dist_df, cols = -Individual1, names_to = "Individual2", values_to = "Distance") %>% data.frame()
adult_colonies_sort <- adult_colonies %>% arrange(Individual1, Individual2)
str(adult_colonies_sort)
plot(density(adult_colonies_sort$Distance, main = "Genetic Distance Distribution", xlab = "Genetic Distance",
             ylab = "Frequency"))
heatmap_data <- adult_colonies_sort %>%
  mutate(Individual1 = factor(Individual1, levels = unique(Individual1)),
         Individual2 = factor(Individual2, levels = unique(Individual2)))
ggplot(heatmap_data, aes(x = Individual1, y = Individual2, fill = Distance)) +
  geom_tile() + scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Genetic Relatedness Heatmap", x = "Individual 1", y = "Individual 2", fill = "Genetic Distance")
threshold = 200
unique(adult_colonies_sort$Individual1)
first_group = 'c16_1'
first_group_data <- adult_colonies_sort %>% filter(Individual1 == first_group) %>%
  mutate(Individual2 = factor(Individual2, levels = Individual2[order(Distance)])) 
p1 <- ggplot(first_group_data, aes(x = Individual2, y = Distance)) +
  geom_point() +
  facet_wrap(~ Individual1, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Genetic Distances for Adult Colonies", x = "Individual2", y = "Genetic Distance")
p1
(first_group_data_names <- adult_colonies_sort %>%
    filter(Individual1 == first_group) %>%
    mutate(Individual2 = factor(Individual2, levels = Individual2[order(Distance)])) %>% 
    filter(Distance < threshold)  %>% reframe(Individual2) %>% as.vector())
## batch run for each Indiv1
individual_list <- unique(adult_colonies_sort$Individual1)
result <- map_df(individual_list, function(first_group) {
  first_group_data <- adult_colonies_sort %>%
    filter(Individual1 == first_group) %>%
    mutate(Individual2 = factor(Individual2, levels = Individual2[order(Distance)])) %>%
    filter(Distance < threshold)
  tibble(
    Individual1 = first_group,
    close_rel = list(first_group_data$Individual2)
  )
})
result$close_rel[result$Individual1 == first_group]
result$close_rel
##old grouping useing 200 threshold
edges <- result %>%
  unnest(close_rel) %>%
  filter(!is.na(close_rel))
graph <- graph_from_data_frame(edges, directed = FALSE)
p0 <- ggraph(graph, layout = "fr") +
  geom_edge_link(aes(edge_alpha = 0.5), show.legend = FALSE) +
  geom_node_point(color = "dodgerblue", size = 5) +
  geom_node_text(aes(label = name), vjust = 1.5, hjust = 0.5, size = 3) +
  theme_void()
p0
unique_genotypes <- V(graph)$name
clusters <- components(graph)
group_list <- split(V(graph)$name, clusters$membership)
length(group_list)
group_list <- lapply(group_list, function(x) ifelse(grepl("_$", x), paste0(x, "?"), x))
group_list1 <- lapply(group_list, function(x) substr(x, 1, nchar(x) - 2))
group_list2 <- lapply(group_list1, function(x) {
  unique_elements <- unique(x)
  paste(unique_elements, collapse = " ")
})
group_list2 <- lapply(group_list2, function(x) unlist(strsplit(x, " ")))
Hobs <- function(x) {
  apply(tab(x), 1, function(ind) {
    heterozygous_loci <- sum(ind == 1, na.rm = TRUE)
    non_missing_loci <- sum(!is.na(ind))
    heterozygous_loci / non_missing_loci
  })
}
hetero <- Hobs(data_genind_adult)
filt_hetero <- hetero[hetero > 0.17]
barplot(hetero, las = 2)

# clones and genetic relatedness------------------------------------------------------------------
# load libraries ----------------------------------------------------------
library(igraph)
library(ggraph)
library(dplyr)
# from colony -------------------------------------------------------------

df1 <- read.table(text = text, header = FALSE, sep = "", fill = TRUE, col.names = c("CloneID", "Prob", "group"))
df1$group1 <- sapply(strsplit(df1$group, ","), function(x) x)
df1$group <- NULL
genotype_data_long <- do.call(rbind, lapply(1:nrow(df1), function(i) data.frame(genotype2 = df1$CloneID[i], Prob = df1$Prob[i], id = df1$group[[i]])))
genotype_data_long$genotype2 = paste0('X', genotype_data_long$genotype2)
genotype_data_long = genotype_data_long %>% dplyr::select(-Prob )
str(genotype_data_long)
mlg_analysis <- mlg(data_genind_adult)
print(paste("MLG Analysis:", mlg_analysis, "- clearly incorrect since there are def reps"))
data_genclone_adult <- as.genclone(data_genind_adult)
data_genclone_adult@mlg <- mlg.filter(data_genclone_adult, dist=bitwise.dist, threshold=0.05)
mlg(data_genclone_adult)
mlg_vec <- data_genclone_adult@mlg
mlg_table <- data.frame(Individual=indNames(data_genclone_adult),MLG=mlg_vec)
clone_groups <- split(mlg_table$Individual,mlg_table$MLG)
clone_groups[sapply(clone_groups,length)>0]
genotype_data_long <- stack(clone_groups[sapply(clone_groups, length) > 0])
colnames(genotype_data_long) <- c("id", "group")
group_labels <- setNames(paste0("X", seq_along(unique(genotype_data_long$group))), unique(genotype_data_long$group))
genotype_data_long$genotype2 <- group_labels[as.character(genotype_data_long$group)]
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

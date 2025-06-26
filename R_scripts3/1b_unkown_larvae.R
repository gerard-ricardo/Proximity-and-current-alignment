##DdOCENT FILTERING PLUS ADDED
# Load required libraries ----------------------------------------------------
library(dartR)
library(adegenet)
library(tidyverse)
library(HardyWeinberg)
library(pegas)
library(dartR.popgen)
library(PopGenReport)
library(tictoc)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggrepel)
library(hierfstat)
library(ape)
library(poppr)
library(dbscan)
library(sp)
library(rgdal)
library(clustertend)
library(cluster)
library(plotly)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2")
# Import data ----------------------------------------------------------------
load("./Rdata/2023_Acro_hyac_gl_dDocent.RData")
# Fix labels ------------------------------------------------------------
data_gl_filtered <- data_gl
data_gl_filtered@other$ind.metrics$genotype <- gsub("(?<=_)5(?=$)", "05", data_gl_filtered@other$ind.metrics$genotype,
                                                    perl = TRUE)
data_gl_filtered@other$ind.metrics$id <- gsub("(_5)(?![0-9])", "_05", data_gl_filtered@other$ind.metrics$id, perl = TRUE) 
data_gl_filtered@ind.names <- gsub("(_5)(?![0-9])", "_05", data_gl_filtered@ind.names, perl = TRUE)
tail(data_gl_filtered@other$ind.metrics, 10)
all_id = data_gl_filtered@other$ind.metrics$genotype %>% unique()
centreid = all_id %>% grep('c', .) 
centreid %>% length()
all_id[centreid]
## c16 should be c15 (there is no c16)
data_gl_filtered@other$ind.metrics$genotype <- gsub("c16", "c15", data_gl_filtered@other$ind.metrics$genotype)
data_gl_filtered@other$ind.metrics$genotype %>% unique()
##remove all larvae with c1-c4 (they were not 100% known)
data_gl_filtered
keep_ids <- data_gl_filtered@other$ind.metrics %>%
  filter(
    (stage == 'larvae' & genotype %in% c('x1', 'x2', 'x3', 'x4', 'unkn1', 'unkn2', 'unkn3', 'unkn4')) |
      (stage == 'adults')
  ) %>%
  pull(id)
length(keep_ids)
keep_ids <- unique(keep_ids)
data_gl_filtered <- data_gl_filtered[data_gl_filtered@ind.names %in% keep_ids, ]
tt = arrange(data_gl_filtered@other$ind.metrics, stage , genotype )
data_gl_filtered
# Filter ------------------------------------------------------------------
# Step 4: Remove Monomorphic Loci --------------------------------------------
data_gl_filtered <- gl.filter.monomorphs(data_gl_filtered, verbose = 2)
# Step 1: Remove Individuals with High Missing Data --------------------------
ind_before <- indNames(data_gl_filtered)
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "ind", threshold = 0.5, v = 3)
ind_after <- indNames(data_gl_filtered)
(removed_individuals <- setdiff(ind_before, ind_after))
# Step 2: Remove Loci with High Missing Data ---------------------------------
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "loc", threshold = 0.5, v = 3)
# Step 3: Filter by Minor Allele Frequency (MAF) i.e rare alleles----------------------------
data_gl_filtered <- gl.filter.maf(data_gl_filtered, threshold = 0.05, v = 3)
mean(data_gl_filtered@other$loc.metrics$maf)
# Step 5: Filter Loci by Read Depth ------------------------------------------
data_gl_filtered <- gl.filter.rdepth(data_gl_filtered, lower = 3, v = 3)
data_gl_filtered <- gl.recalc.metrics(data_gl_filtered, v = 3)
# Step 6: Apply Stricter Filtering on Loci -----------------------------------
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "loc", threshold = 0.95, v = 3)
data_gl_filtered <- gl.filter.rdepth(data_gl_filtered, lower = 20, v = 3)
# Step 7: Remove Loci with Excessively High Depth ----------------------------
mean_depth <- mean(data_gl_filtered$other$loc.metrics$rdepth)
sd_depth <- sd(data_gl_filtered$other$loc.metrics$rdepth)
max_depth <- mean_depth + 3 * sd_depth
min_depth <- mean_depth - 3 * sd_depth
data_gl_filtered <- gl.filter.rdepth(data_gl_filtered, lower = min_depth, upper = max_depth, v = 3)
# Step 8: Remove Secondary Loci ----------------------------------------------
data_gl_filtered <- gl.filter.secondaries(data_gl_filtered, method = "random", v = 3)
gl.report.reproducibility(data_gl_filtered)
data_gl_filtered <- gl.filter.reproducibility(data_gl_filtered, t = 0.95, v = 3)
# Step 9: Filter Loci Deviating from Hardy-Weinberg Equilibrium (HWE) -------
data_gl_filtered <- gl.filter.hwe(
  x = data_gl_filtered,
  subset = "each",
  n.pop.threshold = 1,
  test.type = "Exact",
  mult.comp.adj = TRUE,
  mult.comp.adj.method = "fdr",
  alpha = 0.01,
  n.min = 5,
  verbose = 2
)
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "loc", threshold = 0.99, v = 3)
data_gl_filtered <- gl.filter.maf(data_gl_filtered, threshold = 0.07, v = 3)
mean(data_gl_filtered@other$loc.metrics$maf)
# Final Recalculation of Metrics ---------------------------------------------
data_gl_filtered <- gl.recalc.metrics(data_gl_filtered, v = 3)
data_gl_filtered@other$ind.metrics$genotype <- gsub("(?<=_)5(?=$)", "05", data_gl_filtered@other$ind.metrics$genotype,
                                                    perl = TRUE)
data_gl_filtered@ind.names <- gsub("(_5)(?![0-9])", "_05", data_gl_filtered@ind.names, perl = TRUE)
data_gl_filtered@other$ind.metrics %>% group_by(stage) %>% summarise(n = n()) 
data_genind_pre <- gl2gi(data_gl_filtered)
adults_indices <- which(data_gl_filtered@other$ind.metrics$stage == "adults")
data_gl_filtered_adult <- data_gl_filtered[adults_indices, ]
data_gl_filtered_adult@other$ind.metrics$stage <- droplevels(data_gl_filtered_adult@other$ind.metrics$stage)
data_gl_filtered_adult@ind.names
data_genind_adult <- gl2gi(data_gl_filtered_adult)
# cervus prep -------------------------------------------------------------
# Cervus Platy mapping extraction (working)------------------------------------------------------
##
####NOTE: Seems to be too many canidate reps n offspring file
data1 <- data_genind_pre@tab
head(data1)
data1 = t(data1) %>% data.frame()
colnames(data1) <- gsub("\\.", "_", colnames(data1))
colnames(data1) <- ifelse(grepl("^X", colnames(data1)), colnames(data1), paste0("X", colnames(data1)))
data1$rownames <- rownames(data1)
data2 <- data1 %>% mutate(CloneID = sub("-.*", "", rownames), allele = sub(".*\\.(.)$", "\\1", rownames)) %>% 
  dplyr::select(CloneID, allele, everything()) %>% dplyr::select(-rownames) %>% arrange(., CloneID)
data1_long = data2 %>% tidyr::pivot_longer(-c(CloneID, allele) ,  names_to = "id" ,values_to = "counts") %>% arrange(., CloneID, id) %>% data.frame()
nrow(data1_long)
row_counts = data1_long %>%
  group_by(CloneID, id) %>%
  summarise(row_count = n()) %>%
  ungroup()
min(row_counts$row_count)
data1_long = data1_long %>%
  group_by(CloneID, id) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  data.frame()
nrow(data1_long)
any(data1_long == "", na.rm = TRUE)
repeat_with_separator <- function(allele, counts, sep = "_") {
  if (is.na(allele) | is.na(counts)) {
    return("NA")
  } else {
    return(paste(rep(allele, counts), collapse = sep))
  }
}
data1_long <- data1_long %>%
  mutate(allele_multiplied = mapply(repeat_with_separator, allele, counts, sep = ","))
any(data1_long == "", na.rm = TRUE)
summary_data <- data1_long %>%
  group_by(CloneID, id) %>%
  summarise(
    id = first(id),
    base = paste(allele_multiplied, collapse = ",")
  ) %>%
  arrange(id) %>%
  data.frame()
summary_data$base <- gsub("^,+|,+$", "", summary_data$base)
summary_data <- summary_data %>%
  separate(base, into = c("a", "b"), sep = ",", extra = "drop", fill = "right")
data1_long = summary_data %>% tidyr::pivot_longer(-c(CloneID, id) ,  names_to = "base" ,values_to = "code") %>% arrange(., CloneID, id) %>% data.frame()
data1_long$LocusID = paste0(data1_long$CloneID, data1_long$base)
result <- data1_long %>%
  group_by(CloneID) %>%
  summarise(missing_prop = sum(code == "NA") / n())
hist(result$missing_prop)
quantile(result$missing_prop)
nrow(data1_long)
length(unique(data1_long$CloneID))
good_loci = result %>% filter(missing_prop < 0.01)
nrow(good_loci)
print(paste('You have', nrow(good_loci), 'good loci'))
data1_long <- data1_long %>%
  semi_join(good_loci, by = "CloneID")
length(unique(data1_long$CloneID))
data1_long = data1_long %>% dplyr::select(-c(CloneID , base))
data1_long <- data1_long %>%
  mutate(code = ifelse(code == "NA", "*", code))
data_wide <- data1_long %>% tidyr::pivot_wider(names_from = LocusID, values_from = code, names_prefix = "X") %>% 
  data.frame()
write.csv(data_wide, row.names = FALSE,
          file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Palau genetics mixing/Cervus",
                           "a_hyac_map_letters_code1_unknowns.csv"))
# offspring file ----------------------------------------------------------
df2 = data_gl_filtered$other$ind.metrics
df2$id = paste0('X', df2$id)
mismatches <- setdiff(colnames(data1), df2$id)
df2_lar = subset(df2, stage!= 'adults')
df2_adult = subset(df2, stage!= 'larvae')
df1a = df2_lar %>% dplyr::select(., id, genotype )
offspring_ids <- df2_lar %>% dplyr::select(id)
candidate_parents <- unique(df2_adult$id)
num_candidates <- length(candidate_parents)
offspring_df <- offspring_ids
for (i in 1:num_candidates) {
  offspring_df[[paste0("candidate_", i)]] <- candidate_parents[i]
}
offspring_df
write.csv(offspring_df, row.names = FALSE,
          file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Palau genetics mixing/Cervus",
                           "offspring_a_hyac_both_unknown.csv"))

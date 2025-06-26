## Adds genotype numbers to meta data
# Notes for standard filtering --------------------------------------------
# Notes for dDocent filtering ---------------------------------------------
# load libraries ----------------------------------------------------------
library(poppr)
library(magrittr)
library(adegenet) 
library(pegas)  
library(dplyr)
library(purrr)
library(tidyverse)
library(igraph)
library(ggraph)
# attach clone id (genotype to data) --------------------------------------
data_gl_filtered_adult@other$ind.metrics$id <- as.character(data_gl_filtered_adult@other$ind.metrics$id)
genotype_data_long$id <- as.character(genotype_data_long$id)
str(data_gl_filtered_adult@other$ind.metrics)
data_gl_filtered_adult@other$ind.metrics <- data_gl_filtered_adult@other$ind.metrics %>% left_join(genotype_data_long,  by = "id")
data_gl_filtered_adult@other$ind.metrics$genotype2
tt = data_gl_filtered_adult@other$ind.metrics
# population filtering and objects ----------------------------------------
data_gl_filtered_adult@pop <- as.factor(data_gl_filtered_adult@other$ind.metrics$genotype2)
data_gl_filtered_adult
## temporarily remove missing genotypes
data_genind <- gl2gi(data_gl_filtered)
data_genind_adult <- gl2gi(data_gl_filtered_adult)
genotype_matrix <- data_genind_adult@tab
(callrate <- rowMeans(!is.na(genotype_matrix)))
ind_names <- indNames(data_genind_adult)
genotypes <- data_genind_adult@other$ind.metrics$genotype2
geno_df <- data.frame(individual = ind_names, genotype = genotypes, callrate = callrate, stringsAsFactors = FALSE)
best_geno_df <- geno_df %>% group_by(genotype) %>% slice_max(order_by = callrate, n = 1, with_ties = FALSE) %>% 
  ungroup() %>% data.frame()
best_ind_names <- best_geno_df$individual
best_indices <- match(best_ind_names, indNames(data_genind_adult))
data_genind_adult_unique <- data_genind_adult[best_indices, ]
data_gl_adult_unique = gi2gl(data_genind_adult_unique, parallel = FALSE, verbose = NULL)
data_gl_adult_unique@pop

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
data_gl_filtered@other$ind.metrics$genotype <- gsub("(?<=_)5(?=$)", "_.7", data_gl_filtered@other$ind.metrics$genotype,
                                                    perl = TRUE)
data_gl_filtered@other$ind.metrics$id <- gsub("(_5)(?![0-9])", "_.7", data_gl_filtered@other$ind.metrics$id, perl = TRUE) 
data_gl_filtered@ind.names <- gsub("(_5)(?![0-9])", "_.7", data_gl_filtered@ind.names, perl = TRUE)
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
  filter(!(stage == 'larvae' & genotype %in% c('x1', 'x2', 'x3', 'x4'))) %>% pull(id)
keep_ids <- unique(keep_ids)
data_gl_filtered <- data_gl_filtered[data_gl_filtered@ind.names %in% keep_ids, ]
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
# Blast prep --------------------------------------------------------------

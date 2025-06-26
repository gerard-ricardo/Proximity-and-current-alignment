# COLONY formatting -------------------------------------------------------
##NOTES: This creates files for COLONY. The error rates are first run with priors 'marker_info_...'. Then after it has
data_genind_adult
data_genind_adult
error_rate = 0.01
colony_data <- function(genind_obj) {
  
  input_name <- deparse(substitute(genind_obj))
  
  ## Calculating allelic dropout rate - these estimates are off because of subpop structure
  
  hist(genind_obj@other$loc.metrics$AvgPIC )
  
  genotype_matrix <- as.matrix(tab(genind_obj))
  ncol(genotype_matrix)
  column_names <- colnames(genotype_matrix)
  locus_identifiers <- sapply(strsplit(column_names, "/"), `[`, 1)
  locus_counts <- table(locus_identifiers)
  valid_loci <- names(locus_counts[locus_counts == 2])
  valid_columns <- locus_identifiers %in% valid_loci
  genotype_matrix <- genotype_matrix[, valid_columns]
  ncol(genotype_matrix)
  num_columns_to_keep = ncol(genotype_matrix)
  
  
  
  ##  make marker file ##
  (num_loci <- ncol(genotype_matrix)  /2)
  marker_names <- paste0("locus", 1:num_loci)
  
  marker_types <- rep(0, num_loci)
  
  allelic_dropout_rates <- rep(0.01, num_loci)
  genotyping_error_rates <- rep(error_rate, num_loci)
  
  marker_info <- rbind(
    marker_names,
    marker_types,
    allelic_dropout_rates,
    genotyping_error_rates
  )
  
  
  write.table(marker_info, file = paste0("./data/marker_info_", input_name, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
  
  marker_info
  ncol(marker_info)
  
  ## make offspring file ##
  
  individual_names <- indNames(genind_obj)
  
  colony_data_lines <- character(nrow(genotype_matrix))
  for (i in 1:nrow(genotype_matrix)) {
    individual_data <- c(individual_names[i])
    
    for (j in seq(1, ncol(genotype_matrix), by = 2)) {
      allele1_dosage <- genotype_matrix[i, j]
      allele2_dosage <- genotype_matrix[i, j + 1]
      
      if (is.na(allele1_dosage) || is.na(allele2_dosage)) {
        alleles <- c(0, 0)
      } else if (allele1_dosage == 2 && allele2_dosage == 0) {
        alleles <- c(1, 1)
      } else if (allele1_dosage == 1 && allele2_dosage == 1) {
        alleles <- c(1, 2)
      } else if (allele1_dosage == 0 && allele2_dosage == 2) {
        alleles <- c(2, 2)
      } else {
        alleles <- c(0, 0)
      }
      
      individual_data <- c(individual_data, alleles)
    }
    
    colony_data_lines[i] <- paste(individual_data, collapse = " ")
  }
  
  writeLines(colony_data_lines, con = paste0("./data/col_input_", input_name, ".txt"))
  
  length(unlist(strsplit(colony_data_lines[1], " ")))
  ((length(unlist(strsplit(colony_data_lines[1], " ")))) - 1)/2
}
colony_data(data_genind_adult)  
# adjust marker errors ----------------------------------------------------
data1 <- read.table(file = "./data/marker_info_data_genind_adult.txt", header = TRUE, dec = ",", na.strings = c("", ".", "na"))
data2 <- read.table(file = "C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Palau genetics mixing/Colony/palau14_full_run/palau14_full_run.ErrorRate", 
                    header = TRUE, sep = ",", na.strings = c("", ".", "na"))
error_rates <- data2 %>% select(MarkerID, DropRateEst)
data1_t <- as.data.frame(t(data1))
data1_t <- tibble::rownames_to_column(data1_t, "MarkerID")
updated_data <- data1_t %>% left_join(error_rates, by = "MarkerID") %>% mutate(V2 = ifelse(!is.na(DropRateEst), DropRateEst, V2)) %>% 
  select(-DropRateEst)
data1_updated <- as.data.frame(t(updated_data))
colnames(data1_updated) <- NULL
rownames(data1_updated) <- NULL
write.table(data1_updated, file = './data/marker_info_data_genind_adult_corr.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

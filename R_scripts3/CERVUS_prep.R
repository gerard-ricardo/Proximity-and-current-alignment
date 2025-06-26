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
                           "a_hyac_map_letters_code1_docent1.csv"))
# offspring file ----------------------------------------------------------
df2 = data_gl_filtered$other$ind.metrics
df2$id = paste0('X', df2$id)
mismatches <- setdiff(colnames(data1), df2$id)
df2_lar = subset(df2, stage!= 'adults')
df2_adult = subset(df2, stage!= 'larvae')
df1a = df2_lar %>% dplyr::select(., id, genotype )
df2_adult$id
known_dam <- df2_adult %>%
  group_by(genotype) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(genotype, id) %>%
  rename(known_dam = id) %>% data.frame()
df1b = left_join(df1a, known_dam, by = "genotype")
nrow(df1b)
df1 = dplyr::select(df1b, c(id, known_dam))
df1 <- df1[complete.cases(df1), ]
len = length(unique(df2_adult$id))
cands <- rep(df2_adult$id, length(df1$known_dam)) %>% sort(.)
cands_df <- matrix(cands, nrow = length(df1$known_dam), ncol = len, byrow = FALSE) %>% data.frame()
colnames(cands_df) <- rep("candidate", len)
candidate_indices <- which(colnames(cands_df) == "candidate")
colnames(cands_df)[candidate_indices] <- paste0("candidate_", seq_along(candidate_indices))
offspring_df <- cbind(df1, cands_df)
offspring_df
write.csv(offspring_df, row.names = FALSE, 
          file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Palau genetics mixing/Cervus", 
                           "offspring_a_hyac.csv"))

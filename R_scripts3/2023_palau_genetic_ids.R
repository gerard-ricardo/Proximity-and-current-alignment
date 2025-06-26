library(dplyr)
# 1 Import data -----------------------------------------------------------
load("./Rdata/2023_palau_dart_pos.RData")
unique(data1$Genotype)
duplicate_indices <- duplicated(data1$Genotype) | duplicated(data1$Genotype, fromLast = TRUE)
duplicate_row_indices <- which(duplicate_indices)
print(duplicate_row_indices)
duplicate_rows <- data1[duplicate_indices, ]
duplicate_rows
which(data1$Genotype == '61')
## join the position ID to simple ID
load("./Rdata/2023_palau_IDs.RData")
data6 = left_join(data1, data5, by = 'Genotype')
data7 <-data6[data6$Comment =='larvae',]
table(data7$id)
data7[which(data7$id == '5_30'), ]
##########
data2 <-data1[data1$PlateID ==3,]
spokes <- 1:7
distances <- c(5, 3, 10, 30)
replicates <- 1:2
expected_combinations <- expand.grid(spokes = spokes, 
                                     distance = distances, 
                                     replicate = replicates)
expected_ids <- apply(expected_combinations, 1, function(x) paste(x, collapse = "_"))
actual_ids <- as.character(data2$Genotype)
presence_check <- expected_ids %in% actual_ids
(missing_ids <- expected_ids[!presence_check])
which(data2$Genotype == '4_5_1')

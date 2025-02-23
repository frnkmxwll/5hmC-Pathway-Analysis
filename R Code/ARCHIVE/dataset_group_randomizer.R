###dataset_group_randomizer.R###
# Purpose of this script is to take in a file that is ready for DESeq2 analysis
# and split it into training, validation, and holdout sets.

### INSTALL LIBRARIES
# Setup, uncomment follow
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("DESeq2")
# install.packages("dplyr")
# BiocManager::install("sva")
# install.packages("caTools")

library(caret)
library(dplyr)
library(limma)
library(DESeq2)
library(sva)

### CONFIGURATION
# set working directory
setwd("~/5hmC-Pathway-Analysis/")
counts_name <- "./Output/Raw Data Processing/PMneg_PMpos_06022024_hmr_nocombat/PMneg_PMpos_DESeq2_rawcounts.csv"
normcounts_name <- "./Output/Raw Data Processing/PMneg_PMpos_06022024_hmr_nocombat/PMneg_PMpos_DESeq2_normcounts.csv"
meta_name <- "./Output/Raw Data Processing/PMneg_PMpos_06022024_hmr_nocombat/PMneg_PMpos_DESeq2_conditions.csv"

file_version <- "hmr_nocombat_602020_06022024_CRCHGA_cbpp_seed126"
random_seed = 126

# settings
training_fraction = 60
validation_fraction = 20
holdout_fraction = 20
min_samples = 0.8
min_count = 5
create_sig_training_counts = FALSE
combat_norm = FALSE
save_partition_files = TRUE

# Read in data
counts_data <- read.csv(counts_name, row.names = 1)
normcounts_data <- read.csv(normcounts_name, row.names = 1)
meta <- read.csv(meta_name, row.names = 1)

# Drop the same KT ids from counts_data and normcounts_data
counts_data <- counts_data[, colnames(counts_data) %in% rownames(meta)]
normcounts_data <- normcounts_data[, colnames(normcounts_data) %in% rownames(meta)]

# Define two comparison groups
first_condition = distinct(meta, condition)[1, 1]
second_condition = distinct(meta, condition)[2, 1]

train_prop <- training_fraction / 100
validation_prop <- validation_fraction / 100
holdout_prop <- holdout_fraction / 100

# Combine factors for stratification and identify small groups
meta$stratum <- interaction(meta$condition, meta$primary_present, meta$batch)
group_sizes <- table(meta$stratum)
small_groups <- names(group_sizes[group_sizes <= 2])

# Merge small groups into a single stratum
if (length(small_groups) > 0) {
  meta$stratum <- ifelse(meta$stratum %in% small_groups, "SmallGroup", meta$stratum)
}

# Calculate desired number of samples for each set
total_samples <- nrow(meta)
desired_train_size <- round(total_samples * train_prop)
desired_validation_size <- round(total_samples * validation_prop)
desired_holdout_size <- round(total_samples * holdout_prop)

# Simplified Stratified Partitioning with Proportional Adjustment (All Sets)
set.seed(random_seed)

# Start with the training set and adjust as needed
train_indices <- createDataPartition(meta$stratum, times = 1, p = train_prop, list = FALSE)
excess_train_samples <- length(train_indices) - desired_train_size
if (excess_train_samples > 0) {
  samples_to_move <- sample(train_indices, excess_train_samples)
  train_indices <- setdiff(train_indices, samples_to_move)
}

# Partition the remaining data with adjusted proportions
remaining_indices <- setdiff(1:nrow(meta), train_indices)
validation_indices <- createDataPartition(meta$stratum[remaining_indices], times = 1,
                                          p = desired_validation_size / length(remaining_indices), list = FALSE)
validation_indices <- remaining_indices[validation_indices]
holdout_indices <- setdiff(remaining_indices, validation_indices)

# Adjust holdout if needed to maintain proportions
excess_validation_samples <- length(validation_indices) - desired_validation_size
if (excess_validation_samples > 0) {
  samples_to_move <- sample(validation_indices, excess_validation_samples)
  validation_indices <- setdiff(validation_indices, samples_to_move)
  holdout_indices <- c(holdout_indices, samples_to_move)
}

# Create Data Frames
training_meta <- meta[train_indices, ]
validation_meta <- meta[validation_indices, ]
holdout_meta <- meta[holdout_indices, ]

# Print Proportions (Optional)
total_samples <- nrow(meta)
cat("Training samples: ", nrow(training_meta), "(", round(nrow(training_meta) / total_samples * 100, 2), "%)\n")
cat("Validation samples: ", nrow(validation_meta), "(", round(nrow(validation_meta) / total_samples * 100, 2), "%)\n")
cat("Holdout samples: ", nrow(holdout_meta), "(", round(nrow(holdout_meta) / total_samples * 100, 2), "%)\n")

training_meta_ordered <- training_meta[order(row.names(training_meta)), ]
training_meta_ordered <- training_meta_ordered[order(training_meta_ordered$condition),]
training_meta_ordered <- training_meta_ordered[, -which(names(training_meta_ordered) == "stratum")]

validation_meta_ordered <- validation_meta[order(row.names(validation_meta)), ]
validation_meta_ordered <- validation_meta_ordered[order(validation_meta_ordered$condition),]
validation_meta_ordered <- validation_meta_ordered[, -which(names(validation_meta_ordered) == "stratum")]

holdout_meta_ordered <- holdout_meta[order(row.names(holdout_meta)), ]
holdout_meta_ordered <- holdout_meta_ordered[order(holdout_meta_ordered$condition),]
holdout_meta_ordered <- holdout_meta_ordered[, -which(names(holdout_meta_ordered) == "stratum")]

training_counts <- counts_data[, colnames(counts_data) %in% rownames(training_meta_ordered)]
validation_counts <- counts_data[, colnames(counts_data) %in% rownames(validation_meta_ordered)]
holdout_counts <- counts_data[, colnames(counts_data) %in% rownames(holdout_meta_ordered)]

### FILTER GENES BASED ON MINIMUM COUNTS IN SAMPLES (Before Normalization)
training_counts <- training_counts[rowSums(training_counts >= min_count) >= min_samples * ncol(training_counts), ]

# Apply DESeq2 normalization and VST transformation
apply_combat <- function(counts, meta,combat_norm) {
  # Assume 'counts' is your raw count matrix and 'meta' is your metadata dataframe
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ 1)
    
  # Apply DESeq2 median ratio normalization
  dds <- estimateSizeFactors(dds)
  norm_counts <- counts(dds, normalized = TRUE)
  if (combat_norm){
    # Apply ComBat-Seq batch correction while preserving the condition covariate
      batch_corrected_counts <- ComBat_seq(norm_counts, batch = meta$batch)#, group = meta$condition)
      return(batch_corrected_counts)
    }
    else{
      return(norm_counts)
    }
}

training_normcounts <- apply_combat(training_counts, training_meta_ordered,combat_norm)
validation_normcounts <- apply_combat(validation_counts, validation_meta_ordered,combat_norm)
holdout_normcounts <- apply_combat(holdout_counts, holdout_meta_ordered,combat_norm)

if(save_partition_files){
  # Output resulting files
  folder <- "./Output/Randomization/"
  if (!file.exists(folder)) {
    dir.create("./Output/Randomization/")
  }
  
  output_folder <- paste("./Output/Randomization/", first_condition, "_", second_condition, "_", "DESeq2_", file_version, sep = "")
  dir.create(output_folder)
  
  # Output raw counts
  write.csv(training_meta_ordered, file = paste(output_folder, "/", first_condition, "_", second_condition, "_training_conditions.csv", sep = ""), row.names = TRUE, quote = FALSE)
  write.csv(validation_meta_ordered, file = paste(output_folder, "/", first_condition, "_", second_condition, "_validation_conditions.csv", sep = ""), row.names = TRUE, quote = FALSE)
  write.csv(holdout_meta_ordered, file = paste(output_folder, "/", first_condition, "_", second_condition, "_holdout_conditions.csv", sep = ""), row.names = TRUE, quote = FALSE)
  
  write.csv(training_counts, file = paste(output_folder, "/", first_condition, "_", second_condition, "_training_rawcounts.csv", sep = ""), row.names = TRUE, quote = FALSE)
  write.csv(validation_counts, file = paste(output_folder, "/", first_condition, "_", second_condition, "_validation_rawcounts.csv", sep = ""), row.names = TRUE, quote = FALSE)
  write.csv(holdout_counts, file = paste(output_folder, "/", first_condition, "_", second_condition, "_holdout_rawcounts.csv", sep = ""), row.names = TRUE, quote = FALSE)
  
  # Output VST-transformed counts
  write.csv(training_normcounts, file = paste(output_folder, "/", first_condition, "_", second_condition, "_training_normcounts.csv", sep = ""), row.names = TRUE, quote = FALSE)
  write.csv(validation_normcounts, file = paste(output_folder, "/", first_condition, "_", second_condition, "_validation_normcounts.csv", sep = ""), row.names = TRUE, quote = FALSE)
  write.csv(holdout_normcounts, file = paste(output_folder, "/", first_condition, "_", second_condition, "_holdout_normcounts.csv", sep = ""), row.names = TRUE, quote = FALSE)
}

### Save Config File
config <- c(
  paste("input counts file name:", counts_name), 
  paste("input normcounts file name:", normcounts_name), 
  paste("input conditions file name:", meta_name), 
  paste("file version:", file_version),
  paste("Training fraction:", training_fraction, "%"),
  paste("Validation fraction:", validation_fraction, "%"),
  paste("Holdout fraction:", holdout_fraction, "%"),
  paste("Output folder:", "./Output/Randomization/", first_condition, "_", second_condition, "_", "DESeq2_", file_version, sep = ""),
  paste("Random seed:", random_seed)
)

write.table(config, file = paste("./Output/Randomization/", first_condition, "_", second_condition, "_", "DESeq2_", file_version, "/", first_condition, "_", second_condition, "_config", ".txt", sep = ""), sep = "\t", quote = FALSE, col.names = NA)

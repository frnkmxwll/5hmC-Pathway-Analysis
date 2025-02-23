### Define Functions

### Define color tables
create_color_tables <- function() {
  # # Define contrast groups
  # contrast_groups <- c("condition", as.character(conditions$first_condition), as.character(conditions$second_condition))
  # group1 <- paste(contrast_groups[[3]])
  # group2 <- paste(contrast_groups[[2]])
  # 
  # # Check for an optional fourth group
  # if (length(contrast_groups) > 3) {
  #   group3 <- paste(contrast_groups[[4]])
  # }
  # 
  # # Define condition_table based on contrast groups
  # if (length(contrast_groups) > 3) {
  #   condition_table <- c("#1F78B4B6", "#E31A1CB6", "#CCCCCCA6")
  #   names(condition_table) <- c(group1, group2, group3)
  # } else {
  #   condition_table <- c("#1F78B4B6", "#E31A1CB6")
  #   names(condition_table) <- c(group1, group2)
  # }
  
  # Define other color tables using colorblind-safe palettes
  condition_colors <- c("#1F78B4B6", "#eb5c5d", "#eb5c5d","#CCCCCC","#eb5c5d","#1F78B4B6")
  condition_table <- setNames(condition_colors, c("PM_negative", "PM_positive", "DISEASED", "HEALTHY","AMpos","AMneg"))
  
  # Define other color tables using colorblind-safe palettes
  pm_mets_colors <- c("#1F78B4B6", "#eb5c5d", "#CCCCCC")
  pm_mets_table <- setNames(pm_mets_colors, c("pm_absent", "pm_present", "healthy"))
  
  # Define other color tables using colorblind-safe palettes
  pm_mets_colors_am <- c( "#D6D133","#F0EEA8", "#CCCCCC")
  pm_mets_table_am <- setNames(pm_mets_colors_am, c("pm_absent", "pm_present", "healthy"))
  
  primary_site_colors <- c("#009E73A6", "lightgreen", "lightgreen", "#CCCCCC","#009E73A6")
  primary_site_table <- setNames(primary_site_colors, c("CRC", "ADE", "appendiceal_mucocele", "healthy","HGA"))
  
  path_colors <- c("#009E73A6", "lightgreen", "#CCCCCC")
  path_table <- setNames(path_colors, c("AMN", "Adenocarcinoma", "unknown_or_na"))
  
  primary_present_colors <- c("#FAD2BB", "#FFAE7F", "#CCCCCC")
  primary_present_table <- setNames(primary_present_colors, c("primary_absent", "primary_present", "healthy"))
  
  sex_colors <- c("lightpink", "#CC79A7")
  sex_table <- setNames(sex_colors, c("male", "female"))
  
  other_mets_colors <- c("#E69F00D1", "#009E73A6", "#CCCCCC")
  other_mets_table <- setNames(other_mets_colors, c("other_mets_present", "other_mets_absent", "healthy"))
  
  chemo_colors <- c("#999999A1", "#999999F1")
  chemo_table <- setNames(chemo_colors, c("No", "Yes"))
  
  batch_colors <- c("#a98de8", "#daccee")
  batch_table <- setNames(batch_colors, c("batch_1", "batch_2"))
  
  age_cat_colors <- c("#F5B5AD", "#EF7F70")
  age_cat_table <- setNames(age_cat_colors, c("low", "high"))
  
  race_colors <- c("#9ED4F2", "#56B4E9", "#CDBDF1","#CDBDF1","#a98de8","#a98de8","#9ED4F2", "#56B4E9", "#CDBDF1", "#a98de8", "#d98de8c1", "#d98de8c1")
  race_table <- setNames(race_colors, c("White", "Black_African_American", "American_Indian_Alaska_Native","Asian_Mideast_Indian", "Asian_or_Pacific_Islander","Unknown_Other","zero","one","two","three","four","Other"))
  
  race_cat_colors <- c("#9ED4F2","#CDBDF1")
  race_cat_table <- setNames(race_colors, c("White","Other"))
  
  # **New: Define color table for locationType using muted colors**
  location_type_colors <- c( "#FFAE7F", "#FAD2BB", "#CDBDF1","#56B4E9","#009E73A6", "lightgreen", "#F0EEA8")
  location_type_table <- setNames(location_type_colors, c("1-5kb Upstream", "Promoter", "5' UTR", "Exon", "Intron", "3' UTR", "Intergenic"))
  
  
  # (Using 5 distinct colors for 150, 2000, 5000, 10000, 20000)
  num_top_genes_colors <- c("#eb5c5d","#EF7F70","#E69F00D1","#FFAE7F","#FAD2BB","#D6D133","lightgreen","#009E73A6", "#56B4E9","#1F78B4", "#CDBDF1","#a98de8",  "#d98de8c1","#CC79A7","#999999F1")
  num_top_genes_table <- setNames(num_top_genes_colors, 
                                  c("5","10","15","20","25","50","100","150","250","500","1000","2000","5000","10000","20000"))
  
  # Return a combined list of all color tables, including the new one:
  color_tables <- list(
    condition_table = condition_table,
    pm_mets_table = pm_mets_table,
    pm_mets_table_am = pm_mets_table_am,
    primary_site_table = primary_site_table,
    path_table = path_table,
    primary_present_table = primary_present_table,
    sex_table = sex_table,
    other_mets_table = other_mets_table,
    chemo_table = chemo_table,
    batch_table = batch_table,
    age_cat_table = age_cat_table,
    race_table = race_table,
    race_cat_table = race_cat_table,
    location_type_table = location_type_table,
    # New:
    num_top_genes_table = num_top_genes_table
  )
  
  
  return(color_tables)
}

## Partition Functions
load_and_prepare_data <- function(counts_name, normcounts_name, meta_name, debug_data_summary, debug_progress, rapid, num_raw_genes, 
                                  subset_conditions = NULL, 
                                  excluded_samples = NULL,
                                  random_conds = FALSE, collapse_healthy = FALSE) {
  # Load data
  counts_data <- read.csv(counts_name, row.names = 1)
  normcounts_data <- read.csv(normcounts_name, row.names = 1)
  meta <- read.csv(meta_name, row.names = 1)
  
  # if age_cat column does not exist
  if(!"age_cat" %in% colnames(meta)){
    # add "age_cat" column where >= median is high and < median is low
    meta$age_cat <- ifelse(meta$age >= median(meta$age, na.rm = TRUE), "high", "low")
  }
  
  # if race_cat column does not exist
  if(!"race_cat" %in% colnames(meta)){
    # add race_cat column
    meta$race_cat <- ifelse(meta$race == "White", "White", "Other")
  }
  
  # Debug summary before subsetting
  if (debug_progress) {
    cat("Number of samples before subsetting:\n")
    cat("Meta:", nrow(meta), "\n")
    cat("Counts data:", ncol(counts_data), "\n")
    cat("Normalized counts data:", ncol(normcounts_data), "\n")
  }
  
  # if rapid, take top 1000 rows of counts_data
  if (rapid) {
    counts_data <- counts_data[1:num_raw_genes,]
    normcounts_data <- normcounts_data[1:num_raw_genes,]
  }
  
  # Convert character columns in meta dataframe to factors
  meta[] <- lapply(meta, function(x) if (is.character(x)) factor(x) else x)
  
  # Subset meta based on the conditions provided
  if (!is.null(subset_conditions)) {
    for (condition in subset_conditions) {
      column_name <- condition[[1]]
      value <- unlist(condition[[2]])
      if (column_name %in% colnames(meta)) {
        # keep only rows where column value in column_name are in value
        meta <- meta[meta[, column_name] %in% value, , drop = FALSE]
      } else {
        warning(paste("Column", column_name, "not found in meta data. Skipping this condition."))
      }
    }
  }
  
  # if collapse_healthy == TRUE, collapse healthy samples into one group and all others into DISEASED
  if (collapse_healthy) {
    meta$condition <- ifelse(meta$condition == "HEALTHY", "HEALTHY", "DISEASED")
  }
  
  # replace forbidden characters "/" or "-" with "_" in values
  meta$race <- gsub("/", "_", meta$race)
  meta$race <- gsub("-", "_", meta$race)
  
  # if excluded_samples is not NULL
  if (!is.null(excluded_samples)) {
    meta <- meta[!rownames(meta) %in% excluded_samples, , drop = FALSE]
  }
  
  # if random_conds is TRUE
  if (random_conds) {
    # reassign meta$condition by randomizing order
    meta$condition <- sample(meta$condition)
    meta$peritoneal_mets <- sample(meta$peritoneal_mets)
    meta$primary_site <- sample(meta$primary_site)
  }
  # reorder with respect to condition then rownames
  meta <- meta[order(meta$condition, rownames(meta)), ]
  
  # Prepare data
  counts_data <- counts_data[, colnames(counts_data) %in% rownames(meta)]
  normcounts_data <- normcounts_data[, colnames(normcounts_data) %in% rownames(meta)]
  
  # reorder counts_data and normcounts_data columns order to match rownames order of meta
  counts_data <- counts_data[, rownames(meta)]
  normcounts_data <- normcounts_data[, rownames(meta)]
  
  
  # Debug summary after preparing data
  if (debug_progress) {
    cat("Number of samples after preparing data:\n")
    cat("Meta:", nrow(meta), "\n")
    cat("Counts data:", ncol(counts_data), "\n")
    cat("Normalized counts data:", ncol(normcounts_data), "\n")
  }
  
  if (debug_progress) cat("Data prepared for partitioning.\n")
  
  return(list(counts_data = counts_data, normcounts_data = normcounts_data, meta = meta))
}

define_conditions <- function(meta) {
  first_condition <- distinct(meta, condition)[1, 1]
  second_condition <- distinct(meta, condition)[2, 1]
  return(list(first_condition = first_condition, second_condition = second_condition))
}

partition_data <- function(meta, interaction_params, training_fraction, validation_fraction, holdout_fraction, random_seed, debug_progress) {
  train_prop <- training_fraction / 100
  validation_prop <- validation_fraction / 100
  holdout_prop <- holdout_fraction / 100
  
  meta$stratum <- do.call(interaction, meta[interaction_params])
  group_sizes <- table(meta$stratum)
  small_groups <- names(group_sizes[group_sizes <= 1])
  
  if (length(small_groups) > 0) {
    meta$stratum <- ifelse(meta$stratum %in% small_groups, "SmallGroup", meta$stratum)
  }
  
  total_samples <- nrow(meta)
  desired_train_size <- round(total_samples * train_prop)
  desired_validation_size <- round(total_samples * validation_prop)
  desired_holdout_size <- round(total_samples * holdout_prop)
  
  set.seed(random_seed)
  
  train_indices <- createDataPartition(meta$stratum, times = 1, p = train_prop, list = FALSE)
  excess_train_samples <- length(train_indices) - desired_train_size
  if (excess_train_samples > 0) {
    samples_to_move <- sample(train_indices, excess_train_samples)
    train_indices <- setdiff(train_indices, samples_to_move)
  }
  
  remaining_indices <- setdiff(1:nrow(meta), train_indices)
  validation_indices <- if (validation_fraction == 0) integer(0) else createDataPartition(meta$stratum[remaining_indices], times = 1, p = desired_validation_size / length(remaining_indices), list = FALSE)
  validation_indices <- remaining_indices[validation_indices]
  holdout_indices <- setdiff(remaining_indices, validation_indices)
  
  excess_validation_samples <- length(validation_indices) - desired_validation_size
  if (excess_validation_samples > 0) {
    samples_to_move <- sample(validation_indices, excess_validation_samples)
    validation_indices <- setdiff(validation_indices, samples_to_move)
    holdout_indices <- c(holdout_indices, samples_to_move)
  }
  
  training_meta <- meta[train_indices, ]
  validation_meta <- if (length(validation_indices) > 0) meta[validation_indices, ] else data.frame()
  holdout_meta <- meta[holdout_indices, ]
  
  if (debug_progress) {
    cat("Data partitioned:\n")
    cat("Training samples: ", nrow(training_meta), "\n")
    cat("Validation samples: ", nrow(validation_meta), "\n")
    cat("Holdout samples: ", nrow(holdout_meta), "\n")
  }
  
  return(list(training_meta = training_meta, validation_meta = validation_meta, holdout_meta = holdout_meta))
}

order_and_filter_data <- function(training_meta, validation_meta, holdout_meta, counts_data, median_percentile_cutoff, variance_percentile_cutoff, debug_progress) {
  order_and_filter <- function(meta, counts_data) {
    if (nrow(meta) == 0) return(list(meta_ordered = meta, counts = data.frame()))
    meta_ordered <- meta[order(row.names(meta)), ]
    meta_ordered <- meta_ordered[order(meta_ordered$condition), ]
    counts <- counts_data[, colnames(counts_data) %in% rownames(meta_ordered)]
    # order counts columns in the same order as rows in meta_ordered
    counts <- counts[, rownames(meta_ordered)]
    return(list(meta_ordered = meta_ordered, counts = counts))
  }
  
  normalize_counts <- function(counts) {
    total_counts <- colSums(counts)
    counts_normalized <- sweep(counts, 2, total_counts, FUN = "/")
    return(counts_normalized)
  }
  
  drop_low_median_samples <- function(counts, min_count_percentile, debug_progress) {
    if (debug_progress) {
      cat("Rows before median filtering:", nrow(counts), "\n")
    }
    
    counts_normalized <- normalize_counts(counts)
    row_medians <- apply(counts_normalized, 1, median)
    threshold <- quantile(row_medians, min_count_percentile)
    counts <- counts[row_medians > threshold, ]
    
    if (debug_progress) {
      cat("Rows after median filtering:", nrow(counts), "\n")
    }
    
    return(counts)
  }
  
  drop_low_variance_genes <- function(counts, variance_percentile, debug_progress) {
    if (debug_progress) {
      cat("Rows before variance filtering:", nrow(counts), "\n")
    }
    
    gene_means <- rowMeans(counts)
    gene_centered <- sweep(counts, 1, gene_means)
    gene_variances <- rowSums(gene_centered^2) / (ncol(counts) - 1)
    threshold <- quantile(gene_variances, variance_percentile)
    counts <- counts[gene_variances > threshold, ]
    
    if (debug_progress) {
      cat("Rows after variance filtering:", nrow(counts), "\n")
    }
    
    return(counts)
  }
  
  training <- order_and_filter(training_meta, counts_data)
  training$counts <- drop_low_median_samples(training$counts, median_percentile_cutoff, debug_progress)
  training$counts <- drop_low_variance_genes(training$counts, variance_percentile_cutoff, debug_progress)
  
  # if validation_meta not NULL
  if (!is.null(validation_meta)){
    validation <- order_and_filter(validation_meta, counts_data)
  } else {
    validation <- list(meta_ordered = data.frame(), counts = data.frame())
  }
  if (!is.null(holdout_meta)) {
    holdout <- order_and_filter(holdout_meta, counts_data)
  } else {
    holdout <- list(meta_ordered = data.frame(), counts = data.frame())
  }
  
  if (debug_progress) {
    cat("Data ordered and filtered.\n")
  }
  
  return(list(training = training, validation = validation, holdout = holdout))
}

normalize_counts <- function(training, validation, holdout = NULL, combat_norm = FALSE, length_norm = FALSE, debug_progress = FALSE) {

  normalize_dataset <- function(counts, meta) {
    if (nrow(counts) == 0 || ncol(counts) == 0) {
      return(matrix(nrow = 0, ncol = 0))
    }
    
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ 1)
    dds <- estimateSizeFactors(dds)
    norm_counts <- counts(dds, normalized = TRUE)
    
    if (debug_progress) {
      cat("DESeq2 normalization applied to account for library size differences.\n")
    }
    
    if (length_norm) {
      feature_lengths <- as.numeric(sapply(strsplit(rownames(norm_counts), "_"), function(x) as.numeric(x[3]) - as.numeric(x[2])))
      norm_counts <- t(t(norm_counts) / feature_lengths)
      
      if (debug_progress) {
        cat("Length normalization applied.\n")
      }
    }
    
    if (combat_norm) {
      norm_counts <- ComBat_seq(norm_counts, batch = meta$batch)
      
      if (debug_progress) {
        cat("ComBat-seq batch correction applied.\n")
      }
    }
    
    return(norm_counts)
  }

  training_normcounts <- normalize_dataset(training$counts, training$meta_ordered)
  
  if(is.null(validation) || nrow(validation$counts) == 0){
    validation_normcounts <- matrix(nrow = 0, ncol = 0)
  } else {
    validation_normcounts <- normalize_dataset(validation$counts, validation$meta_ordered)
  }
  if (is.null(holdout) || nrow(holdout$counts) == 0) {
    holdout_normcounts <- matrix(nrow = 0, ncol = 0)
  } else {
    holdout_normcounts <- normalize_dataset(holdout$counts, holdout$meta_ordered)
  }
  
  if (debug_progress) {
    cat("Normalization and batch correction applied to all datasets.\n")
  }
  
  return(list(training_normcounts = training_normcounts, 
              validation_normcounts = validation_normcounts, 
              holdout_normcounts = holdout_normcounts))
}

## DESeq2 Functions
perform_deseq2 <- function(training_data, design_formula, cutoff_type, padj_cutoff, pvalue_cutoff, lfc_cutoff, num_top_genes, debug_progress, debug_data_summary) {
  # Register settings for DESeq2 parallelization

  
  if(Sys.info()[['sysname']]=="Linux"){
    BiocParallel::register(MulticoreParam(workers=20), default = TRUE)
    # Run DESeq2
    # Create DESeq2 object
    dds <- DESeqDataSetFromMatrix(countData = training_data$counts,
                                  colData = training_data$meta_ordered,
                                  design = design_formula)
    dds <- DESeq2::DESeq(dds, parallel=TRUE)
  }
  else{
    # Run DESeq2
    # Create DESeq2 object
    dds <- DESeqDataSetFromMatrix(countData = training_data$counts,
                                  colData = training_data$meta_ordered,
                                  design = design_formula)
    dds <- DESeq2::DESeq(dds)#, parallel=FALSE)
  }

  # Debug statement to display DESeq2 summary if debug_data_summary is TRUE
  if (debug_data_summary) {
    cat("DESeq2 Summary:\n")
    print(summary(dds))
  }
  
  # Define contrast groups based on metadata (within the function)
  groups <- unique(training_data$meta_ordered["condition"])
  contrast_groups <- c("condition", as.character(groups[1, 1]), as.character(groups[2, 1]))
  
  # Get results with specified contrast
  res <- results(dds, contrast = contrast_groups)
  
  # Convert to tibble and add gene names
  res_table <- as.data.frame(res)
  res_table$gene <- rownames(res_table)
  
  # Filter significant genes based on cutoff type
  if (cutoff_type == 0) {
    sig <- res_table %>% filter(padj <= padj_cutoff, abs(log2FoldChange) >= lfc_cutoff)
  } else if (cutoff_type == 1) {
    sig <- res_table %>% filter(abs(log2FoldChange) >= lfc_cutoff, pvalue <= pvalue_cutoff)
  } else if (cutoff_type == 2) {
    # Sort by stat in descending order
    res_table_sorted <- res_table %>% arrange(desc(abs(stat)))
    # Filter res_table_tb down to num_top_genes first rows only
    sig <- res_table_sorted %>% slice_head(n = num_top_genes)
  }
  
  if (debug_data_summary) {
    cat("DESeq2 significant genes:\n")
    print(dim(sig))
  }
  
  return(sig)
}

deseq2_filter <- function(training_data, design_formula, cutoff_type, padj_cutoff, pvalue_cutoff, lfc_cutoff, num_top_genes,x,y) {
  # Create DESeq2 object
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = training_data$counts,
                                        colData = training_data$meta_ordered,
                                        design = design_formula)
  if(Sys.info()[['sysname']]=="Linux"){
    BiocParallel::register(MulticoreParam(workers=20), default = TRUE)
    # Run DESeq2
    dds <- DESeq2::DESeq(dds, quiet=TRUE, parallel=TRUE)
  }
  else{
    # Run DESeq2
    dds <- DESeq2::DESeq(dds, quiet=TRUE)#, parallel=FALSE)
  }
  
  # Define contrast groups based on metadata (within the function)
  groups <- unique(training_data$meta_ordered["condition"])
  contrast_groups <- c("condition", as.character(groups[1, 1]), as.character(groups[2, 1]))
  
  # Get results with specified contrast
  res <- DESeq2::results(dds, contrast = contrast_groups)
  
  # Convert to tibble and add gene names
  res_table <- as.data.frame(res)
  res_table$gene <- rownames(res_table)
  
  # Filter significant genes based on cutoff type
  if (cutoff_type == 0) {
    sig <- res_table[res_table$padj <= padj_cutoff & abs(res_table$log2FoldChange) >= lfc_cutoff, ]
  } else if (cutoff_type == 1) {
    sig <- res_table[abs(res_table$log2FoldChange) >= lfc_cutoff & res_table$pvalue <= pvalue_cutoff, ]
  } else if (cutoff_type == 2) {
    # Sort by stat in descending order
    res_table_sorted <- res_table[order(-abs(res_table$stat)), ]
    # Filter res_table_tb down to num_top_genes first rows only
    sig <- res_table_sorted[1:num_top_genes, ]
  }
  
  sig_indices <- match(sig$gene, colnames(t(training_data$counts)))
  
  return(sig_indices)
}

deseq2_nestglm_filter <- function(training_data, design_formula, cutoff_type, padj_cutoff, pvalue_cutoff, lfc_cutoff, num_top_genes,
                                  norm_counts_data_training, conditions_vector_training, alpha_values, number_of_runs, nfolds_val, 
                                  standardize_bit, class_weights, meta_training, feat_occ_thresholds, random_seed, n_outer_folds, n_inner_folds,
                                  perform_deseq2_only,specific_folds_inner,x_local,y_local) {
  
  ### Subsetting to only the local data
  y_local = y_local
  x_local = x_local
  training_meta = training_data$meta_ordered
  raw_training_counts = training_data$counts
  norm_counts_data_training = norm_counts_data_training
  conditions_vector_training = conditions_vector_training
  class_weights = class_weights
  
  # Find the column names in raw_training_counts that are also in x_local_rows
  x_local_rows <- rownames(x_local)
  common_columns <- intersect(colnames(raw_training_counts), x_local_rows)
  # find indices of common columns in raw_training_counts
  x_local_indices <- match(common_columns, colnames(raw_training_counts))
  
  
  raw_training_counts = raw_training_counts[, common_columns, drop = FALSE]
  training_meta = training_meta[common_columns,,drop = FALSE]
  norm_counts_data_training = norm_counts_data_training[, common_columns, drop = FALSE]
  conditions_vector_training = conditions_vector_training[x_local_indices]
  class_weights = class_weights[x_local_indices]

  # Filter training_counts to keep only columns that are in the x_local row indices
  
    # Debug logging
  log_debug <- function(message, data = NULL) {
    cat(message, "\n")
    flush.console()
    if (!is.null(data)) {
      if (is.vector(data) || is.factor(data)) {
        print(data)
      } else {
        str(data)
      }
    }
  }
  optimal_auc_tra_outter <- 0
  selected_predictors <- NULL
  best_filter_alpha <<- NULL
  best_filter_threshold <<- NULL
  
  log_debug("Starting deseq2_nestglm_filter")
  
  # Create DESeq2 object
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = raw_training_counts, colData = training_meta, design = design_formula)
  log_debug("DESeq2 object created")
  
  if (Sys.info()[['sysname']] == "Linux") {
    BiocParallel::register(MulticoreParam(workers = 2), default = TRUE)
    dds <- DESeq2::DESeq(dds, quiet = TRUE, parallel = TRUE)
  } else {
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
  }
  log_debug("DESeq2 analysis completed")
  
  groups <- unique(training_data$meta_ordered["condition"])
  contrast_groups <- c("condition", as.character(groups[1, 1]), as.character(groups[2, 1]))
  res <- DESeq2::results(dds, contrast = contrast_groups)
  log_debug("DESeq2 results obtained")
  
  res_table <- as.data.frame(res)
  res_table$gene <- rownames(res_table)
  log_debug("DESeq2 results table created")
  
  if (cutoff_type == 0) {
    sig <- res_table[res_table$padj <= padj_cutoff & abs(res_table$log2FoldChange) >= lfc_cutoff, ]
  } else if (cutoff_type == 1) {
    sig <- res_table[abs(res_table$log2FoldChange) >= lfc_cutoff & res_table$pvalue <= pvalue_cutoff, ]
  } else if (cutoff_type == 2) {
    res_table_sorted <- res_table[order(-abs(res_table$stat)), ]
    sig <- res_table_sorted[1:num_top_genes, ]
  }
  log_debug("Significant genes filtered")
  # Match significant gene names with the row names of norm_counts_data_training
  sig_indices <- match(sig$gene, rownames(norm_counts_data_training))
  
  
  if(perform_deseq2_only==0){
    # Log the significant gene indices for debugging
    log_debug(paste("performing repeated glmnet on: ",length(sig_indices)," sig_indices"))
    
    # Ensure sig_indices do not contain NA values
    sig_indices <- sig_indices[!is.na(sig_indices)]
    # Subset norm_counts_data_training to only include the significant genes
    norm_counts_data_training_sig <- t(norm_counts_data_training[sig_indices, , drop = FALSE])
    conditions_vector_training_fac <- as.factor(make.names(conditions_vector_training))
    
    #alpha_results <- future.apply::future_lapply(alpha_values, function(alpha_value) {
    alpha_results <- lapply(alpha_values, function(alpha_value) {
    
      #compiled_results <- future.apply::future_lapply(1:number_of_runs, function(i) {
      compiled_results <- lapply(1:number_of_runs, function(i) {
        #set.seed(random_seed)
        #stratification_factor <- interaction(training_meta$stratum)
        #customIndex <- createFolds(stratification_factor, k = nfolds_val, list = TRUE)
        #foldid <- rep(NA, length(stratification_factor))
        # assign consistent foldids to make alpha performances directly comparable.
        
        #for (fold in seq_along(customIndex)) {
        #  foldid[customIndex[[fold]]] <- fold
        #}
        
        #registerDoParallel(2)
        cvfit <- cv.glmnet(
          x = norm_counts_data_training_sig,
          y = conditions_vector_training_fac,
          family = "binomial",
          type.measure = "auc",
          nfolds = nfolds_val,
          alpha = alpha_value,
          weights = class_weights,
          #foldid = foldid, #ensures same fold arrangement for all alphas
          standardize = standardize_bit,
          parallel=FALSE
        )
        tmp_coeffs <- coef(cvfit, s = "lambda.1se")
        
        tryCatch({
          new_results <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
        }, error = function(e) {
          log_debug("Error: Failed to extract Dimnames:", e$message)
          return(data.frame())
        })
        
        if (nrow(new_results) == 0) {
          log_debug(paste("No coefs for alpha:", alpha_value))
          return(data.frame())
        }
        return(data.frame(run_number = rep(i, nrow(new_results)), new_results))
      })#, future.packages = c("glmnet", "dplyr", "caret"), future.seed = TRUE)
      
      compiled_results <- do.call(rbind, compiled_results)
      set.seed(random_seed)
      
      if (is.null(compiled_results)) {
        log_debug("Skipping alpha...")
        return(data.frame())
      }
      coefficient_table <- compiled_results %>%
        filter(coefficient != 0) %>%
        group_by(name) %>%
        summarise(
          freq_coef = n(),
          mean_coef = mean(coefficient),
          min_coef = min(coefficient),
          max_coef = max(coefficient)
        )
      
      coefficient_table$abs_mean <- abs(coefficient_table$mean_coef)
      coefficient_table <- coefficient_table %>%
        arrange(desc(freq_coef), desc(abs_mean)) %>%
        filter(name != "(Intercept)")
      
      #thresh_results <- future.apply::future_lapply(feat_occ_thresholds, function(threshold) {
      thresh_results <- lapply(feat_occ_thresholds, function(threshold) {
      
        if (nrow(coefficient_table) <= 1) {
          log_debug(paste("Threshold", threshold, "has 1 or fewer rows. Skipping."))
          return(NULL)
        }
        
        selected_predictors <- coefficient_table %>%
          filter(freq_coef >= threshold * number_of_runs) %>%
          pull(name)
        
        if (is.null(selected_predictors)) {
          log_debug("Printing selected_predictors")
          print(selected_predictors)
          log_debug("Printing coefficient_table")
          print(head(coefficient_table))
          log_debug(paste("Threshold", threshold, "no or fewer than 2 features selected. Skipping."))
          return(NULL)
        }
        norm_counts_data_training <- t(norm_counts_data_training)
        norm_counts_data_training_top <- norm_counts_data_training[, selected_predictors, drop = FALSE]

        repeated_fit <- tryCatch({
          nestcv.glmnet(
            conditions_vector_training_fac,
            norm_counts_data_training_top,
            family = "binomial",
            type.measure = "auc",
            weights = class_weights,
            #alphaSet = alpha_values,
            cv.cores = 1,
            n_outer_folds = n_outer_folds, 
            n_inner_folds = n_inner_folds,
            finalCV = TRUE,
            min_1se = 1,
            outer_method = "cv",
            verbose = FALSE,
            filterFUN = NULL, 
            filter_options = NULL
          )
        }, error = function(e) {
          log_debug("Error in nestcv.glmnet:", e$message)
          return(NULL)
        })
        
        if (is.null(repeated_fit)) {
          return(NULL)
        }
        
        metrics_row <- data.frame(
          filter_alpha = alpha_value,
          filter_threshold = threshold,
          auc_tra_outter = signif(pROC::auc(repeated_fit$roc), 3)
        )
        return(list(metrics_row = metrics_row, repeated_fit = repeated_fit))
      })#, future.packages = c("glmnet", "dplyr", "caret", "ROCR"), future.seed = TRUE)
      return(list(alpha_value = alpha_value, thresh_results = thresh_results))
    })#, future.packages = c("glmnet", "dplyr", "caret", "ROCR"), future.seed = TRUE)
    
    for (alpha_result in alpha_results) {
      alpha_value <- alpha_result$alpha_value
      thresh_results <- alpha_result$thresh_results
      for (thresh_result in thresh_results) {
        if (is.null(thresh_result)) next
        if (thresh_result$metrics_row$auc_tra_outter > optimal_auc_tra_outter) {
          optimal_auc_tra_outter <- thresh_result$metrics_row$auc_tra_outter
          final_fit <- thresh_result$repeated_fit$final_fit
          alpha_lambda <- thresh_result$repeated_fit$final_param
          best_filter_alpha <<- thresh_result$metrics_row$filter_alpha
          best_filter_threshold <<- thresh_result$metrics_row$filter_threshold
          selected_predictors <- final_fit$glmnet.fit$beta@Dimnames[[1]]
          log_debug("Updated best model!")
        }
      }
    }
    
    if (is.null(selected_predictors)) {
      log_debug("Error: selected_predictors is NULL after processing")
      return(NULL)
    }
    selected_indices <- match(selected_predictors, rownames(norm_counts_data_training))
  }else{
    selected_indices <- sig_indices
  }
  
  log_debug("Sending predictors indices:", length(selected_indices))

  return(selected_indices)
}


### glmnet Functions
# ranger filter for nestcv:
# Function to create a random forest and get important variables
custom_ranger_filter <- function(x, y, num.trees = NULL, mtry = NULL, importance = "impurity", min.node.size = NULL,max.depth = NULL, class.weights = NULL, seed = 1, nfilter = 25) {
  # Check if mtry is NULL, set default value

  # Set seed for reproducibility
  set.seed(seed)
  
  # Fit the random forest model
  rf_model <- ranger(
    x = x,
    y = as.factor(y),
    num.trees = num.trees,
    max.depth = max.depth,
    min.node.size = min.node.size,
    mtry = mtry, 
    importance = importance,
    class.weights = class.weights,
    classification = TRUE,
    write.forest = FALSE,
    verbose = FALSE
    #num.threads = 1
  )
  
  # Get variable importance
  importance_values <- rf_model$variable.importance
  #cat("Variable importance values: ", importance_values, "\n")
  
  # Get the indexes of the most important variables
  important_vars <- order(importance_values, decreasing = TRUE)[1:nfilter]
  # debug output important vars
  #cat("Important variables: ", important_vars, "\n")
  return(important_vars)
}

apply_mapping <- function(meta, map) {
  meta[] <- lapply(meta, function(x) {
    x <- as.character(x)
    for (i in 1:length(map)) {
      x[x == names(map)[i]] <- map[[i]]
    }
    return(x)
  })
  return(meta)
}

train_glmnet_model <- function(norm_counts_data_training,conditions_vector_training,alpha_values, number_of_runs, nfolds_val, standardize_bit,class_weights, meta_training,
                               debug_progress,norm_counts_data_validation,conditions_vector_validation,conditions_vector_training_fac,conditions_vector_validation_fac,
                               feat_occ_thresholds, num_top_genes,rep_val,root_folder) {
  optimal_auc_val <- -Inf
  best_model <- NULL
  optimal_alpha <- NA
  optimal_threshold <- NA
  best_model_metrics <- NULL
  roc_data_tra_best <- list()
  roc_data_val_best <- list()
  model_metrics <- data.frame()
  
  with_progress({
    handlers("txtprogressbar")
    runs_prog <- progressr::progressor(number_of_runs)
    alpha_results <- future.apply::future_lapply(alpha_values, function(alpha_value) {
      compiled_results <- future.apply::future_lapply(1:number_of_runs, function(i) {
        runs_prog(message = sprintf("Processing run %g", i))
        set.seed(i)
        stratification_factor <- interaction(meta_training$stratum)
        customIndex <- createFolds(stratification_factor, k = nfolds_val, list = TRUE)
        foldid <- rep(NA, length(stratification_factor))
        for (fold in seq_along(customIndex)) {
          foldid[customIndex[[fold]]] <- fold
        }
        cvfit <- cv.glmnet(
          x = norm_counts_data_training,
          y = conditions_vector_training,
          family = "binomial",
          type.measure = "auc",
          nfolds = nfolds_val,
          alpha = alpha_value,
          weights = class_weights,
          foldid = foldid,
          standardize = standardize_bit
        )
        tmp_coeffs <- coef(cvfit, s = "lambda.min")
        new_results <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
        data.frame(run_number = rep(i, nrow(new_results)), new_results)
      }, future.packages = c("glmnet", "dplyr", "caret", "progress"), future.seed = TRUE)
      
      compiled_results <- do.call(rbind, compiled_results)
      
      set.seed(random_seed) #revert from incrementing seeds in above loop.
      
      coefficient_table <- compiled_results %>%
        filter(coefficient != 0) %>%
        group_by(name) %>%
        summarise(
          freq_coef = n(),
          mean_coef = mean(coefficient),
          min_coef = min(coefficient),
          max_coef = max(coefficient)
        )
      
      coefficient_table$abs_mean <- abs(coefficient_table$mean_coef)
      coefficient_table <- coefficient_table %>%
        arrange(desc(freq_coef), desc(abs_mean)) %>%
        filter(name != "(Intercept)")
      
      #if (debug_progress) {
      #  print("Coefficient table:")
      #  print(coefficient_table, n = 5)
      #}
      
      # Initialize roc_data_tra and roc_data_val before the loop
      roc_data_tra <- list()
      roc_data_val <- list()
      
      # Parallelize the threshold loop
      thresh_results <- future.apply::future_lapply(feat_occ_thresholds, function(threshold) {
        
        selected_predictors <- coefficient_table %>%
          filter(freq_coef >= threshold * number_of_runs) %>%
          pull(name)
        
        # Skip iteration if no features are above the threshold
        if (length(selected_predictors) < 2) {
          if (debug_progress) {
            cat("No features selected for threshold:", threshold, "\n")
          }
          return(NULL)
        }
        
        #if (debug_progress) {
        #  print(paste("Threshold:", threshold))
        #  print("Selected predictors:")
        #  print(selected_predictors)
        #}
        
        norm_counts_data_training_top <- norm_counts_data_training[, selected_predictors, drop = FALSE]
        norm_counts_data_validation_top <- norm_counts_data_validation[, selected_predictors, drop = FALSE]
        
        stratification_factor <- interaction(meta_training$stratum)
        customIndex <- createFolds(stratification_factor, k = nfolds_val, returnTrain = TRUE)
        
        train_control <- trainControl(
          method = "repeatedcv",
          number = nfolds_val,
          repeats = rep_val,
          summaryFunction = twoClassSummary,
          classProbs = TRUE,
          index = customIndex,
          verboseIter = FALSE
        )
        
        glmnet_model <- train(
          x = norm_counts_data_training_top,
          y = conditions_vector_training_fac,
          method = "glmnet",
          trControl = train_control,
          metric = "ROC",
          weights = class_weights,
          standardize = standardize_bit
        )
        
        probabilities_val <- predict(glmnet_model, newdata = as.matrix(norm_counts_data_validation_top), type = "prob")
        probabilities_tra <- predict(glmnet_model, newdata = as.matrix(norm_counts_data_training_top), type = "prob")
        
        pred_val <- prediction(probabilities_val[, 2], conditions_vector_validation)
        pred_tra <- prediction(probabilities_tra[, 2], conditions_vector_training)
        
        roc_perf_val <- performance(pred_val, "tpr", "fpr")
        roc_perf_tra <- performance(pred_tra, "tpr", "fpr")
        
        auc_tra <- as.numeric(performance(pred_tra, measure = "auc")@y.values)
        auc_val <- as.numeric(performance(pred_val, measure = "auc")@y.values)
        
        cm_val <- confusionMatrix(factor(ifelse(probabilities_val[, 2] > 0.5, "X1", "X0")), conditions_vector_validation_fac)
        cm_tra <- confusionMatrix(factor(ifelse(probabilities_tra[, 2] > 0.5, "X1", "X0")), conditions_vector_training_fac)
        
        #if(debug_progress) {
        #  print("Model trained with following performance:")
        #  print("Training AUC:")
        #  print(auc_tra)
        #  print("Validation AUC:")
        #  print(auc_val)
        #}
        
        metrics_row <- data.frame(
          alpha = alpha_value,
          threshold = threshold,
          num_features = length(selected_predictors),
          auc_tra = auc_tra,
          sensitivity_tra = cm_tra$byClass["Sensitivity"],
          specificity_tra = cm_tra$byClass["Specificity"],
          recall_tra = cm_tra$byClass["Recall"],
          precision_tra = cm_tra$byClass["Precision"],
          accuracy_tra = cm_tra$overall["Accuracy"],
          f1_tra = (2*(cm_tra$byClass["Precision"] * cm_tra$byClass["Recall"])/(cm_tra$byClass["Recall"] + cm_tra$byClass["Precision"])),
          auc_val = auc_val,
          sensitivity_val = cm_val$byClass["Sensitivity"],
          specificity_val = cm_val$byClass["Specificity"],
          recall_val = cm_val$byClass["Recall"],
          precision_val = cm_val$byClass["Precision"],
          accuracy_val = cm_val$overall["Accuracy"],
          f1_val = (2*(cm_val$byClass["Precision"] * cm_val$byClass["Recall"])/(cm_val$byClass["Recall"] + cm_val$byClass["Precision"])),
          number_of_runs = number_of_runs,
          num_top_genes = num_top_genes,
          nfolds_val = nfolds_val,
          root_folder = root_folder
        )
        
        list(metrics_row = metrics_row, roc_perf_tra = roc_perf_tra, auc_tra = auc_tra, roc_perf_val = roc_perf_val, auc_val = auc_val, glmnet_model = glmnet_model)
        
      }, future.packages = c("glmnet", "dplyr", "caret", "ROCR"), future.seed = TRUE)
      
      return(list(alpha_value = alpha_value, thresh_results = thresh_results, roc_data_tra = roc_data_tra, roc_data_val = roc_data_val))
    }, future.packages = c("glmnet", "dplyr", "caret", "progress", "ROCR"), future.seed = TRUE)
  })
  
  # Filter out NULL results
  thresh_results <- do.call(rbind, Filter(Negate(is.null), results))
  
  for (alpha_result in alpha_results) {
    alpha_value <- alpha_result$alpha_value
    thresh_results <- alpha_result$thresh_results
    roc_data_tra <- alpha_result$roc_data_tra
    roc_data_val <- alpha_result$roc_data_val
    
    for (result in thresh_results) {
      if (is.null(result)) next  # Skip null results
      
      metrics_row <- result$metrics_row
      
      #if (debug_progress) {
      #  print("metrics_row:")
      #  print(metrics_row)
      #}
      print("Result:")
      print(result)
      
      roc_perf_tra <- result$roc_perf_tra
      auc_tra <- result$auc_tra
      roc_perf_val <- result$roc_perf_val
      auc_val <- result$auc_val
      glmnet_model <- result$glmnet_model
      
      print("auc_tra")
      print(auc_tra)
      print("glmnet_model")
      print(glmnet_model)
      
      model_metrics <- rbind(model_metrics, metrics_row)
      
      #if (debug_progress) {
      #  print(paste("Updating roc_data_tra for threshold:", metrics_row$threshold))
      #}
      
      roc_data_tra[[as.character(metrics_row$threshold)]] <- list(roc_perf_tra, auc_tra)
      roc_data_val[[as.character(metrics_row$threshold)]] <- list(roc_perf_val, auc_val)
      
      if (auc_val > optimal_auc_val) {
        optimal_auc_val <- auc_val
        best_model <- glmnet_model
        optimal_alpha <- metrics_row$alpha
        optimal_threshold <- metrics_row$threshold
        best_model_metrics <- metrics_row
        roc_data_tra_best <- roc_data_tra
        roc_data_val_best <- roc_data_val
        best_params <- glmnet_model$bestTune
        best_coefficients <- as.matrix(coef(glmnet_model$finalModel, s = best_params$lambda))
      }
    }
  }
  
  return(list(
    best_model = best_model,
    optimal_alpha = optimal_alpha,
    optimal_threshold = optimal_threshold,
    best_model_metrics = best_model_metrics,
    best_params = best_params,
    best_coefficients = best_coefficients,
    roc_data_tra_best = roc_data_tra_best,
    roc_data_val_best = roc_data_val_best
  ))
}

train_nestedcv_model <- function(norm_counts_data_training,conditions_vector_training,alpha_values, number_of_runs, nfolds_val, standardize_bit,class_weights, meta_training,
                               debug_progress,norm_counts_data_validation,conditions_vector_validation,conditions_vector_training_fac,conditions_vector_validation_fac,
                               feat_occ_thresholds, num_top_genes,rep_val,root_folder,random_seed,n_outer_folds,n_inner_folds) {
  optimal_auc_tra_outter <- -Inf
  best_model <- NULL
  selected_predictors <- NULL
  optimal_alpha <- NA
  optimal_threshold <- NA
  best_model_metrics <- NULL
  roc_data_tra_best <- list()
  roc_data_val_best <- list()
  model_metrics <- data.frame()
  
  with_progress({
    handlers("txtprogressbar")
    runs_prog <- progressr::progressor(number_of_runs)
    alpha_results <- future.apply::future_lapply(alpha_values, function(alpha_value) {
      compiled_results <- future.apply::future_lapply(1:number_of_runs, function(i) {
        runs_prog(message = sprintf("Processing run %g", i))
        set.seed(i)
        stratification_factor <- interaction(meta_training$stratum)
        customIndex <- createFolds(stratification_factor, k = nfolds_val, list = TRUE)
        foldid <- rep(NA, length(stratification_factor))
        for (fold in seq_along(customIndex)) {
          foldid[customIndex[[fold]]] <- fold
        }
        cvfit <- cv.glmnet(
          x = norm_counts_data_training,
          y = conditions_vector_training,
          family = "binomial",
          type.measure = "auc",
          nfolds = nfolds_val,
          alpha = alpha_value,
          weights = class_weights,
          foldid = foldid,
          standardize = standardize_bit
        )
        tmp_coeffs <- coef(cvfit, s = "lambda.min")
        tryCatch({
          new_results <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
        }, error = function(e) {
          cat("Error: Failed to extract Dimnames - ", e$message, "\n")
          return(data.frame())
        })
        
        if (nrow(new_results) == 0) {
          cat(paste("No coefs for alpha: ",alpha_value))
          return(data.frame())  # Return an empty data frame if new_results is empty
        }
        return(data.frame(run_number = rep(i, nrow(new_results)), new_results))
      }, future.packages = c("glmnet", "dplyr", "caret", "progress"), future.seed = TRUE)
      
      compiled_results <- do.call(rbind, compiled_results)
      
      set.seed(random_seed) #revert from incrementing seeds in above loop.
      
      if (is.null(compiled_results)) {
        cat(paste("Skipping alpha..."))
        return(data.frame())  # Return an empty data frame if new_results is empty
      }
      coefficient_table <- compiled_results %>%
        filter(coefficient != 0) %>%
        group_by(name) %>%
        summarise(
          freq_coef = n(),
          mean_coef = mean(coefficient),
          min_coef = min(coefficient),
          max_coef = max(coefficient)
        )
      
      coefficient_table$abs_mean <- abs(coefficient_table$mean_coef)
      coefficient_table <- coefficient_table %>%
        arrange(desc(freq_coef), desc(abs_mean)) %>%
        filter(name != "(Intercept)")
      
      #if (debug_progress) {
      #  print("Coefficient table:")
      #  print(coefficient_table, n = 5)
      #}
      
      # Initialize roc_data_tra and roc_data_val before the loop
      roc_data_tra <- list()
      roc_data_val <- list()
      
      # Parallelize the threshold loop
      thresh_results <- future.apply::future_lapply(feat_occ_thresholds, function(threshold) {
        
        # Check if coefficient_table has 2 or fewer rows
        if (nrow(coefficient_table) <= 2) {
          if (debug_progress) {
            cat("Coefficient table has 2 or fewer rows. Skipping threshold:", threshold, "\n")
          }
          return(NULL)
        }
        
        selected_predictors <- coefficient_table %>%
          filter(freq_coef >= threshold * number_of_runs) %>%
          pull(name)
        
        # Skip iteration if no features are above the threshold
        if (is.null(selected_predictors) || length(selected_predictors) <2 ) {
          if (debug_progress) {
            cat("No or fewer than 2 features selected for threshold:", threshold, "\n")
          }
          return(NULL)
        }
        
        if (debug_progress) {
          print(paste("Starting on alpha:", alpha_value, " with threshold:",threshold ))
          print(paste("Selected model with ",length(selected_predictors), " predictors."))
        }
        
        norm_counts_data_training_top <- norm_counts_data_training[, selected_predictors, drop = FALSE]
        stratification_factor <- interaction(meta_training$stratum)
        customIndex <- caret::createFolds(stratification_factor, k = n_outer_folds)#, returnTrain = TRUE)
        
        #cat("norm_counts_data_training_typ: ",norm_counts_data_training_top)
        
        repeated_fit <- nestcv.glmnet(
          conditions_vector_training_fac,
          norm_counts_data_training_top,
          family = "binomial",
          weights = class_weights,
          alphaSet = alpha_values,
          cv.cores = 4,
          n_outer_folds = n_outer_folds, 
          n_inner_folds = n_inner_folds,
          min_1se = 0,
          finalCV = TRUE,
          outer_method = c("cv"),
          outer_folds = customIndex,
          #penalty.factor = rep(1,ncol(training_normcounts_matrix)),
          verbose = FALSE,
          filterFUN = NULL, 
          filter_options = NULL
        )
        metrics_row <- data.frame(
          filter_alpha = alpha_value,
          filter_threshold = threshold,
          num_features = length(repeated_fit$final_vars),
          auc_tra_inner = signif(pROC::auc(innercv_roc(repeated_fit)),3),
          auc_tra_outter = signif(pROC::auc(repeated_fit$roc),3), # Outter fold auc
          filter_runs = number_of_runs,
          filter_nfolds = nfolds_val,
          n_inner_folds = n_inner_folds,
          n_outer_folds = n_outer_folds
        )
        return(list(metrics_row = metrics_row, repeated_fit = repeated_fit))
      }, future.packages = c("glmnet", "dplyr", "caret", "ROCR"), future.seed = TRUE)
      return(list(alpha_value = alpha_value, thresh_results = thresh_results))
    }, future.packages = c("glmnet", "dplyr", "caret", "progress", "ROCR"), future.seed = TRUE)
  })
  
  # Filter out NULL results
  #thresh_results <- do.call(rbind, Filter(Negate(is.null), results))
  
  for (alpha_result in alpha_results) {
    alpha_value <- alpha_result$alpha_value
    thresh_results <- alpha_result$thresh_results
    
    for (thresh_result in thresh_results) {
      if (is.null(thresh_result)) next  # Skip null results
      #if (debug_progress) {
      #  print(paste("Updating roc_data_tra for threshold:", metrics_row$threshold))
      #}
      if (thresh_result$metrics_row$auc_tra_outter > optimal_auc_tra_outter) {
        filter_alpha = alpha_value
        filter_threshold = thresh_result$metrics_row$filter_threshold
        num_features = thresh_result$metrics_row$num_features
        auc_tra_inner <- thresh_result$metrics_row$auc_tra_inner
        auc_tra_outter <- thresh_result$metrics_row$auc_tra_outter # outter fold
        number_of_runs = number_of_runs
        filter_nfolds = thresh_result$metrics_row$filter_nfolds
        n_inner_folds = thresh_result$metrics_row$n_inner_folds
        n_outer_folds = thresh_result$metrics_row$n_outer_folds
        repeated_fit <- thresh_result$repeated_fit
        final_fit <- thresh_result$repeated_fit$final_fit
        alpha_lambda <- repeated_fit$final_param
        coefficients <- repeated_fit$final_coef
        roc_data_tra_inner = innercv_roc(repeated_fit) # Inner CV
        roc_data_val_outter = repeated_fit$roc # outter CV
      }
    }
  }
  ### Test best model on validation set:
  #selected_predictors = rownames(coefficients)[-1]
  
  # Try to extract the Dimnames
  try({
    selected_predictors = final_fit$glmnet.fit$beta@Dimnames[[1]]
  }, silent = TRUE)
  
  # Check if there was an error and print the error message
  if (is.null(selected_predictors)) {
    cat("Error: Failed to extract Dimnames\n")
    pred_val <- NULL
    roc_perf_val <- NULL
    auc_val <- NULL
  }
  
  norm_counts_data_validation_top <- norm_counts_data_validation[, selected_predictors, drop = FALSE]
  # Ensure only selected predictors present in the validation dataset are chosen
  common_predictors <- intersect(selected_predictors, colnames(norm_counts_data_validation))
  
  # Check if the lengths of the selected predictors and the common predictors match
  if (length(selected_predictors) == length(common_predictors)) {
    cat("Matched selected_preds and cols in norm_counts")
    #cat(paste("selected_predictors:",selected_predictors))
    #cat(paste("colnames(norm_counts_data_validation_top):",colnames(norm_counts_data_validation_top)))
    probabilities_val <- predict(final_fit, newx = as.matrix(norm_counts_data_validation_top), s = alpha_lambda[1], type = "response")
    pred_val <- prediction(probabilities_val, conditions_vector_validation)
    roc_perf_val <- performance(pred_val, "tpr", "fpr")
    auc_val <- as.numeric(performance(pred_val, measure = "auc")@y.values)
  }else{
    cat("Mismatched selected_preds and cols in norm_counts")
    #cat(paste("selected_predictors:",selected_predictors))
    #cat(paste("colnames(norm_counts_data_validation_top):",colnames(norm_counts_data_validation_top)))
    pred_val <- NULL
    roc_perf_val <- NULL
    auc_val <- NULL
  }

  return(
    list(
      part_seed = random_seed,
      filter_alpha = filter_alpha,
      filter_threshold = filter_threshold,
      num_features = num_features,
      auc_tra_inner = auc_tra_inner,
      auc_tra_outter = auc_tra_outter,
      auc_val = auc_val,
      number_of_runs = number_of_runs,
      filter_nfolds = filter_nfolds,
      n_inner_folds = n_inner_folds,
      n_outer_folds = n_outer_folds,
      alpha_lambda = alpha_lambda,
      coefficients = coefficients,
      roc_data_tra_inner = roc_data_tra_inner,
      roc_data_tra_outter = roc_data_val_outter,
      roc_data_val = roc_perf_val,
      repeated_fit = repeated_fit,
      final_fit = final_fit
    )
  )
}

# Function to flatten each result into a single data frame row from nested cv
flatten_result <- function(result) {
  # Helper function to safely extract elements or return NULL
  safe_extract <- function(x, field) {
    if (!is.null(x[[field]])) {
      return(x[[field]])
    } else {
      return(NA)
    }
  }
  
  # Extract coefficients names and convert to a single string
  coef_names <- if (!is.null(result$coefficients)) {
    paste(rownames(result$coefficients), collapse = "; ")
  } else {
    NA
  }
  
  # Extract alpha_lambda values and convert to a single string
  alpha_lambda_str <- if (!is.null(result$alpha_lambda)) {
    paste(result$alpha_lambda, collapse = "; ")
  } else {
    NA
  }
  
  # Create a data frame with the other elements up to coefficients and the concatenated coefficient names
  df <- data.frame(
    part_seed = safe_extract(result, "part_seed"),
    filter_alpha = safe_extract(result, "filter_alpha"),
    filter_threshold = safe_extract(result, "filter_threshold"),
    num_features = safe_extract(result, "num_features"),
    auc_tra_inner = safe_extract(result, "auc_tra_inner"),
    auc_tra_outter = safe_extract(result, "auc_tra_outter"),
    auc_val = safe_extract(result, "auc_val"),
    number_of_runs = safe_extract(result, "number_of_runs"),
    filter_nfolds = safe_extract(result, "filter_nfolds"),
    n_inner_folds = safe_extract(result, "n_inner_folds"),
    n_outer_folds = safe_extract(result, "n_outer_folds"),
    alpha_lambda = alpha_lambda_str,
    coefficients = coef_names,
    stringsAsFactors = FALSE
  )
  
  return(df)
}


subset_normalized_counts <- function(norm_counts_data, sig_results_genes) {
  if (nrow(norm_counts_data) == 0 && ncol(norm_counts_data) == 0) {
    print("Holdout matrix is empty.")
    return(NULL)
  }
  # Transpose the data for proper indexing
  norm_counts_data_t <- t(as.matrix(norm_counts_data))
  
  # Identify genes that are not in the normalized counts data
  missing_genes <- setdiff(sig_results_genes, colnames(norm_counts_data_t))
  
  # Debug statements
  if (length(missing_genes) > 0) {
    cat("Some sig_results_genes are not in norm_counts_data columns:", missing_genes, "\n")
  } else {
    cat("All sig_results_genes found in norm_counts_data columns.\n")
  }
  
  cat("Subset genes found:", sum(sig_results_genes %in% colnames(norm_counts_data_t)), "out of", length(sig_results_genes), "\n")
  
  # Subset the normalized counts data
  subset_data <- norm_counts_data_t[, sig_results_genes, drop = FALSE]
  
  # Check dimensions of the subset data
  cat("Dimensions of subset data:", dim(subset_data), "\n")
  
  return(subset_data)
}

# Function to create conditions vector
create_conditions_vector <- function(meta_data, column_name) {
  if (nrow(meta_data) == 0) {
    return(NULL)
  }
  
  # Get the unique class names in the specified column
  class_names <- unique(meta_data[[column_name]])
  
  # Initialize an empty vector
  conditions_vector <- c()
  
  # Assign numeric labels to each class
  for (i in seq_along(class_names)) {
    conditions_vector <- c(conditions_vector, rep(i-1, sum(meta_data[[column_name]] == class_names[i])))
  }
  
  return(conditions_vector)
}

evaluate_model <- function(best_model, best_params, norm_counts_data_training, norm_counts_data_validation, norm_counts_data_holdout, conditions_vector_training, conditions_vector_validation, conditions_vector_holdout) {
  print("Evaluating model")
  
  coefficients_matrix <- coef(best_model$finalModel, s = best_params$lambda)
  #non_zero_coefficients <- coefficients_matrix[coefficients_matrix != 0]
  #non_zero_coefficients_list <- as.data.frame(as.matrix(non_zero_coefficients))
  selected_predictors <- rownames(coefficients_matrix)[-1]
  
  #print("coefficients_matrix")
  #print(coefficients_matrix)
  #print("non_zero_coefficients")
  #print(non_zero_coefficients)
  #print("non_zero_coefficients_list")
  #print(non_zero_coefficients_list)
  #print("selected_predictors:")
  #print(selected_predictors)
  
  norm_counts_data_training_top <- norm_counts_data_training[, selected_predictors, drop = FALSE]
  norm_counts_data_validation_top <- norm_counts_data_validation[, selected_predictors, drop = FALSE]
  norm_counts_data_holdout_top <- norm_counts_data_holdout[, selected_predictors, drop = FALSE]
  
  print("Subsetting files...")
  
  probabilities_val <- predict(best_model$finalModel, newx = as.matrix(norm_counts_data_validation_top), s = best_params$lambda, type = "response")
  probabilities_hol <- predict(best_model$finalModel, newx = as.matrix(norm_counts_data_holdout_top), s = best_params$lambda, type = "response")
  probabilities_tra <- predict(best_model$finalModel, newx = as.matrix(norm_counts_data_training_top), s = best_params$lambda, type = "response")
  
  pred_val <- prediction(probabilities_val, conditions_vector_validation)
  pred_hol <- prediction(probabilities_hol, conditions_vector_holdout)
  pred_tra <- prediction(probabilities_tra, conditions_vector_training)
  
  roc_perf_val <- performance(pred_val, "tpr", "fpr")
  roc_perf_hol <- performance(pred_hol, "tpr", "fpr")
  roc_perf_tra <- performance(pred_tra, "tpr", "fpr")
  
  auc_tra <- as.numeric(performance(pred_tra, measure = "auc")@y.values)
  auc_val <- as.numeric(performance(pred_val, measure = "auc")@y.values)
  auc_hol <- as.numeric(performance(pred_hol, measure = "auc")@y.values)
  
  perf_sens <- performance(pred_val, measure = "sens")
  perf_spec <- performance(pred_val, measure = "spec")
  
  sens_values <- unlist(perf_sens@y.values)
  spec_values <- unlist(perf_spec@y.values)
  thresholds <- unlist(perf_sens@x.values)
  
  youden_j <- sens_values + spec_values - 1
  optimal_index <- which.max(youden_j)
  opt_threshold <- thresholds[optimal_index]
  cat("Optimal threshold (Youden's J):", opt_threshold, "\n")
  
  predicted_labels_tra <- ifelse(probabilities_tra >= opt_threshold, 1, 0)
  predicted_labels_val <- ifelse(probabilities_val >= opt_threshold, 1, 0)
  predicted_labels_hol <- ifelse(probabilities_hol >= opt_threshold, 1, 0)
  
  # Consistent storage for holdout ROC data with threshold keys
  roc_perf_hol_list <- list()
  roc_perf_hol_list[[as.character(opt_threshold)]] <- list(roc_perf_hol, auc_hol)
  roc_perf_tra_list <- list()
  roc_perf_tra_list[[as.character(opt_threshold)]] <- list(roc_perf_tra, auc_tra)
  roc_perf_val_list <- list()
  roc_perf_val_list[[as.character(opt_threshold)]] <- list(roc_perf_val, auc_val)
  
  list(
    roc_perf_tra = roc_perf_tra_list,
    roc_perf_val = roc_perf_val_list,
    roc_perf_hol = roc_perf_hol_list,
    auc_tra = auc_tra,
    auc_val = auc_val,
    auc_hol = auc_hol,
    opt_threshold = opt_threshold,
    non_zero_coefficients_list = selected_predictors,
    best_coefficients = coefficients_matrix
  )
}


## Visualization Functions
# Function to calculate AUC using the trapezoidal rule
calculate_auc <- function(fpr, tpr) {
  auc <- sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2)
  return(auc)
}

# Function to interpolate TPR values at common FPR points
interpolate_roc <- function(fpr, tpr, common_fpr) {
  approx(fpr, tpr, xout = common_fpr)$y
}

# Function to calculate the mean and standard deviation of TPR at common FPR points
calculate_mean_sd_tpr <- function(interpolated_roc_list) {
  all_tpr <- sapply(interpolated_roc_list, function(roc) roc$tpr)
  mean_tpr <- rowMeans(all_tpr)
  sd_tpr <- apply(all_tpr, 1, sd)
  list(mean_tpr = mean_tpr, sd_tpr = sd_tpr)
}

# Load required libraries

plot_nested_roc_curves_ARCHIVE <- function(roc_data_df, title, debug = FALSE) {
  library(doParallel)
  library(foreach)
  library(plotly)
  library(pROC)
  
  if(debug) cat("Starting plot_nested_roc_curves function\n")
  
  p <- plot_ly()
  
  # Convert each column of the data frame into a list of roc data
  roc_data_list <- apply(roc_data_df, 2, function(col) list(
    predictor = eval(parse(text = col["predictor"])),
    response = eval(parse(text = col["response"]))
  ))
  
  if(debug) {
    cat("Created roc_data_list with", length(roc_data_list), "elements\n")
    cat("Structure of first element in roc_data_list:\n")
    print(str(roc_data_list[[1]]))
  }
  
  # Combine all predictors and responses
  all_predictors <- unlist(lapply(roc_data_list, function(x) x$predictor))
  all_responses <- unlist(lapply(roc_data_list, function(x) x$response))
  
  # Create a single pooled ROC object
  pooled_roc <- roc(all_responses, all_predictors, quiet = TRUE)
  
  if(debug) cat("Created pooled ROC object\n")
  
  # Calculate pooled AUC
  pooled_auc <- auc(pooled_roc)
  
  if(debug) cat("Calculated pooled AUC:", pooled_auc, "\n")
  
  ### Apply smoothing by interpolating specificities with a given length
  # Create individual ROC curves for plotting and SD calculation
  roc_list <- lapply(roc_data_list, function(roc_data) {
    roc(roc_data$response, roc_data$predictor, quiet = TRUE)
  })
  
  if(debug) {
    cat("Created roc_list with", length(roc_list), "elements\n")
    cat("Structure of first element in roc_list:\n")
    print(str(roc_list[[1]]))
  }
  
  # Define common specificity points
  specificities <- seq(0, 1, length.out = 30)
  
  # Calculate sensitivities for each curve at the common specificity points
  sensitivities_list <- lapply(roc_list, function(roc_obj) {
    coords(roc_obj, x = specificities, input = "specificity", ret = "sensitivity")
  })
  
  if(debug) {
    cat("Created sensitivities_list with", length(sensitivities_list), "elements\n")
    cat("Structure of first element in sensitivities_list:\n")
    print(str(sensitivities_list[[1]]))
  }
  
  # Convert list to matrix
  sensitivities_matrix <- do.call(cbind, sensitivities_list)
  
  if(debug) {
    cat("Created sensitivities_matrix\n")
    cat("Dimensions of sensitivities_matrix:", dim(sensitivities_matrix), "\n")
    cat("First few rows and columns of sensitivities_matrix:\n")
    print(sensitivities_matrix[1:5, 1:5])
  }
  
  # Calculate mean and SD of sensitivities
  mean_sensitivities <- rowMeans(sensitivities_matrix, na.rm = TRUE)
  sd_sensitivities <- apply(sensitivities_matrix, 1, sd, na.rm = TRUE)
  
  if(debug) {
    cat("Calculated mean and SD of sensitivities\n")
    cat("Length of mean_sensitivities:", length(mean_sensitivities), "\n")
    cat("Length of sd_sensitivities:", length(sd_sensitivities), "\n")
  }
  
  # Calculate lower and upper bounds (2 SD)
  lower_tpr <- pmax(0, mean_sensitivities - 2 * sd_sensitivities)
  upper_tpr <- pmin(1, mean_sensitivities + 2 * sd_sensitivities)
  
  # Calculate FPR (1 - specificity)
  fpr <- 1 - specificities
  
  if(debug) {
    cat("FPR range:", range(fpr), "\n")
    cat("Mean TPR range:", range(mean_sensitivities), "\n")
    cat("Lower TPR range:", range(lower_tpr), "\n")
    cat("Upper TPR range:", range(upper_tpr), "\n")
  }
  
  # Plot individual ROC curves
  for (i in seq_along(roc_list)) {
    roc_obj <- roc_list[[i]]
    auc_value <- auc(roc_obj)
    
    p <- add_trace(p, x = 1 - roc_obj$specificities, y = roc_obj$sensitivities, 
                   type = 'scatter', mode = 'lines',
                   line = list(color = "rgba(169, 169, 169, 0.33)", width = 1.25),
                   #name = paste("Model", i, "AUC", round(auc_value, 3)))
                   showlegend = FALSE)
    
    if(debug && i %% 10 == 0) cat("Added trace for Model", i, "with AUC", auc_value, "\n")
  }
  
  # Add a single grey line to represent all individual models in the legend
  p <- add_trace(p, x = c(0, 1), y = c(0, 1), type = 'scatter', mode = 'lines',
                 line = list(color = "rgba(169, 169, 169, 0.33)", width = 1.25),
                 name = "Individual Models")
  
  if(debug) cat("Added representative line for individual models to legend\n")
  
  # Add the pooled ROC curve
  # Add line from origin to first point of pooled ROC curve
  p <- add_trace(p, x = c(0, fpr[length(fpr)]), y = c(0, mean_sensitivities[length(mean_sensitivities)]), 
                 type = 'scatter', mode = 'lines',
                 line = list(color = "rgb(75, 0, 130)", width = 3),
                 showlegend = FALSE)
  
  p <- add_trace(p, x = fpr, y = mean_sensitivities, type = 'scatter', mode = 'lines',
                 line = list(color = "rgb(75, 0, 130)", width = 3),
                 name = "Pooled ROC")
  
  if(debug) cat("Added pooled ROC curve\n")
  
  # Add the ribbon for confidence bands
  p <- add_ribbons(p, x = fpr, ymin = lower_tpr, ymax = upper_tpr, 
                   line = list(color = "transparent"), 
                   fillcolor = "rgba(153, 102, 204, 0.3)", 
                   name = "±2 SD")
  
  if(debug) cat("Added confidence interval ribbon\n")
  
  # Add the diagonal reference line
  p <- add_trace(p, x = c(0, 1), y = c(0, 1), type = 'scatter', mode = 'lines',
                 line = list(color = 'rgb(0, 0, 0)', dash = 'dash'), showlegend = FALSE)
  
  # Calculate CI based on SD of AUCs
  auc_values <- sapply(roc_list, auc)
  auc_sd <- sd(auc_values)
  ci_lower <- max(0, pooled_auc - 2 * auc_sd)
  ci_upper <- min(1, pooled_auc + 2 * auc_sd)
  
  # Add the Pooled AUC label with 95% CI based on SD
  p <- p %>% layout(
    annotations = list(
      text = paste("Pooled AUC =", round(pooled_auc, 3), " (n =", length(roc_list), ")", 
                   "\n(95% CI:", round(ci_lower, 3), "-", round(ci_upper, 3), ")"),
      x = 0.95, y = 0.05, xref = 'paper', yref = 'paper',
      showarrow = FALSE, font = list(size = 12, color = 'black'),
      xanchor = 'right', yanchor = 'bottom'
    )
  )
  
  if(debug) cat("Added annotation with pooled AUC and CI\n")
  
  # Add layout settings
  p <- p %>% layout(
    title = list(text = title, font = list(size = 20)),
    xaxis = list(title = 'False Positive Rate', titlefont = list(size = 18), showgrid = FALSE, tickfont = list(size = 16)),
    yaxis = list(title = 'True Positive Rate', titlefont = list(size = 18), showgrid = FALSE, tickfont = list(size = 16),
                 range = c(0, 1)),
    plot_bgcolor = 'rgba(0, 0, 0, 0)',
    xaxis2 = list(overlaying = 'x', showline = FALSE, showgrid = FALSE, zeroline = FALSE),
    yaxis2 = list(overlaying = 'y', showline = FALSE, showgrid = FALSE, zeroline = FALSE)
  )
  
  if(debug) cat("Added layout settings\n")
  
  if(debug) cat("Finished plot_nested_roc_curves function\n")
  
  # Return the plot object
  return(p)
}

plot_nested_roc_curves <- function(roc_data_df, title, debug = FALSE) {
  library(doParallel)
  library(foreach)
  library(plotly)
  library(pROC)
  
  if(debug) cat("Starting plot_nested_roc_curves function\n")
  
  p <- plot_ly()
  
  # Convert each column of the data frame into a list of roc data
  roc_data_list <- apply(roc_data_df, 2, function(col) list(
    predictor = eval(parse(text = col["predictor"])),
    response = eval(parse(text = col["response"]))
  ))
  
  if(debug) {
    cat("Created roc_data_list with", length(roc_data_list), "elements\n")
    cat("Structure of first element in roc_data_list:\n")
    print(str(roc_data_list[[1]]))
  }
  
  # Combine all predictors and responses
  all_predictors <- unlist(lapply(roc_data_list, function(x) x$predictor))
  all_responses <- unlist(lapply(roc_data_list, function(x) x$response))
  
  # Create a single pooled ROC object
  pooled_roc <- roc(all_responses, all_predictors, quiet = TRUE)
  
  if(debug) cat("Created pooled ROC object\n")
  
  # Calculate pooled AUC
  pooled_auc <- auc(pooled_roc)
  
  if(debug) cat("Calculated pooled AUC:", pooled_auc, "\n")
  
  ### Apply smoothing by interpolating TPR at common FPR points
  # Create individual ROC curves for plotting and SD calculation
  roc_list <- lapply(roc_data_list, function(roc_data) {
    roc(roc_data$response, roc_data$predictor, quiet = TRUE)
  })
  
  if(debug) {
    cat("Created roc_list with", length(roc_list), "elements\n")
    cat("Structure of first element in roc_list:\n")
    print(str(roc_list[[1]]))
  }
  
  # Define common FPR points
  fpr_common <- seq(0, 1, length.out = 100)
  
  # Interpolate TPR values at common FPR points for each ROC curve
  tpr_interpolated_list <- lapply(roc_list, function(roc_obj) {
    # Extract FPR and TPR from roc_obj
    fpr <- 1 - roc_obj$specificities
    tpr <- roc_obj$sensitivities
    # Interpolate TPR at common FPR points
    tpr_interpolated <- approx(fpr, tpr, xout = fpr_common, yleft = 0, yright = 1)$y
    return(tpr_interpolated)
  })
  
  if(debug) {
    cat("Created tpr_interpolated_list with", length(tpr_interpolated_list), "elements\n")
    cat("Structure of first element in tpr_interpolated_list:\n")
    print(str(tpr_interpolated_list[[1]]))
  }
  
  # Convert list to matrix
  tpr_matrix <- do.call(cbind, tpr_interpolated_list)
  
  if(debug) {
    cat("Created tpr_matrix\n")
    cat("Dimensions of tpr_matrix:", dim(tpr_matrix), "\n")
    cat("First few rows and columns of tpr_matrix:\n")
    print(tpr_matrix[1:5, 1:5])
  }
  
  # Calculate mean and SD of TPR
  mean_tpr <- rowMeans(tpr_matrix, na.rm = TRUE)
  sd_tpr <- apply(tpr_matrix, 1, sd, na.rm = TRUE)
  
  if(debug) {
    cat("Calculated mean and SD of TPR\n")
    cat("Length of mean_tpr:", length(mean_tpr), "\n")
    cat("Length of sd_tpr:", length(sd_tpr), "\n")
  }
  
  # Calculate lower and upper bounds (2 SD)
  lower_tpr <- pmax(0, mean_tpr - 2 * sd_tpr)
  upper_tpr <- pmin(1, mean_tpr + 2 * sd_tpr)
  
  # Assign fpr_common to fpr for consistency
  fpr <- fpr_common
  
  if(debug) {
    cat("FPR range:", range(fpr), "\n")
    cat("Mean TPR range:", range(mean_tpr), "\n")
    cat("Lower TPR range:", range(lower_tpr), "\n")
    cat("Upper TPR range:", range(upper_tpr), "\n")
  }
  
  # Plot individual ROC curves
  for (i in seq_along(roc_list)) {
    roc_obj <- roc_list[[i]]
    auc_value <- auc(roc_obj)
    
    p <- add_trace(p, x = 1 - roc_obj$specificities, y = roc_obj$sensitivities, 
                   type = 'scatter', mode = 'lines',
                   line = list(color = "rgba(169, 169, 169, 0.33)", width = 1.25),
                   showlegend = FALSE)
    
    if(debug && i %% 10 == 0) cat("Added trace for Model", i, "with AUC", auc_value, "\n")
  }
  
  # Add a single grey line to represent all individual models in the legend
  p <- add_trace(p, x = c(0, 1), y = c(0, 1), type = 'scatter', mode = 'lines',
                 line = list(color = "rgba(169, 169, 169, 0.33)", width = 1.25),
                 name = "Individual Models")
  
  if(debug) cat("Added representative line for individual models to legend\n")
  
  # Add the pooled ROC curve
  p <- add_trace(p, x = fpr, y = mean_tpr, type = 'scatter', mode = 'lines',
                 line = list(color = "rgb(75, 0, 130)", width = 3),
                 name = "Pooled ROC")
  
  if(debug) cat("Added pooled ROC curve\n")
  
  # Add the ribbon for confidence bands
  p <- add_ribbons(p, x = fpr, ymin = lower_tpr, ymax = upper_tpr, 
                   line = list(color = "transparent"), 
                   fillcolor = "rgba(153, 102, 204, 0.3)", 
                   name = "±2 SD")
  
  if(debug) cat("Added confidence interval ribbon\n")
  
  # Add the diagonal reference line
  p <- add_trace(p, x = c(0, 1), y = c(0, 1), type = 'scatter', mode = 'lines',
                 line = list(color = 'rgb(0, 0, 0)', dash = 'dash'), showlegend = FALSE)
  
  # Calculate CI based on SD of AUCs
  auc_values <- sapply(roc_list, auc)
  auc_sd <- sd(auc_values)
  ci_lower <- max(0, pooled_auc - 2 * auc_sd)
  ci_upper <- min(1, pooled_auc + 2 * auc_sd)
  
  # Add the Pooled AUC label with 95% CI based on SD
  p <- p %>% layout(
    annotations = list(
      text = paste("Pooled AUC =", round(pooled_auc, 3), " (n =", length(roc_list), ")", 
                   "\n(95% CI:", round(ci_lower, 3), "-", round(ci_upper, 3), ")"),
      x = 0.95, y = 0.05, xref = 'paper', yref = 'paper',
      showarrow = FALSE, font = list(size = 12, color = 'black'),
      xanchor = 'right', yanchor = 'bottom'
    )
  )
  
  if(debug) cat("Added annotation with pooled AUC and CI\n")
  
  # Add layout settings
  p <- p %>% layout(
    title = list(text = title, font = list(size = 20)),
    xaxis = list(title = 'False Positive Rate', titlefont = list(size = 18), showgrid = FALSE, tickfont = list(size = 16)),
    yaxis = list(title = 'True Positive Rate', titlefont = list(size = 18), showgrid = FALSE, tickfont = list(size = 16),
                 range = c(0, 1)),
    plot_bgcolor = 'rgba(0, 0, 0, 0)',
    xaxis2 = list(overlaying = 'x', showline = FALSE, showgrid = FALSE, zeroline = FALSE),
    yaxis2 = list(overlaying = 'y', showline = FALSE, showgrid = FALSE, zeroline = FALSE)
  )
  
  if(debug) cat("Added layout settings\n")
  
  if(debug) cat("Finished plot_nested_roc_curves function\n")
  
  # Return the plot object
  return(p)
}

# # Function to interpolate TPR values at common FPR points
# interpolate_rocr_roc <- function(fpr, tpr, common_fpr) {
#   approx(fpr, tpr, xout = common_fpr)$y
# }
# 
# # Function to calculate the mean and standard deviation of TPR at common FPR points for ROCR performance list
# calculate_mean_sd_rocr_tpr <- function(performance_list, common_fpr) {
#   interpolated_tpr_list <- lapply(performance_list, function(perf) {
#     fpr <- perf@x.values[[1]]
#     tpr <- perf@y.values[[1]]
#     interpolate_rocr_roc(fpr, tpr, common_fpr)
#   })
#   all_tpr <- do.call(cbind, interpolated_tpr_list)
#   mean_tpr <- rowMeans(all_tpr)
#   sd_tpr <- apply(all_tpr, 1, sd)
#   list(mean_tpr = mean_tpr, sd_tpr = sd_tpr)
# }

# plot_rocr_roc_curves <- function(performance_list, title, debug = TRUE) {
#   library(plotly)
#   library(ROCR)
#   
#   if(debug) cat("Starting plot_rocr_roc_curves function\n")
#   
#   p <- plot_ly()
#   
#   # Define common specificity points
#   specificities <- seq(0, 1, length.out = 33)
#   common_fpr <- 1 - specificities
#   
#   if(debug) cat("Defined common specificity points\n")
#   if(debug) cat("Common FPR range:", common_fpr, "\n")
#   if(debug) cat("specificities:", specificities, "\n")
#   
#   # Extract sensitivities at common specificity points for each curve
#   sensitivities_list <- lapply(performance_list, function(perf) {
#     fpr <- perf@x.values[[1]]
#     tpr <- perf@y.values[[1]]
#     approx(fpr, tpr, xout = common_fpr)$y
#   })
#   
#   if(debug) {
#     cat("Extracted sensitivities for each curve\n")
#     cat("Number of curves:", length(sensitivities_list), "\n")
#   }
#   
#   # Convert list to matrix
#   sensitivities_matrix <- do.call(cbind, sensitivities_list)
#   
#   if(debug) {
#     cat("Created sensitivities matrix\n")
#     cat("Dimensions of sensitivities_matrix:", dim(sensitivities_matrix), "\n")
#   }
#   
#   # Calculate mean and SD of sensitivities
#   mean_tpr <- rowMeans(sensitivities_matrix, na.rm = TRUE)
#   sd_tpr <- apply(sensitivities_matrix, 1, sd, na.rm = TRUE)
#   
#   if(debug) cat("Calculated mean and SD of sensitivities\n")
#   
#   # Calculate lower and upper bounds (2 SD)
#   ribbon_tpr_lower <- pmax(0, mean_tpr - 2 * sd_tpr)
#   ribbon_tpr_upper <- pmin(1, mean_tpr + 2 * sd_tpr)
#   
#   # Calculate pooled AUC
#   pooled_auc <- calculate_auc(common_fpr, mean_tpr)
#   
#   if(debug) cat("Calculated pooled AUC:", pooled_auc, "\n")
#   
#   # Calculate individual AUCs and their SD
#   individual_aucs <- sapply(performance_list, function(perf) {
#     calculate_auc(perf@x.values[[1]], perf@y.values[[1]])
#   })
#   auc_sd <- sd(individual_aucs)
#   
#   # Calculate CI based on SD of AUCs
#   ci_lower <- max(0, pooled_auc - 2 * auc_sd)
#   ci_upper <- min(1, pooled_auc + 2 * auc_sd)
#   
#   
#   if(debug) {
#     cat("Calculated CI based on SD of AUCs\n")
#     cat("CI:", ci_lower, "-", ci_upper, "\n")
#   }
#   
#   # Plot individual ROC curves
#   for (i in seq_along(performance_list)) {
#     perf <- performance_list[[i]]
#     fpr <- perf@x.values[[1]]
#     tpr <- perf@y.values[[1]]
#     auc <- individual_aucs[i]
#     
#     p <- add_trace(p, x = fpr, y = tpr, type = 'scatter', mode = 'lines',
#                    line = list(color = "rgba(169, 169, 169, 0.33)", width = 1),
#                    name = paste("Model", i, "AUC", round(auc, 3)), showlegend = FALSE)
#     
#     if(debug && i %% 10 == 0) cat("Added trace for Model", i, "with AUC", auc, "\n")
#   }
#   
#   # Add a single grey line to represent all individual models in the legend
#   p <- add_trace(p, x = c(0, 1), y = c(0, 1), type = 'scatter', mode = 'lines',
#                  line = list(color = "rgba(169, 169, 169, 0.33)", width = 1.25),
#                  name = "Individual Models")
#   
#   if(debug) cat("Added representative line for individual models to legend\n")
#   
#   # Add the pooled ROC curve
#   p <- add_trace(p, x = common_fpr, y = mean_tpr, type = 'scatter', mode = 'lines',
#                  line = list(color = "rgb(75, 0, 130)", width = 3),
#                  name = "Pooled ROC")
#   
#   if(debug) cat("Added pooled ROC curve\n")
#   
#   # Add the ribbon for 2 standard deviations
#   p <- add_ribbons(p, x = common_fpr, ymin = ribbon_tpr_lower, ymax = ribbon_tpr_upper, 
#                    line = list(color = "transparent"), 
#                    fillcolor = "rgba(153, 102, 204, 0.3)",
#                    name = "±2 SD")
#   
#   if(debug) cat("Added confidence interval ribbon\n")
#   
#   # Add the diagonal reference line
#   p <- add_trace(p, x = c(0, 1), y = c(0, 1), type = 'scatter', mode = 'lines',
#                  line = list(color = 'rgb(0, 0, 0)', dash = 'dash'), showlegend = FALSE)
#   
#   # Add the Pooled AUC label with 95% CI based on SD
#   p <- p %>% layout(
#     annotations = list(
#       text = paste("Pooled AUC =", round(pooled_auc, 3), " (n =", length(performance_list), ")", 
#                    "\n(95% CI:", round(ci_lower, 3), "-", round(ci_upper, 3), ")"),
#       x = 0.95, y = 0.05, xref = 'paper', yref = 'paper',
#       showarrow = FALSE, font = list(size = 12, color = 'black'),
#       xanchor = 'right', yanchor = 'bottom'
#     )
#   )
#   
#   if(debug) cat("Added annotation with pooled AUC and CI\n")
#   
#   # Add layout settings
#   p <- p %>% layout(
#     title = list(text = title, font = list(size = 20)),
#     xaxis = list(title = 'False Positive Rate', titlefont = list(size = 18), showgrid = FALSE, tickfont = list(size = 16),
#                  range = c(0, 1)),
#     yaxis = list(title = 'True Positive Rate', titlefont = list(size = 18), showgrid = FALSE, tickfont = list(size = 16),
#                  range = c(0, 1)),
#     plot_bgcolor = 'rgba(0, 0, 0, 0)',
#     xaxis2 = list(overlaying = 'x', showline = FALSE, showgrid = FALSE, zeroline = FALSE),
#     yaxis2 = list(overlaying = 'y', showline = FALSE, showgrid = FALSE, zeroline = FALSE)
#   )
#   
#   if(debug) cat("Added layout settings\n")
#   
#   if(debug) cat("Finished plot_rocr_roc_curves function\n")
#   
#   # Return the plot object
#   return(p)
# }
# 
# plot_rocr_roc_curves_DEPR <- function(performance_list, title) {
#   p <- plot_ly()
#   
#   # Define common FPR points for interpolation
#   common_fpr <- seq(0, 1, length.out = 100)
#   
#   # Calculate the mean and standard deviation of TPR at common FPR points
#   mean_sd_tpr <- calculate_mean_sd_rocr_tpr(performance_list, common_fpr)
#   mean_tpr <- mean_sd_tpr$mean_tpr
#   sd_tpr <- mean_sd_tpr$sd_tpr
#   
#   # Extract FPR and TPR values for the ribbon
#   ribbon_tpr_upper <- mean_tpr + 2 * sd_tpr
#   ribbon_tpr_lower <- mean_tpr - 2 * sd_tpr
#   
#   # Calculate AUC for the upper and lower ribbons
#   auc_upper <- calculate_auc(common_fpr, ribbon_tpr_upper)
#   auc_lower <- calculate_auc(common_fpr, ribbon_tpr_lower)
#   mean_auc <- calculate_auc(common_fpr, mean_tpr)
#   num_models <- length(performance_list)
#   
#   for (i in seq_along(performance_list)) {
#     perf <- performance_list[[i]]
#     fpr <- perf@x.values[[1]]
#     tpr <- perf@y.values[[1]]
#     
#     # Calculate the AUC using the trapezoidal rule
#     auc <- calculate_auc(fpr, tpr)
#     
#     line_color <- "rgba(169, 169, 169, 0.33)"  # 10% transparent dark purple
#     line_width <- 1
#     
#     p <- add_trace(p, x = fpr, y = tpr, type = 'scatter', mode = 'lines',
#                    line = list(color = line_color, width = line_width),
#                    name = paste("Model", i, "AUC", round(auc, 3)))  # Use AUC value in the label
#   }
#   
#   # Add the average ROC curve
#   p <- add_trace(p, x = common_fpr, y = mean_tpr, type = 'scatter', mode = 'lines',
#                  line = list(color = "rgb(75, 0, 130)", width = 3),
#                  name = "Average ROC")
#   
#   # Add the ribbon for 2 standard deviations
#   p <- add_ribbons(p, x = common_fpr, ymin = ribbon_tpr_lower, ymax = ribbon_tpr_upper, 
#                    line = list(color = "transparent"), 
#                    fillcolor = "rgba(153, 102, 204, 0.3)", # Light, transparent purple
#                    name = "2 Standard Deviations")
#   
#   # Add the diagonal reference line
#   p <- add_trace(p, x = c(0, 1), y = c(0, 1), type = 'scatter', mode = 'lines',
#                  line = list(color = 'rgb(0, 0, 0)', dash = 'dash'), showlegend = FALSE)  # Diagonal line color same as dark blue
#   
#   # Add the Avg. AUC label with 95% CI
#   p <- p %>% layout(
#     annotations = list(
#       text = paste("Mean AUC =", round(mean_auc, 3), " (n =", num_models, ")", "\n(95% CI:", round(auc_lower, 3), "-", round(auc_upper, 3), ")"),
#       x = 0.95, y = 0.05, xref = 'paper', yref = 'paper',
#       showarrow = FALSE, font = list(size = 12, color = 'black'),
#       xanchor = 'right', yanchor = 'bottom'
#     )
#   )
#   
#   # Add layout settings with thicker axis lines and larger fonts
#   p <- p %>% layout(
#     title = list(text = title, font = list(size = 20)),
#     xaxis = list(title = 'False Positive Rate', titlefont = list(size = 18), showgrid = FALSE, tickfont = list(size = 16),
#                  # set x axis lim from 0 to 1
#                  range = c(0, 1)),
#     yaxis = list(title = 'True Positive Rate', titlefont = list(size = 18), showgrid = FALSE, tickfont = list(size = 16),
#                  range = c(0, 1)),
#     plot_bgcolor = 'rgba(0, 0, 0, 0)',
#     xaxis2 = list(overlaying = 'x', showline = FALSE, showgrid = FALSE, zeroline = FALSE),
#     yaxis2 = list(overlaying = 'y', showline = FALSE, showgrid = FALSE, zeroline = FALSE)
#     # set x and y axis lim from 0 to 1
#     
#   )
#   
#   # Return the plot object
#   return(p)
# }

# Function to initialize results dataframes
initialize_results_dfs <- function() {
  metrics_df <- data.frame(
    model_id = integer(),
    random_seed = integer(),
    alpha = numeric(),
    threshold = numeric(),
    num_features = integer(),
    auc_tra = numeric(),
    auc_val = numeric(),
    auc_hol = numeric(),
    sensitivity_tra = numeric(),
    specificity_tra = numeric(),
    recall_tra = numeric(),
    precision_tra = numeric(),
    accuracy_tra = numeric(),
    f1_tra = numeric(),
    sensitivity_val = numeric(),
    specificity_val = numeric(),
    recall_val = numeric(),
    precision_val = numeric(),
    accuracy_val = numeric(),
    f1_val = numeric(),
    sensitivity_hol = numeric(),
    specificity_hol = numeric(),
    recall_hol = numeric(),
    precision_hol = numeric(),
    accuracy_hol = numeric(),
    f1_hol = numeric(),
    optimal_alpha = numeric(),
    optimal_threshold = numeric(),
    num_top_genes = integer(),
    nfolds_val = integer(),
    root_folder = character()
  )
  
  coefficients_df <- data.frame(
    model_id = integer(),
    feature = character(),
    coefficient = numeric()
  )
  
  return(list(metrics_df = metrics_df, coefficients_df = coefficients_df))
}

# Function to create a vertical bar plot with sorted bars and Plotly white theme
create_vertical_bar_plot <- function(feature_counts, median_model_features, percentage_cutoff = 50) {
  # Filter the data by the percentage cut-off
  feature_counts <- feature_counts[feature_counts$percentage >= percentage_cutoff,]
  
  # Sort the data by percentage in descending order
  feature_counts <- feature_counts[order(-feature_counts$percentage),]
  
  # Define the plotly plot object
  p <- plot_ly()
  
  # Add the bar traces
  p <- add_trace(p, data = feature_counts,
                 x = ~percentage,
                 y = ~reorder(GeneSymbol, percentage),
                 type = 'bar',
                 orientation = 'h',
                 marker = list(color = ifelse(feature_counts$feature %in% median_model_features, 'purple', 'grey'),
                               colorscale = 'Viridis'),
                 name = "Feature Count")
  
  # Add layout settings with thicker axis lines and larger fonts
  p <- p %>% layout(
    title = list(text = "Percentage of Models Each Feature Appeared In", font = list(size = 20)),
    xaxis = list(title = 'Percentage of Models', titlefont = list(size = 18), showgrid = FALSE, tickfont = list(size = 16)),
    yaxis = list(title = 'Feature', titlefont = list(size = 18), showgrid = FALSE, tickfont = list(size = 16)),
    plot_bgcolor = 'rgba(0, 0, 0, 0)',
    template = 'plotly_white'
  )
  
  # Return the plot object
  return(p)
}

create_boxplot <- function(coefficients_df, median_model_features, feature_counts, percentage_cutoff = 50) {
  # Filter the coefficients based on the percentage cutoff
  filtered_features <- feature_counts$feature[feature_counts$percentage >= percentage_cutoff]
  coefficients_df <- coefficients_df[coefficients_df$feature %in% filtered_features, ]
  
  # Calculate average coefficient for each feature and sort by descending order
  average_coefficients <- aggregate(coefficient ~ feature + GeneSymbol, data = coefficients_df, FUN = mean)
  average_coefficients <- average_coefficients[order(-average_coefficients$coefficient), ]
  
  # Reorder coefficients_df based on the sorted features
  coefficients_df <- coefficients_df[order(match(coefficients_df$feature, average_coefficients$feature)), ]
  
  # Set factors to ensure the order is maintained in the plot
  coefficients_df$feature <- factor(coefficients_df$feature, levels = average_coefficients$feature)
  coefficients_df$GeneSymbol <- factor(coefficients_df$GeneSymbol, levels = average_coefficients$GeneSymbol)
  
  # Create the box plot
  box_plot <- plot_ly(data = coefficients_df,
                      x = ~GeneSymbol,
                      y = ~coefficient,
                      type = 'box',
                      boxpoints = 'none',
                      #boxpoints = 'all',
                      #jitter = 0.3,
                      #pointpos = -1.8,
                      marker = list(color = ifelse(coefficients_df$feature %in% median_model_features, 'purple', 'grey'))) %>%
    layout(title = "Coefficient Distribution by Feature",
           xaxis = list(title = "GeneSymbol"),
           yaxis = list(title = "Coefficient"),
           template = 'plotly_white')
  
  return(box_plot)
}

# Function to annotate features with additional debugging
annotate_features <- function(df) {
  # Split the features into chromosome, start, and end positions
  chr_positions <- strsplit(as.character(df$feature), "_")
  chr <- sapply(chr_positions, function(x) x[1])
  start <- as.numeric(sapply(chr_positions, function(x) x[2]))
  end <- as.numeric(sapply(chr_positions, function(x) x[3]))
  
  # Create a GRanges object
  gr <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = end))
  
  # Load the TxDb annotation database
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  
  # Annotate the peaks with gene information
  peakAnno <- annotatePeak(gr, tssRegion = c(-3000, 1000), TxDb = txdb, annoDb = "org.Hs.eg.db", overlap = "all")
  
  # Extract the gene symbols and location types
  gene_symbols <- as.character(peakAnno@anno$SYMBOL)
  location_types <- as.character(peakAnno@anno$annotation)
  
  # Replace missing gene symbols with "N/A"
  gene_symbols[is.na(gene_symbols)] <- "N/A"
  
  # Debug: Print regions with N/A gene symbols and matched regions
  na_indices <- which(gene_symbols == "N/A")
  if (length(na_indices) > 0) {
    cat("Regions with N/A gene symbols:\n")
    for (i in na_indices) {
      cat("Feature:", df$feature[i], "\n")
      #cat("Matched regions:\n")
      #match <- peakAnno@anno[which(peakAnno@anno$SYMBOL == ""), ]
      #print(match)
    }
  }
  
  # Simplify location types
  simplify_location_type <- function(location) {
    if (grepl("Intron", location)) {
      return("Intron")
    } else if (grepl("Exon", location)) {
      return("Exon")
    } else if (grepl("Promoter", location)) {
      return("Promoter")
    } else if (grepl("Distal Intergenic", location)) {
      return("Distal Intergenic")
    } else if (grepl("3' UTR", location)) {
      return("3' UTR")
    } else if (grepl("5' UTR", location)) {
      return("5' UTR")
    } else if (grepl("Downstream", location)) {
      return("Downstream")
    } else {
      return("Unknown")
    }
  }
  
  location_types <- sapply(location_types, simplify_location_type)
  
  # Handle multiple gene symbols
  unique_symbols <- unique(gene_symbols)
  symbol_count <- sapply(unique_symbols, function(symbol) sum(gene_symbols == symbol))
  symbol_index <- setNames(rep(0, length(unique_symbols)), unique_symbols)
  
  for (i in seq_along(gene_symbols)) {
    symbol <- gene_symbols[i]
    if (symbol_count[symbol] > 1) {
      symbol_index[symbol] <- symbol_index[symbol] + 1
      gene_symbols[i] <- paste0(symbol, "(", symbol_index[symbol], ")")
    } else {
      gene_symbols[i] <- symbol
    }
  }
  
  # Get gene ranges
  gene_ranges <- unlist(genes(txdb, single.strand.genes.only = FALSE))
  
  # Find overlaps
  overlaps <- findOverlaps(gr, gene_ranges)
  overlapping_peaks <- queryHits(overlaps)
  
  # Initialize gene_dist with NA
  gene_dist <- rep(NA, length(gr))
  
  # Set distance to 0 for overlapping peaks
  gene_dist[overlapping_peaks] <- 0
  
  # Calculate distances for non-overlapping peaks
  non_overlapping_peaks <- setdiff(seq_along(gr), overlapping_peaks)
  if (length(non_overlapping_peaks) > 0) {
    distances_to_start <- distanceToNearest(gr[non_overlapping_peaks], gene_ranges)
    distances_to_end <- distanceToNearest(gr[non_overlapping_peaks], gene_ranges, ignore.strand = TRUE)
    
    min_distances <- pmin(mcols(distances_to_start)$distance, mcols(distances_to_end)$distance)
    
    gene_dist[non_overlapping_peaks] <- min_distances
  }
  
  # Add GeneSymbol, LocationType, and gene_dist to the dataframe
  df$GeneSymbol <- gene_symbols
  df$LocationType <- location_types
  df$gene_dist <- gene_dist
  
  return(df)
}

annotate_features2 <- function(df, annots = c('hg19_basicgenes','hg19_custom_genehancer')) { #'hg19_lncrna_gencode',
  # Force R to use standard notation for large numbers
  options(scipen = 999)
  
  # Split the features into chromosome, start, and end positions
  chr_positions <- strsplit(as.character(df$feature), "_")
  annotated_df <- df
  # intialize chr, start and end columns in annotated_df
  annotated_df$chr <- sapply(chr_positions, function(x) x[1])
  annotated_df$start <- as.numeric(sapply(chr_positions, function(x) x[2]))
  annotated_df$end <- as.numeric(sapply(chr_positions, function(x) x[3]))

  
  # message if any start or end == NA
  if (any(is.na(annotated_df$start)) || any(is.na(annotated_df$end))) {
    print("Some start or end positions are NA.")
  }
  
  # Print rows where start or end are na
  if (any(is.na(annotated_df$start))) {
    na_annotated_df <- annotated_df[is.na(annotated_df$start),]
  }
  # Create a GRanges object
  gr <- GRanges(seqnames = annotated_df$chr, ranges = IRanges(start = annotated_df$start, end = annotated_df$end))
  
  ### Load custom annotations datbase
  read_annotations(con = '~/reference_genomes/GeneHancer_AnnotSV_hg19_v5.20.txt', genome = 'hg19', name = 'genehancer', format = 'bed', extraCols = c(symbol="character",type="character"))
  #print(annotatr_cache$get('hg19_custom_genehancer'))
  
  # Select annotations for intersection with regions
  
  
  # Build the annotations (a single GRanges object)
  annotations <- build_annotations(genome = 'hg19', annotations = annots)
  
  # Intersect the regions we read in with the annotations
  dm_annotated <- annotate_regions(
    regions = gr,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE
  )
  # Coerce to a data.frame
  df_dm_annotated <- data.frame(dm_annotated)
  
  df_dm_annotated$start <- as.integer(df_dm_annotated$start)
  df_dm_annotated$end <- as.integer(df_dm_annotated$end)
  
  # drop dplicte rows where annot.symbol annot.type and start are the same
  df_dm_annotated <- df_dm_annotated[!duplicated(df_dm_annotated[,c('annot.symbol','annot.type','start')]),]
  #browser()
  # sort annot.type in the following order hg19_genes_exons > hg19_genes_introns > hg19_genes_promoters > hg19_genes_5UTRs > hg19_genes_3UTRs > hg19_genes_1to5kb > hg19_lncrna_gencode > hg19_custom_genehancer
  df_dm_annotated$annot.type <- factor(df_dm_annotated$annot.type, levels = c('hg19_genes_5UTRs','hg19_genes_3UTRs','hg19_genes_exons','hg19_genes_introns','hg19_genes_promoters','hg19_genes_1to5kb','hg19_lncrna_gencode','hg19_custom_genehancer'))
  # Sort the dataframe by annot.type
  df_dm_annotated <- df_dm_annotated[order(df_dm_annotated$annot.type),]
  
  # Drop rows with annot.symbol starting with "MIR" if there are other rows with the same seqnames, start, and end
  # Identify groups with the same seqnames, start, and end
  #grouped <- df_dm_annotated %>% 
  #  group_by(seqnames, start, end) %>%
  #  filter(any(!grepl("^MIR", annot.symbol)))
  # Remove rows with annot.symbol starting with "MIR" in the identified groups
  # df_dm_annotated <- grouped %>%
  #   filter(!grepl("^MIR", annot.symbol))
  
  # Group rows by seqnames, start and end, Add a column called all_symbols that represents list of all annot.symbol for same seqnames, start and end
  df_dm_annotated <- df_dm_annotated %>%
    group_by(seqnames, start, end) %>%
    mutate(all_symbols = paste(annot.symbol, collapse = ", "))
  
  # In all_symbols column, drop duplicate symbols within each row
  df_dm_annotated$all_symbols <- sapply(df_dm_annotated$all_symbols, function(x) paste(unique(unlist(strsplit(x, ", "))), collapse = ", "))
  
  # Now drop duplicates based on chr start and end, keeping the first item that appears based on order abov
  df_dm_annotated <- df_dm_annotated[!duplicated(df_dm_annotated[,c('seqnames','start','end')]),]
  
  # Drop rows where annot.symbol is NA
  df_dm_annotated <- df_dm_annotated[!is.na(df_dm_annotated$annot.symbol),]
  
  annotated_df <- merge(annotated_df, df_dm_annotated, by.x = c("chr", "start", "end"), by.y = c("seqnames", "start", "end"), all.x = TRUE)
  
  annotated_df$start <- as.integer(annotated_df$start)
  annotated_df$end <- as.integer(annotated_df$end)
  
  # Replace na annot.symbol with Intergenic: ", chr[i], "_", start[i], "_", end[i])" in df_annotated
  annotated_df$annot.symbol[is.na(annotated_df$annot.symbol)] <- paste(annotated_df$chr[is.na(annotated_df$annot.symbol)], "_", annotated_df$start[is.na(annotated_df$annot.symbol)], "_", annotated_df$end[is.na(annotated_df$annot.symbol)], sep = "")
  
  # if all_symbols == NA then replace with Intergenic: ", chr[i], "_", start[i], "_", end[i])" in df_annotated
  annotated_df$all_symbols[is.na(annotated_df$all_symbols)
  ] <- paste(annotated_df$chr[is.na(annotated_df$all_symbols)], "_", annotated_df$start[is.na(annotated_df$all_symbols)], "_", annotated_df$end[is.na(annotated_df$all_symbols)], sep = "")
  
  
  # Simplify location types
  simplify_location_type <- function(location) {
    if (grepl("genes_promoters", location)) {
      return("Promoter")
    } else if (grepl("genes_1to5kb", location)) {
      return("1-5kb Upstream")
    } else if (grepl("genes_5UTRs", location)) {
      return("5' UTR")
    } else if (grepl("genes_3UTRs", location)) {
      return("3' UTR")
    } else if (grepl("genes_exons", location)) {
      return("Exon")
    } else if (grepl("genes_introns", location)) {
      return("Intron")
    }  else if (grepl("lncrna_gencode", location)) {
      return("lncRNA")
    } else if (grepl("genehancer", location)) {
      return("Enhancer")
    } else {
      return("Intergenic")
    }
  }

  annotated_df$locationType <- sapply(annotated_df$annot.type, simplify_location_type)

  # Add GeneSymbol, LocationType, to df and return it
  df$GeneSymbol <- annotated_df$annot.symbol
  df$LocationType <- annotated_df$locationType

  return(annotated_df)
}

print_confusion_matrix <- function(results_list, median_model_idx, norm_counts_data, conditions_vector, set_name) {
  # Extract the median model, best parameters, and optimal threshold
  median_model <- results_list[[median_model_idx]]$best_model
  best_params <- results_list[[median_model_idx]]$best_params
  opt_threshold <- results_list[[median_model_idx]]$opt_threshold
  
  # Subset data based on selected predictors
  selected_predictors <- rownames(best_params)[-1]
  norm_counts_data_top <- norm_counts_data[, selected_predictors, drop = FALSE]
  
  # Predict probabilities
  probabilities <- predict(median_model$finalModel, newx = as.matrix(norm_counts_data_top), s = best_params$lambda, type = "response")
  
  # Generate predicted labels based on the optimal threshold
  predicted_labels <- ifelse(probabilities >= opt_threshold, 1, 0)
  
  # Print confusion matrix
  cat(paste("Confusion Matrix for", set_name, "Set:\n"))
  print(confusionMatrix(factor(predicted_labels), factor(conditions_vector)))
}

### Functions for plotting netedCV variable stability plot:
# Function to calculate pooled SEM
pooled_sem <- function(sem, n) {
  weighted_variance <- sum((sem^2) * n) / sum(n)
  return(sqrt(weighted_variance))
}

# Function to create the stability plot (REQUIRES FINAL MODEL TRAIN)
create_stability_plot_ARCHIVE <- function(results_list, min_freq_threshold, max_runs) {
  # insert debugging breakpoint
  browser()
  
  library(tibble)
  stability_df <- data.frame()
  num_seeds <- length(results_list)
  
  # Convert threshold to percentage
  min_freq_threshold <- (min_freq_threshold / max_runs) * 100
  
  # Build the coefficients table
  for (item in results_list) {
    var_stability_df <- var_stability(item$repeated_fit, percent = FALSE, sort = TRUE) %>%
      rownames_to_column(var = "feature")
    var_stability_df$model_id <- item$part_seed
    stability_df <- rbind(stability_df, var_stability_df)
  }
  
  
  # drop rows where stability_df final == no
  # I do not want to include this, as I want to calculate % of times the feature appears in final model
  #stability_df <- stability_df[stability_df$final == "yes",]
  
  # Summarize by feature
  # Group by feature and calculate pooled statistics, count_in_final, and percent_in_final
  stability_df_grouped <- stability_df %>%
    group_by(feature) %>%
    summarise(
      mean = mean(mean, na.rm = TRUE),
      pooled_sem = pooled_sem(sem, frequency),
      frequency = sum(frequency, na.rm = TRUE),
      #count_in_final = sum(frequency, na.rm = TRUE), # Assuming frequency represents count_in_final
      #percent_in_final = (sum(frequency, na.rm = TRUE) / (num_seeds * max_runs)) * 100,
      .groups = 'drop'
    )
  
  # Annotate features (assuming annotate_features is defined)
  stability_df_grouped <- annotate_features2(stability_df_grouped)
  
  
  # Convert frequency to percentage
  stability_df_grouped <- stability_df_grouped %>%
    mutate(frequency_percentage = (frequency / (max_runs)) * 100)
  
  # Sort by abs(mean) and filter by frequency percentage
  stability_df_grouped <- stability_df_grouped %>%
    arrange(abs(mean)) %>%
    filter(frequency_percentage > min_freq_threshold)
  
  # Prepare the y-axis labels
  stability_df_grouped <- stability_df_grouped %>%
    #mutate(label = paste(annot.symbol, "(", locationType, ")", sep = ""))
    mutate(label = paste(all_symbols, sep = " "))
  
  # Add this new code block here
  stability_df_grouped <- stability_df_grouped %>%
    group_by(label) %>%
    mutate(label = if(n() > 1) paste0(label, "_", row_number()) else label) %>%
    ungroup()
  
  stability_df_grouped$label <- factor(stability_df_grouped$label, levels = stability_df_grouped$label)
  
  # Print the grouped data frame to debug
  print(stability_df_grouped)
  
  # Create the plot
  all_model_stability <- plot_ly(
    data = stability_df_grouped,
    x = ~mean,
    y = ~label,
    type = 'scatter',
    mode = 'markers',
    marker = list(
      size = ~12.5 * sqrt(pmax(pmin(frequency_percentage, max(stability_df_grouped$frequency_percentage)) / max(stability_df_grouped$frequency_percentage))),
      color = ~pmax(pmin(frequency_percentage, max(stability_df_grouped$frequency_percentage)), min(stability_df_grouped$frequency_percentage)),
      colorscale = 'Viridis',
      showscale = TRUE,
      colorbar = list(title = 'Frequency (%)')
    ),
    error_x = ~list(
      array = pooled_sem,
      color = 'black',
      thickness = 1
    )
  ) %>%
    layout(
      title = 'Variable Importance Plot',
      xaxis = list(title = 'Mean (Variable Importance)'),
      yaxis = list(title = 'Gene Symbol (Location Type)'),
      margin = list(b = 100)
    )
  
  return(list(all_model_stability, stability_df_grouped))
}

# DOES NOT REQUIRE FINAL MODEL
create_stability_plot <- function(results_list, min_freq_threshold, max_runs) {
  # Insert a manual version of var_stability logic
  manual_var_stability <- function(obj, percent = FALSE, sort = TRUE) {
    # obj should be a nestcv.glmnet object or similar
    
    # Extract coefficients from outer folds
    coef_mat <- cv_coef(obj) 
    # coef_mat: rows = features, columns = folds
    
    # Compute statistics
    means <- rowMeans(coef_mat, na.rm = TRUE)
    sds <- apply(coef_mat, 1, sd, na.rm = TRUE)
    sem <- sds / sqrt(ncol(coef_mat))
    frequency <- rowSums(coef_mat != 0)
    
    df <- data.frame(
      mean = means,
      sd = sds,
      sem = sem,
      frequency = frequency,
      row.names = rownames(coef_mat),
      stringsAsFactors = FALSE
    )
    
    # Sorting by absolute mean importance if requested
    if (sort) {
      df <- df[order(abs(df$mean), decreasing = TRUE), , drop = FALSE]
    }
    
    return(df)
  }
  
  library(tibble)
  stability_df <- data.frame()
  num_seeds <- length(results_list)
  
  # Convert threshold to percentage
  min_freq_threshold <- (min_freq_threshold / max_runs) * 100
  
  # Build the coefficients table
  for (item in results_list) {
    # Instead of var_stability, we call our manual_var_stability()
    var_stability_df <- manual_var_stability(item$repeated_fit, percent = FALSE, sort = TRUE) %>%
      rownames_to_column(var = "feature")
    var_stability_df$model_id <- item$part_seed
    stability_df <- rbind(stability_df, var_stability_df)
  }
  
  # Summarize by feature
  stability_df_grouped <- stability_df %>%
    group_by(feature) %>%
    summarise(
      mean = mean(mean, na.rm = TRUE),
      pooled_sem = pooled_sem(sem, frequency),  # Assuming pooled_sem is defined elsewhere
      frequency = sum(frequency, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Annotate features (assuming annotate_features2 is defined)
  stability_df_grouped <- annotate_features2(stability_df_grouped)
  
  # Convert frequency to percentage
  stability_df_grouped <- stability_df_grouped %>%
    mutate(frequency_percentage = (frequency / (max_runs)) * 100)
  
  # Sort by abs(mean) and filter by frequency percentage
  stability_df_grouped <- stability_df_grouped %>%
    arrange(abs(mean)) %>%
    filter(frequency_percentage > min_freq_threshold)
  
  # Prepare the y-axis labels
  stability_df_grouped <- stability_df_grouped %>%
    mutate(label = paste(all_symbols, sep = " ")) %>%
    group_by(label) %>%
    mutate(label = if(n() > 1) paste0(label, "_", row_number()) else label) %>%
    ungroup()
  
  stability_df_grouped$label <- factor(stability_df_grouped$label, levels = stability_df_grouped$label)
  
  # Print the grouped data frame to debug
  print(stability_df_grouped)
  
  # Create the plot
  all_model_stability <- plot_ly(
    data = stability_df_grouped,
    x = ~mean,
    y = ~label,
    type = 'scatter',
    mode = 'markers',
    marker = list(
      size = ~12.5 * sqrt(pmax(pmin(frequency_percentage, max(stability_df_grouped$frequency_percentage)) / max(stability_df_grouped$frequency_percentage))),
      color = ~pmax(pmin(frequency_percentage, max(stability_df_grouped$frequency_percentage)), min(stability_df_grouped$frequency_percentage)),
      colorscale = 'Viridis',
      showscale = TRUE,
      colorbar = list(title = 'Frequency (%)')
    ),
    error_x = ~list(
      array = pooled_sem,
      color = 'black',
      thickness = 1
    )
  ) %>%
    layout(
      title = 'Variable Importance Plot',
      xaxis = list(title = 'Mean (Variable Importance)'),
      yaxis = list(title = 'Gene Symbol (Location Type)'),
      margin = list(b = 100)
    )
  
  return(list(all_model_stability, stability_df_grouped))
}

## Generate prediction probability boxplots:

### create dataframe for ROC by group
create_prediction_dataframe_meta <- function(results_list) {
  # Initialize an empty list to hold dataframes per result
  prediction_dfs <- list()
  
  # Iterate through each result in the results_list
  for (i in seq_along(results_list)) {
    result <- results_list[[i]]
    # Collect all sample indices from outer folds
    sample_indices <- unlist(result$repeated_fit$outer_folds)
    
    # Append sample indices to prediction_df
    temp_df <- data.frame(sample_index = sample_indices)
    
    # Add metadata columns by looking up respective indices in result$meta_ordered
    metadata <- result$meta_ordered[temp_df$sample_index, ]
    temp_df <- cbind(temp_df, metadata)
    
    # Collect all probabilities
    probabilities <- c()
    for (fold_pred in result$repeated_fit$outer_result) {
      fold_probabilities <- plogis(fold_pred$preds$predyp)
      probabilities <- c(probabilities, fold_probabilities)
    }
    
    # Add probabilities to the temp_df
    temp_df <- temp_df %>%
      mutate(probability = probabilities)
    
    # if average probability when condition == DISEASED or PM_positive is < 0.5, then flip the probabilities
    if (mean(temp_df$probability[temp_df$condition %in% c("DISEASED", "PM_positive")]) < 0.5) {
      temp_df$probability <- 1 - temp_df$probability
    }
    
    # Add the result index or model id
    temp_df$model_id <- i
    
    # Add temp_df to the list
    prediction_dfs[[i]] <- temp_df
  }
  
  return(prediction_dfs)
}

# Define the function
create_prediction_dataframe <- function(results_list) {
  # Initialize an empty dataframe for predictions
  prediction_df <- data.frame()
  
  # Iterate through each result in the results_list
  for (result in results_list) {
    # Collect all sample indices from outer folds
    sample_indices <- unlist(result$repeated_fit$outer_folds)
    
    # Append sample indices to prediction_df
    temp_df <- data.frame(sample_index = sample_indices)
    
    # Add metadata columns by looking up respective indices in result$meta_ordered
    metadata <- result$meta_ordered[temp_df$sample_index, ]
    temp_df <- cbind(temp_df, metadata)
    
    # Collect all probabilities
    probabilities <- c()
    for (fold_pred in result$repeated_fit$outer_result) {
      fold_probabilities <- plogis(fold_pred$preds$predyp)
      probabilities <- c(probabilities, fold_probabilities)
    }
    
    # Add probabilities to the temp_df
    temp_df <- temp_df %>%
      mutate(probability = probabilities)
    
    # Append temp_df to prediction_df
    prediction_df <- bind_rows(prediction_df, temp_df)
  }
  
  return(prediction_df)
}

library(plotly)
library(RColorBrewer)
library(pROC)

plot_average_roc_per_value <- function(roc_data_per_value, title, color_table) {
  library(plotly)
  library(RColorBrewer)
  library(pROC)
  
  p <- plot_ly()
  
  for (value_name in names(roc_data_per_value)) {
    roc_list <- roc_data_per_value[[value_name]]
    if (length(roc_list) == 0) next
    
    # Define common FPR points for interpolation
    fpr_common <- seq(0, 1, length.out = 100)
    
    # Interpolate TPR values at common FPR points for each ROC curve
    tpr_interpolated_list <- lapply(roc_list, function(roc_obj) {
      # Extract FPR and TPR from roc_obj
      fpr <- 1 - roc_obj$specificities
      tpr <- roc_obj$sensitivities
      # Interpolate TPR at common FPR points
      tpr_interpolated <- approx(fpr, tpr, xout = fpr_common, yleft = 0, yright = 1)$y
      return(tpr_interpolated)
    })
    
    # Convert list to matrix
    tpr_matrix <- do.call(cbind, tpr_interpolated_list)
    
    # Calculate mean and 95% CI of TPR at each FPR point using quantiles
    mean_tpr <- rowMeans(tpr_matrix, na.rm = TRUE)
    lower_tpr <- apply(tpr_matrix, 1, quantile, probs = 0.025, na.rm = TRUE)
    upper_tpr <- apply(tpr_matrix, 1, quantile, probs = 0.975, na.rm = TRUE)
    
    # Calculate AUCs for each ROC curve
    auc_values <- sapply(roc_list, function(roc_obj) {
      auc(roc_obj)
    })
    
    # Calculate mean AUC and 95% CI using quantiles
    mean_auc <- mean(auc_values, na.rm = TRUE)
    ci_lower_auc <- quantile(auc_values, 0.025, na.rm = TRUE)
    ci_upper_auc <- quantile(auc_values, 0.975, na.rm = TRUE)
    
    # Get the color for the current value_name from the provided color_table
    color <- color_table[[value_name]]
    if (is.null(color)) {
      # If no color provided for this value_name, use a default color
      color <- "#000000"  # Black as default
    }
    
    # Add the mean ROC curve with AUC and 95% CI in the legend
    p <- add_trace(p, 
                   x = fpr_common, 
                   y = mean_tpr, 
                   type = 'scatter', 
                   mode = 'lines',
                   line = list(color = color, width = 2),
                   name = paste0(
                     value_name, 
                     " (AUC = ", round(mean_auc, 3), 
                     " [95% CI: ", round(ci_lower_auc, 3), " - ", round(ci_upper_auc, 3), "])"
                   ),
                   hoverinfo = 'text',
                   text = paste(
                     "Value:", value_name,
                     "<br>FPR:", round(fpr_common, 3),
                     "<br>Mean TPR:", round(mean_tpr, 3),
                     "<br>95% CI TPR: ", round(lower_tpr, 3), " - ", round(upper_tpr, 3),
                     "<br>Mean AUC:", round(mean_auc, 3),
                     "<br>95% CI AUC:", round(ci_lower_auc, 3), "-", round(ci_upper_auc, 3)
                   )
    )
    
    # Add the ribbon for 95% confidence bands
    p <- add_ribbons(p, 
                     x = fpr_common, 
                     ymin = lower_tpr, 
                     ymax = upper_tpr, 
                     line = list(color = 'transparent'), 
                     fillcolor = adjustcolor(color, alpha.f = 0.2), 
                     showlegend = FALSE,
                     hoverinfo = 'skip'  # Avoid showing ribbon info on hover
    )
  }
  
  # Add the diagonal reference line
  p <- add_trace(p, 
                 x = c(0, 1), 
                 y = c(0, 1), 
                 type = 'scatter', 
                 mode = 'lines',
                 line = list(color = 'rgb(169, 169, 169)', dash = 'dash'), 
                 showlegend = FALSE)
  
  # Customize layout
  p <- p %>% layout(
    title = list(text = title, font = list(size = 20)),
    xaxis = list(
      title = 'False Positive Rate', 
      titlefont = list(size = 18), 
      showgrid = FALSE, 
      tickfont = list(size = 16)
    ),
    yaxis = list(
      title = 'True Positive Rate', 
      titlefont = list(size = 18), 
      showgrid = FALSE, 
      tickfont = list(size = 16),
      range = c(0, 1)
    ),
    legend = list(font = list(size = 12)),
    hovermode = 'closest',
    plot_bgcolor = 'rgba(0, 0, 0, 0)'
  )
  
  return(p)
}


### End of ROC by group

generate_boxplots_ggplot_pm <- function(data, color_tables) {
  # Load necessary libraries
  library(ggplot2)
  library(ggpubr)      # For stat_compare_means
  library(gridExtra)
  library(dplyr)
  
  # Ensure 'condition' is a factor
  data$condition <- as.factor(data$condition)
  
  # Use colors from color_tables for the 'condition' variable
  fill_colors <- color_tables$condition_table
  
  # Ensure that the names in fill_colors match the levels in data$condition
  # Reorder fill_colors to match levels of condition
  fill_colors <- fill_colors[levels(data$condition)]
  
  # Function to add sample counts to facet labels for a specific facet variable
  add_sample_counts <- function(data, facet_variable) {
    facet_counts <- data %>%
      group_by(!!sym(facet_variable)) %>%
      summarise(n = n()) %>%
      mutate(facet_label = paste0(!!sym(facet_variable), " (n=", n, ")"))
    
    # Update the data with new facet labels
    data <- data %>%
      left_join(facet_counts, by = facet_variable) %>%
      mutate(facet_label = factor(facet_label, levels = facet_counts$facet_label))
    
    return(data)
  }
  
  # Calculate total sample size for subtitle
  total_samples <- nrow(data)
  
  # Function to generate and customize faceted plots with median labels
  plot_facet <- function(data, facet_variable, title) {
    data <- add_sample_counts(data, facet_variable)
    
    p <- ggplot(data, aes(x = condition, y = probability, fill = condition)) +
      geom_boxplot(color = 'black', outlier.size = 0.5, alpha = 0.8) +
      # Add median labels
      stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)),
                   position = position_dodge(width = 0.75), 
                   vjust = 0.5, size = 3, color = "white") +
      facet_wrap(~ facet_label) +
      scale_fill_manual(values = fill_colors) +
      labs(
        title = title,
        subtitle = paste0("Total Samples: n=", total_samples),
        y = "Probability",
        x = "Condition"  # You can remove or modify this if you want to hide the axis title
      ) +
      # Add p-value annotations using stat_compare_means
      stat_compare_means(
        method = "wilcox.test",
        label = "p.format",
        label.y = 0.95,  # Position labels at 95% of y-axis
        hide.ns = TRUE   # Hide non-significant p-values
      ) +
      scale_y_continuous(
        limits = c(0, 1),              # Set y-axis limits from 0 to 1
        breaks = seq(0, 1, by = 0.2),  # Define y-axis ticks at every 0.2
        expand = expansion(mult = c(0, 0))  # Remove extra space
      ) +
      theme_minimal() +
      theme(
        strip.background = element_rect(fill = "grey80"),
        strip.text = element_text(margin = margin(t = 10, b = 10)),
        axis.text.x = element_blank(),      # Hide x-axis text labels
        axis.title.x = element_blank(),     # Optional: Hide x-axis title
        plot.margin = margin(b = 20),
        legend.position = "none"  # Remove legend if not needed
      )
    
    return(p)
  }
  
  # Function to generate and customize the main plot (no faceting) with median labels
  plot_main <- function(data, title) {
    p <- ggplot(data, aes(x = condition, y = probability, fill = condition)) +
      geom_boxplot(color = 'black', outlier.size = 0.5, alpha = 0.8) +
      # Add median labels
      stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)),
                   position = position_dodge(width = 0.75), 
                   vjust = 0.5, size = 3, color = "white") +
      scale_fill_manual(values = fill_colors) +
      labs(
        title = title,
        subtitle = paste0("Total Samples: n=", total_samples),
        y = "Probability",
        x = "Condition"  # You can remove or modify this if you want to hide the axis title
      ) +
      # Add p-value annotations using stat_compare_means
      stat_compare_means(
        method = "wilcox.test",
        label = "p.format",
        label.y = 0.95,  # Position labels at 95% of y-axis
        hide.ns = TRUE   # Hide non-significant p-values
      ) +
      scale_y_continuous(
        limits = c(0, 1),              # Set y-axis limits from 0 to 1
        breaks = seq(0, 1, by = 0.2),  # Define y-axis ticks at every 0.2
        expand = expansion(mult = c(0, 0))  # Remove extra space
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),      # Hide x-axis text labels
        axis.title.x = element_blank(),     # Optional: Hide x-axis title
        plot.margin = margin(b = 20),
        legend.position = "none"  # Remove legend if not needed
      )
    
    return(p)
  }
  
  # Generate the main plot (no faceting)
  p1 <- plot_main(
    data = data,
    title = "Peritoneal Metastasis Status"
  )
  
  # Generate faceted plots with median labels
  p2 <- plot_facet(
    data = data,
    facet_variable = "primary_present",
    title = "Peritoneal Metastasis Status and Presence of Primary"
  )
  
  p3 <- plot_facet(
    data = data,
    facet_variable = "primary_site",
    title = "Peritoneal Metastasis Status and Primary Site"
  )
  
  # Combine the plots into a single plot with facets using gridExtra
  combined_plot <- gridExtra::grid.arrange(
    p1, p2, p3, 
    ncol = 1, 
    heights = c(1, 1, 1)
  )
  
  # Return the combined plot
  return(combined_plot)
}



generate_boxplots_ggplot_healthy <- function(data, color_tables) {
  # Load necessary libraries
  library(ggplot2)
  library(ggpubr)      # For stat_compare_means
  library(gridExtra)
  library(dplyr)
  
  # Ensure variables are factors
  data <- data %>%
    mutate(
      condition = as.factor(condition),
      peritoneal_mets = as.factor(peritoneal_mets),
      primary_present = as.factor(primary_present),
      primary_site = as.factor(primary_site),
      age_cat = as.factor(age_cat),
      sex = as.factor(sex),
      batch = as.factor(batch)
    )
  
  # Define pairwise comparisons within 'condition'
  comparisons <- combn(levels(data$condition), 2, simplify = FALSE)
  
  # Set y-axis limits from 0 to 1
  y_limit <- 1
  
  # Calculate total sample size for subtitle
  total_samples <- nrow(data)
  
  # Function to generate and customize individual boxplots with centered, white median annotations
  generate_plot <- function(data, fill_variable, fill_colors, title) {
    p <- ggplot(data, aes_string(x = "condition", y = "probability", fill = fill_variable)) +
      geom_boxplot(color = 'black', outlier.size = 0.5, alpha = 0.8) +
      # Add centered, white median labels
      stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)),position = position_dodge(width = 0.8),
                   vjust = 0.5, size = 3, color = "white") +
      scale_fill_manual(values = fill_colors) +
      labs(
        title = title,
        subtitle = paste0("Total Samples: n=", total_samples),
        y = "Probability",
        x = "Condition"
      ) +
      # Add p-value annotations using stat_compare_means
      stat_compare_means(
        method = "wilcox.test",
        comparisons = comparisons,
        label = "p.format",
        label.y = 0.95,  # Position labels at 95% of y-axis
        hide.ns = TRUE   # Hide non-significant p-values
      ) +
      scale_y_continuous(
        limits = c(0, y_limit),
        breaks = seq(0, 1, by = 0.2),
        expand = expansion(mult = c(0, 0))  # Remove extra space
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(b = 20),
        legend.position = "none"  # Remove legend if not needed
      )
    
    return(p)
  }
  
  # Generate each plot without faceting
  # Plot 1: Condition
  p1 <- generate_plot(
    data = data,
    fill_variable = "condition",
    fill_colors = color_tables$condition_table,
    title = "Condition"
  )
  
  # Plot 2: Peritoneal Metastasis Status
  p2 <- generate_plot(
    data = data,
    fill_variable = "peritoneal_mets",
    fill_colors = color_tables$pm_mets_table,
    title = "Peritoneal Metastasis Status"
  )
  
  # Plot 3: Primary Tumor Presence
  p3 <- generate_plot(
    data = data,
    fill_variable = "primary_present",
    fill_colors = color_tables$primary_present_table,
    title = "Primary Tumor Presence"
  )
  
  # Plot 4: Primary Site
  p4 <- generate_plot(
    data = data,
    fill_variable = "primary_site",
    fill_colors = color_tables$primary_site_table,
    title = "Primary Site"
  )
  
  # Plot 5: Age Category
  p5 <- generate_plot(
    data = data,
    fill_variable = "age_cat",
    fill_colors = color_tables$age_cat_table,
    title = "Age Category"
  )
  
  # Plot 6: Sex
  p6 <- generate_plot(
    data = data,
    fill_variable = "sex",
    fill_colors = color_tables$sex_table,
    title = "Sex"
  )
  
  # Plot 7: Batch
  p7 <- generate_plot(
    data = data,
    fill_variable = "batch",
    fill_colors = color_tables$batch_table,
    title = "Batch"
  )
  
  # Combine all plots into a grid layout
  combined_plot <- gridExtra::grid.arrange(
    p1, p2, p3, p4, p5, p6, p7,
    ncol = 2,
    top = "Combined Probability Boxplots",
    bottom = "Y-Axis Range: 0 to 1"
  )
  
  # Return the combined plot
  return(combined_plot)
}



### Functions for generating heatmaps
generate_heatmaps <- function(training_data, 
                              validation_data = NULL, 
                              holdout_data = NULL, 
                              annotated_features, 
                              meta_training, 
                              meta_validation = NULL, 
                              meta_holdout = NULL, 
                              stability_df_grouped, 
                              color_tables = NULL) {
  # Load the ComplexHeatmap package
  library(ComplexHeatmap)
  library(circlize)  # Also load circlize for colorRamp2 function
  
  # drop columns in training data that are not in annotated_features$GeneSymbol
  training_data_top <- training_data[rownames(training_data) %in% annotated_features$feature, ]
  
  preprocess_data <- function(data, features) {
    if (is.null(data)) return(NULL)
    
    # Log2 transform the normalized counts
    log2_data <- data #log2(data + 1)  # Add 1 to avoid log(0)
    
    # Optionally, you might want to center each gene/feature
    centered_data <- t(scale(t(log2_data), scale = TRUE, center = TRUE))
    
    # Subset to features of interest if necessary
    if (!is.null(features)) {
      centered_data <- centered_data[rownames(centered_data) %in% features$feature, ]
    }
    
    return(centered_data)
  }
  
  # Preprocess the data
  counts_data_training_top_t <- preprocess_data(training_data_top, annotated_features)
  counts_data_validation_top_t <- preprocess_data(validation_data, annotated_features)
  counts_data_holdout_top_t <- preprocess_data(holdout_data, annotated_features)
  
  # Check if training data is NULL
  if (is.null(counts_data_training_top_t) || is.null(meta_training)) {
    stop("Training data and meta_training must not be NULL")
  }
  
  # Helper function to select metadata columns
  select_meta_columns <- function(meta) {
    if (is.null(meta)) return(NULL)
    meta %>%
      dplyr::select(condition, primary_site, primary_present, peritoneal_mets, batch, sex, race, age_cat)
  }
  
  # Select relevant metadata columns
  annotation_training <- select_meta_columns(meta_training)
  annotation_validation <- select_meta_columns(meta_validation)
  annotation_holdout <- select_meta_columns(meta_holdout)
  
  # Helper function to create heatmap annotation
  create_annotation <- function(annotation, condition_table, primary_site_table, primary_present_table, peritoneal_mets_table, batch_table,sex_table,race_table,age_cat_table) {
    if (is.null(annotation)) return(NULL)
    HeatmapAnnotation(
      condition = annotation$condition,
      primary_present = annotation$primary_present,
      primary_site = annotation$primary_site,
      sex = annotation$sex,
      batch = annotation$batch,
      age_cat = annotation$age_cat,
      race = annotation$race,
      col = list(
        condition = color_tables$condition_table,
        primary_site = color_tables$primary_site_table,
        sex = color_tables$sex_table,
        primary_present = color_tables$primary_present_table,
        batch = color_tables$batch_table,
        age_cat = color_tables$age_cat_table,
        race = color_tables$race_table
      )
    )
  }
  
  # Create heatmap annotations
  ha_training <- create_annotation(annotation_training, condition_table, primary_site_table, primary_present_table, peritoneal_mets_table, batch_table,sex_table,race_table,age_cat_table)
  ha_validation <- create_annotation(annotation_validation, condition_table, primary_site_table, primary_present_table, peritoneal_mets_table, batch_table,sex_table,race_table,age_cat_table)
  ha_holdout <- create_annotation(annotation_holdout, condition_table, primary_site_table, primary_present_table, peritoneal_mets_table, batch_table,sex_table,race_table,age_cat_table)
  
  robust_dist = function(x, y) {
    qx = quantile(x, c(0.05, 0.95))
    qy = quantile(y, c(0.05, 0.95))
    l = x > qx[1] & x < qx[2] & y > qy[1] & y < qy[2]
    x = x[l]
    y = y[l]
    sqrt(sum((x - y)^2))
  }
  
  # Helper function to create heatmaps
  create_heatmap <- function(data, annotation, stability_df_grouped) {
    if (is.null(data)) return(NULL)
    
    # Create a named vector for row labels
    row_labels <- stability_df_grouped$all_symbols[match(rownames(data), stability_df_grouped$feature)]
    names(row_labels) <- rownames(data)
    
    Heatmap(data, 
            top_annotation = annotation,
            col = colorRamp2(c(min(data, na.rm = TRUE), min(data, na.rm = TRUE)/4, min(data, na.rm = TRUE)/8, 0, max(data, na.rm = TRUE)/8, max(data, na.rm = TRUE)/4, max(data, na.rm = TRUE)), 
                             c("#4676b5", "#82b1d3", "#dbeff6", "white", "#fee395", "#fc9961", "#d73027")),
            show_row_names = TRUE,
            row_labels = row_labels,
            show_column_names = FALSE,
            clustering_distance_rows = "euclidean", #"euclidean", #robust_dist,#
            clustering_distance_columns = "euclidean", #"pearson", #robust_dist #"euclidean", # euclidean pearson
            clustering_method_columns = "ward.D2"
            )
  }
  
  # Create heatmaps
  heatmap_tra <- create_heatmap(counts_data_training_top_t, ha_training, annotated_features)
  heatmap_val <- create_heatmap(counts_data_validation_top_t, ha_validation, annotated_features)
  heatmap_hol <- create_heatmap(counts_data_holdout_top_t, ha_holdout, annotated_features)
  
  # Combine heatmaps, ignoring NULLs
  heatmaps <- list(heatmap_tra, heatmap_val, heatmap_hol)
  combined_heatmap <- Reduce(`+`, heatmaps[!sapply(heatmaps, is.null)])
  
  return(combined_heatmap)
}

### Generate boxplot by gene
generate_boxplots_by_gene <- function(training_data, meta_training, annotated_features,color_tables) {
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(ggpubr)
  
  # Subset training data to include only the annotated features
  training_data_top <- training_data[rownames(training_data) %in% annotated_features$feature, ]
  
  # Transpose the data and convert to long format
  long_data <- as.data.frame(t(training_data_top)) %>%
    tibble::rownames_to_column("sample") %>%
    pivot_longer(cols = -sample, names_to = "gene", values_to = "expression")
  
  # Add condition information
  meta_training$sample <- rownames(meta_training)
  long_data <- long_data %>%
    left_join(meta_training %>% dplyr::select(sample, condition), by = "sample")
  
  # Add gene symbols and mean values
  gene_info <- annotated_features %>%
    dplyr::select(feature, all_symbols, mean) %>%
    dplyr::rename(gene = feature, gene_symbol = all_symbols)
  
  long_data <- long_data %>%
    left_join(gene_info, by = "gene")
  
  # Sort genes based on mean values
  gene_order <- gene_info %>%
    arrange(desc(mean)) %>%
    pull(gene_symbol)
  
  # Convert gene_symbol to a factor with levels ordered by mean
  long_data$gene_symbol <- factor(long_data$gene_symbol, levels = gene_order)
  
  # Ensure condition is a factor
  long_data$condition <- factor(long_data$condition)
  
  # Create a function to generate plot for a subset of data
  create_plot <- function(data, title) {
    # Ensure condition is a factor
    data$condition <- factor(data$condition)
    condition_levels <- levels(data$condition)
    
    # Map conditions to colors
    condition_colors <- sapply(condition_levels, function(cond) {
      if(cond %in% c("DISEASED", "PM_positive")) {
        color_tables$condition_table[["DISEASED"]]  # Red color
      } else if (cond == "HEALTHY") {
        color_tables$condition_table[["HEALTHY"]] # Grey color
      } else if (cond == "PM_negative") {
        color_tables$condition_table[["PM_negative"]]  # Blue color
      } else {
        "black"  # Default color for any other conditions
      }
    })
    names(condition_colors) <- condition_levels
    
    # Calculate whisker ends for each gene and condition
    whisker_ends <- data %>%
      group_by(gene_symbol, condition) %>%
      summarize(
        q1 = quantile(expression, 0.25),
        q3 = quantile(expression, 0.75),
        iqr = q3 - q1,
        whisker_top = min(max(expression), q3 + 1.5 * iqr),
        .groups = "drop"
      ) %>%
      group_by(gene_symbol) %>%
      summarize(max_whisker = max(whisker_top), .groups = "drop")
    
    # Calculate the maximum y-range
    max_y <- max(whisker_ends$max_whisker)
    
    # Create the base plot
    p <- ggplot(data, aes(x = gene_symbol, y = expression, fill = condition)) +
      geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, outlier.shape = NA) +
      #geom_jitter(aes(color = condition), 
      #            position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2),
      #            size = 0.5, alpha = 0.3,
      #            stroke = 0.3) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "top") +
      labs(x = "Gene", y = "Expression", fill = "Condition", color = "Condition", title = title) +
      scale_fill_manual(values = condition_colors) +
      scale_color_manual(values = condition_colors) +
      guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
      scale_y_continuous(limits = c(0, max_y))  # Set the y-axis limits
    
    return(p)
  }
  
  # Create faceted plot function
  create_faceted_plot <- function(data, title, is_positive = TRUE) {
    library(ggplot2)
    library(dplyr)
    library(ggpubr)
    library(gridExtra)
    
    # Dataframe to store statistical results
    stats_results_list <- list()
    
    if (nrow(data) == 0) {
      stop("Error: Input data is empty. Please check your data source.")
    }
    
    # Ensure 'condition' is a factor
    data$condition <- as.factor(data$condition)
    condition_levels <- levels(data$condition)
    
    # Map conditions to colors
    condition_colors <- sapply(condition_levels, function(cond) {
      if(cond %in% c("DISEASED", "PM_positive")) {
        color_tables$condition_table[["DISEASED"]]  # Red color
      } else if (cond == "HEALTHY") {
        color_tables$condition_table[["HEALTHY"]] # Grey color
      } else if (cond == "PM_negative") {
        color_tables$condition_table[["PM_negative"]]  # Blue color
      } else {
        "black"  # Default color for any other conditions
      }
    })
    names(condition_colors) <- condition_levels
    
    n_conditions <- n_distinct(data$condition)
    
    # Calculate Q3 + 1.5*IQR for each gene
    label_y_data <- data %>%
      group_by(gene) %>%
      summarise(label_y = quantile(expression, 0.75) + 1.5 * IQR(expression)) %>%
      ungroup()
    
    # Get unique genes from the data
    genes <- unique(data$gene)
    
    # Create a list to store individual plots
    plot_list <- list()
    
    # Create subplots for each gene
    for (i in seq_along(genes)) {
      current_gene <- genes[i]
      
      # Subset data for the current gene
      gene_data <- data %>% filter(gene == current_gene)
      
      if (nrow(gene_data) == 0) {
        warning(paste("No data found for gene:", current_gene))
        next
      }
      
      # Calculate max y value
      max_y <- max(gene_data$expression, na.rm = TRUE)
      
      # Get label_y value for the current gene
      current_label_y <- label_y_data %>% 
        filter(gene == current_gene) %>% 
        pull(label_y)
      
      if (length(current_label_y) == 0) {
        current_label_y <- max_y
      }
      
      # Calculate the highest label position
      if (n_conditions >= 2) {
        comparisons <- combn(levels(gene_data$condition), 2, simplify = FALSE)
        step_increase <- 0.05
        highest_label_y <- current_label_y + (step_increase * current_label_y * (length(comparisons) - 1))
      } else {
        highest_label_y <- current_label_y
      }
      
      # Set ymax to be slightly above the highest label position
      ymax <- highest_label_y * 1.15
      
      # Get the gene_symbol for the current gene
      gene_symbol <- unique(gene_data$gene_symbol)
      
      # Calculate mean values for each condition
      mean_values <- gene_data %>%
        group_by(condition) %>%
        summarise(mean_expression = mean(expression, na.rm = TRUE)) %>%
        ungroup()
      
      # Perform statistical tests and store results
      if (n_conditions >= 2) {
        test_results <- compare_means(expression ~ condition, data = gene_data,
                                      method = "wilcox.test", p.adjust.method = "BH")
        test_results$gene <- gene_symbol  # Add gene information
        test_results$mean_group1 <- mean_values$mean_expression[mean_values$condition == test_results$group1]
        test_results$mean_group2 <- mean_values$mean_expression[mean_values$condition == test_results$group2]
        
        stats_results_list[[i]] <- test_results
      }
      
      # Create the plot for the current gene
      p <- ggplot(gene_data, aes(x = condition, y = expression, fill = condition)) +
        geom_boxplot(width = 0.7, outlier.shape = NA) +
        #geom_jitter(aes(color = condition), position = position_jitter(width = 0.2), size = 0.5, alpha = 0.3, stroke = 0.3) +
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = "none",
              plot.title = element_text(face = "bold")) +
        labs(x = NULL, y = "5hmC Normalized Counts", title = gene_symbol) +
        scale_fill_manual(values = condition_colors) +
        scale_color_manual(values = condition_colors) +
        coord_cartesian(ylim = c(NA, ymax))  # Set the upper limit of y-axis to ymax
      
      # Add mean value annotations
      p <- p + geom_text(data = mean_values, 
                         aes(x = condition, y = -Inf, vjust = 0,label = sprintf("%.1f", mean_expression)),
                         position = position_dodge(width = 0.8), size = 3, fontface = "bold")
      
      # Add stat_compare_means only if there are at least two conditions
      if (n_conditions >= 2) {
        p <- p + stat_compare_means(
          comparisons = comparisons,
          label = "p = {p.adj}",
          method = "wilcox.test",
          step.increase = step_increase,
          label.y = current_label_y,
          tip.length = 0
        )
      }
      
      plot_list[[i]] <- p
    }
    
    if (length(plot_list) == 0) {
      stop("Error: No plots were created. Please check your data and gene symbols.")
    }
    
    # Combine statistical results into a dataframe
    if (length(stats_results_list) > 0) {
      stats_results <- do.call(rbind, stats_results_list)
    } else {
      stats_results <- data.frame()
    }
    
    # Combine all plots
    combined_plot <- ggpubr::ggarrange(plotlist = plot_list, ncol = 4, nrow = ceiling(length(plot_list) / 4))
    
    # Create a separate legend
    legend_plot <- ggplot(data, aes(x = condition, y = expression, fill = condition)) +
      geom_boxplot() +
      scale_fill_manual(values = condition_colors) +
      theme_void() +
      theme(legend.position = "bottom")
    
    legend <- get_legend(legend_plot)
    
    # Add overall title and legend
    final_plot <- ggpubr::annotate_figure(
      combined_plot,
      top = ggpubr::text_grob(title, face = "bold", size = 14),
      bottom = legend
    )
    
    # Return both the plot and the statistical results
    return(list(plot = final_plot, stats_results = stats_results))
  }
  
  # Split data based on mean values
  negative_mean <- long_data %>% filter(mean < 0)
  positive_mean <- long_data %>% filter(mean >= 0)
  
  # Create plots and collect statistical results
  negative_results <- create_faceted_plot(negative_mean, "Genes with Negative Mean Expression")
  positive_results <- create_faceted_plot(positive_mean, "Genes with Positive Mean Expression")
  
  plot_negative <- create_plot(negative_mean, "Genes with Negative Mean Expression")
  plot_positive <- create_plot(positive_mean, "Genes with Positive Mean Expression")
  
  # Combine statistical results
  all_stats_results <- bind_rows(negative_results$stats_results, positive_results$stats_results)
  
  # Return a list of plots and statistical results
  return(list(
    negative = plot_negative, 
    positive = plot_positive, 
    faceted_negative = negative_results$plot, 
    faceted_positive = positive_results$plot,
    stats_results = all_stats_results
  ))
}

library(plotly)
library(dplyr)
library(pROC)
library(caret)
library(RColorBrewer)

plot_roc_curves_per_feature_nestedcv_logreg_ARCHIVE <- function(normcounts_data, meta_data, stability_df_grouped, title_prefix, n_outer_folds = 5) {
  # Ensure condition is a factor with PM_negative as the reference level
  meta_data$condition <- factor(meta_data$condition)#, levels = c("PM_negative", "PM_positive"))
  
  # Split features based on mean coefficient sign
  negative_features <- stability_df_grouped$feature[stability_df_grouped$mean < 0]
  positive_features <- stability_df_grouped$feature[stability_df_grouped$mean > 0]
  
  # Function to perform nested CV for a single feature
  nested_cv_for_feature <- function(feature) {
    feature_data <- normcounts_data[feature, ]
    model_data <- data.frame(
      condition = meta_data$condition,
      feature_value = as.numeric(feature_data)
    )
    
    # Outer fold
    outer_folds <- createFolds(model_data$condition, k = n_outer_folds)
    outer_results <- lapply(outer_folds, function(test_index) {
      train <- model_data[-test_index, ]
      test <- model_data[test_index, ]
      
      # Train model
      model <- glm(condition ~ feature_value, data = train, family = "binomial")
      
      # Predict on test set
      predictions <- predict(model, newdata = test, type = "response")
      list(true = test$condition, pred = predictions)
    })
    
    # Combine results from all outer folds
    true_labels <- do.call(c, lapply(outer_results, `[[`, "true"))
    pred_probs <- do.call(c, lapply(outer_results, `[[`, "pred"))
    
    # Create ROC object
    roc_obj <- roc(true_labels, pred_probs)
    
    return(roc_obj)
  }
  
  # Function to create plot for a set of features
  create_plot <- function(features, title, color_palette) {
    p <- plot_ly()
    
    # Generate color palette
    num_colors <- min(9, length(features))  # RColorBrewer palettes have a maximum of 9 colors
    colors <- colorRampPalette(brewer.pal(num_colors, color_palette))(length(features))
    
    for (i in seq_along(features)) {
      feature <- features[i]
      roc_obj <- nested_cv_for_feature(feature)
      
      # Get the label for this feature
      feature_label <- stability_df_grouped$label[stability_df_grouped$feature == feature]
      
      # Add trace for this feature
      p <- add_trace(p, x = 1 - roc_obj$specificities, y = roc_obj$sensitivities, 
                     type = 'scatter', mode = 'lines',
                     name = paste(feature_label, "AUC:", round(auc(roc_obj), 3)),
                     line = list(color = colors[i]),
                     hoverinfo = 'text',
                     text = paste("Feature:", feature_label, 
                                  "<br>FPR:", round(1 - roc_obj$specificities, 2), 
                                  "<br>TPR:", round(roc_obj$sensitivities, 2)))
    }
    
    # Add the diagonal reference line
    p <- add_trace(p, x = c(0, 1), y = c(0, 1), type = 'scatter', mode = 'lines',
                   line = list(color = 'rgb(0, 0, 0)', dash = 'dash'), 
                   showlegend = FALSE)
    
    # Add layout settings
    p <- p %>% layout(
      title = list(text = title, font = list(size = 20)),
      xaxis = list(title = 'False Positive Rate', titlefont = list(size = 18), 
                   showgrid = FALSE, tickfont = list(size = 16)),
      yaxis = list(title = 'True Positive Rate', titlefont = list(size = 18), 
                   showgrid = FALSE, tickfont = list(size = 16),
                   range = c(0, 1)),
      legend = list(x = 1.05, y = 1, xanchor = 'left', yanchor = 'top'),
      hovermode = 'closest'
    )
    
    return(p)
  }
  
  # Create plots
  negative_plot <- create_plot(negative_features, paste(title_prefix, "- Negative Mean Coefficients"), "PuBu")
  positive_plot <- create_plot(positive_features, paste(title_prefix, "- Positive Mean Coefficients"), "YlOrRd")
  
  return(list(negative_plot = negative_plot, positive_plot = positive_plot))
}

library(plotly)
library(RColorBrewer)
library(caret)
library(pROC)

plot_roc_curves_per_feature_nestedcv_logreg <- function(normcounts_data, meta_data, stability_df_grouped, title_prefix, n_outer_folds = 5) {
  # Ensure condition is a factor with the correct reference level
  meta_data$condition <- factor(meta_data$condition, levels = c("PM_negative", "PM_positive"))
  
  # Split features based on mean coefficient sign
  negative_features <- stability_df_grouped$feature[stability_df_grouped$mean < 0]
  positive_features <- stability_df_grouped$feature[stability_df_grouped$mean > 0]
  
  # Function to perform nested CV for a single feature
  nested_cv_for_feature <- function(feature) {
    # Initialize lists to store ROC curves and AUCs
    roc_list <- list()
    auc_list <- numeric()
    
    for (seed in 1:50) {
      set.seed(seed)
      
      feature_data <- normcounts_data[feature, ]
      model_data <- data.frame(
        condition = meta_data$condition,
        feature_value = as.numeric(feature_data)
      )
      
      # Ensure there are both classes present
      if(length(unique(model_data$condition)) < 2){
        next  # Skip this seed if only one class is present
      }
      
      # Outer fold
      outer_folds <- createFolds(model_data$condition, k = n_outer_folds, returnTrain = FALSE)
      
      outer_results <- tryCatch({
        lapply(outer_folds, function(test_index) {
          train <- model_data[-test_index, ]
          test <- model_data[test_index, ]
          
          # Ensure both classes are present in the training set
          if(length(unique(train$condition)) < 2){
            return(NULL)  # Skip this fold
          }
          
          # Train model
          model <- glm(condition ~ feature_value, data = train, family = "binomial")
          
          # Predict on test set
          predictions <- predict(model, newdata = test, type = "response")
          
          list(true = test$condition, pred = predictions)
        })
      }, error = function(e) {
        return(NULL)  # Skip this seed if an error occurs
      })
      
      # Remove NULL results (folds that were skipped)
      outer_results <- outer_results[!sapply(outer_results, is.null)]
      
      # Continue only if there are valid results
      if(length(outer_results) == 0){
        next  # Skip to the next seed
      }
      
      # Collect results
      true_labels <- do.call(c, lapply(outer_results, `[[`, "true"))
      pred_probs <- do.call(c, lapply(outer_results, `[[`, "pred"))
      
      # Check if both classes are present in true_labels
      if(length(unique(true_labels)) < 2){
        next  # Skip this seed
      }
      
      # Create ROC object
      roc_obj <- tryCatch({
        roc(true_labels, pred_probs, levels = c("PM_negative", "PM_positive"), direction = "<")
      }, error = function(e) {
        return(NULL)  # Skip this seed if ROC computation fails
      })
      
      if(is.null(roc_obj)){
        next  # Skip this seed
      }
      
      roc_list[[length(roc_list) + 1]] <- roc_obj
      auc_list <- c(auc_list, auc(roc_obj))
    }
    
    # Check if we have at least one valid ROC curve
    if(length(roc_list) == 0){
      return(NULL)  # Skip this feature if no valid ROC curves are available
    }
    
    # Define a sequence of FPR values
    fpr_values <- seq(0, 1, length.out = 100)
    
    # Initialize a matrix to store interpolated TPR values
    tpr_matrix <- matrix(NA, nrow = length(fpr_values), ncol = length(roc_list))

# Populate the matrix with interpolated TPR values
for(i in seq_along(roc_list)) {
  roc_obj <- roc_list[[i]]
  # Use linear interpolation to get TPR at specified FPR values
  tpr_interp <- tryCatch({
    # Ensure that the ROC object covers the entire FPR range
    # 'coords' can be used alternatively, but 'approx' provides more control
    approx(x = 1 - roc_obj$specificities, y = roc_obj$sensitivities, xout = fpr_values, rule = 2)$y
  }, error = function(e) {
    rep(NA, length(fpr_values))  # If interpolation fails, fill with NA
  })
  tpr_matrix[, i] <- tpr_interp
}

# Compute mean TPR, ignoring NA values
mean_tpr <- rowMeans(tpr_matrix, na.rm = TRUE)

# Compute mean AUC
mean_auc <- mean(auc_list, na.rm = TRUE)

# Compute 95% Confidence Interval for AUC
ci_lower <- quantile(auc_list, 0.025, na.rm = TRUE)
ci_upper <- quantile(auc_list, 0.975, na.rm = TRUE)

# Create mean ROC object with CI
mean_roc_obj <- list(
  specificities = 1 - fpr_values,  # This will represent FPR
  sensitivities = mean_tpr,
  auc = mean_auc,
  ci_lower = ci_lower,
  ci_upper = ci_upper
)
class(mean_roc_obj) <- "roc"

return(mean_roc_obj)
  }
  
  # Function to create plot for a set of features
  create_plot <- function(features, title, color_palette) {
    p <- plot_ly()
    
    # Generate color palette
    num_colors <- min(9, length(features))  # RColorBrewer palettes have a maximum of 9 colors
    colors <- colorRampPalette(brewer.pal(num_colors, color_palette))(length(features))
    
    for (i in seq_along(features)) {
      feature <- features[i]
      roc_obj <- nested_cv_for_feature(feature)
      
      # Skip plotting if no valid ROC was computed
      if(is.null(roc_obj)){
        next
      }
      
      # Get the label for this feature
      feature_label <- stability_df_grouped$label[stability_df_grouped$feature == feature]
      
      # Add trace for the mean ROC curve of this feature
      p <- add_trace(p,
                     x = 1 - roc_obj$specificities,  # FPR on X-axis
                     y = roc_obj$sensitivities,      # TPR on Y-axis
                     type = 'scatter',
                     mode = 'lines',
                     name = paste0(
                       feature_label, " Mean AUC: ", round(roc_obj$auc, 3),
                       " (95% CI: ", round(roc_obj$ci_lower, 3), " - ", round(roc_obj$ci_upper, 3), ")"
                     ),
                     line = list(color = colors[i]),
                     hoverinfo = 'text',
                     text = paste(
                       "Feature:", feature_label,
                       "<br>FPR:", round(1 - roc_obj$specificities, 2),
                       "<br>TPR:", round(roc_obj$sensitivities, 2),
                       "<br>AUC:", round(roc_obj$auc, 3),
                       "<br>95% CI:", round(roc_obj$ci_lower, 3), "-", round(roc_obj$ci_upper, 3)
                     ))
    }
    
    # Add the diagonal reference line (FPR vs. TPR)
    p <- add_trace(p, x = c(0, 1), y = c(0, 1), type = 'scatter', mode = 'lines',
                   line = list(color = 'rgb(0, 0, 0)', dash = 'dash'), 
                   showlegend = FALSE)
    
    # Add layout settings
    p <- p %>% layout(
      title = list(text = title, font = list(size = 20)),
      xaxis = list(title = 'False Positive Rate', titlefont = list(size = 18), 
                   showgrid = FALSE, tickfont = list(size = 16)),
      yaxis = list(title = 'True Positive Rate', titlefont = list(size = 18), 
                   showgrid = FALSE, tickfont = list(size = 16),
                   range = c(0, 1)),
      legend = list(x = 1.05, y = 1, xanchor = 'left', yanchor = 'top'),
      hovermode = 'closest'
    )
    
    return(p)
  }
  
  # Create plots
  negative_plot <- create_plot(negative_features, paste(title_prefix, "- Negative Mean Coefficients"), "PuBu")
  positive_plot <- create_plot(positive_features, paste(title_prefix, "- Positive Mean Coefficients"), "YlOrRd")
  
  return(list(negative_plot = negative_plot, positive_plot = positive_plot))
}


#### Processing enhancer annotated table
library(dplyr)
library(tidyr)


process_csv <- function(file_path = NULL, df = NULL) {
  # Read the CSV file if a file path is provided, otherwise use the provided dataframe
  if (!is.null(file_path)) {
    df <- read_csv(file_path, col_names = TRUE, show_col_types = FALSE)
    input_is_file <- TRUE
  } else if (is.null(df)) {
    stop("Either file_path or df must be provided")
  } else {
    input_is_file <- FALSE
  }
  
  # Process the dataframe
  result <- df %>%
    # First, rename only the columns that don't cause conflicts
    rename(
      mean_coeff = mean,
      gene_symbol = annot.symbol,
      location_type = locationType,
      enhancers = label,
      enhancer_targets = symbols_gh_full
    ) %>%
    # Then, handle the conflicting columns separately and ensure enhancers is character type
    mutate(
      frequency = round(frequency_percentage, 1),  # Round to 1 decimal place
      mean_coeff = round(mean_coeff, 4),  # Round to 4 decimal places
      strand = as.character(annot.strand),  # Create a new column instead of renaming
      location_type = ifelse(location_type == "Enhancer", "Intergenic", location_type),
      enhancers = as.character(enhancers),  # Ensure enhancers is character type
      enhancers = sapply(strsplit(enhancers, ", "), function(x) {
        gh_values <- grep("^GH\\d+", x, value = TRUE)
        if (length(gh_values) > 0) paste(gh_values, collapse = ", ") else NA_character_
      }),
      enhancer_targets = as.character(enhancer_targets),  # Ensure enhancer_targets is character type
      enhancer_targets = sapply(strsplit(enhancer_targets, ", "), function(x) {
        non_na_values <- x[x != "NA" & !is.na(x)]
        paste(non_na_values[!grepl("^chr(\\d+|X)_", non_na_values)], collapse = ", ")
      }),
      gene_symbol = ifelse(location_type == "Intergenic", "", gene_symbol)  # Set gene_symbol to empty string if location_type is Intergenic
    ) %>%
    filter(!is.na(enhancers)) %>%
    group_by(chr, start, end) %>%
    summarise(
      frequency = dplyr::first(frequency),  # Explicitly use dplyr::first()
      mean_coeff = dplyr::first(mean_coeff),
      strand = dplyr::first(strand),
      gene_symbol = dplyr::first(gene_symbol),
      location_type = dplyr::first(location_type),
      enhancers = dplyr::first(enhancers),
      enhancer_targets = paste(unique(unlist(strsplit(enhancer_targets, ", "))), collapse = ", "),
      .groups = "drop"
    ) %>%
    dplyr::select(frequency, mean_coeff, chr, start, end, strand, gene_symbol, location_type, enhancers, enhancer_targets)
  
  # If input was a file, save the result to the same folder
  if (input_is_file) {
    # Extract directory and file name
    dir_path <- dirname(file_path)
    file_name <- basename(file_path)
    
    # Create new file name
    new_file_name <- paste0("processed_", file_name)
    
    # Combine directory and new file name
    output_path <- file.path(dir_path, new_file_name)
    
    # Save the result
    write_csv(result, output_path)
    
    cat("Processed file saved as:", output_path, "\n")
  }
  
  return(result)
}
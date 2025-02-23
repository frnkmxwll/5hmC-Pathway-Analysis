###dataset_group_randomizer.R###
# Purpose: Split a DESeq2-ready file into training, validation, and holdout sets.

### INSTALL LIBRARIES IF NECESSARY ----
# install.packages("caret")
# install.packages("dplyr")
# install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("limma")
# BiocManager::install("sva")
# install.packages("doParallel")
# install.packages("foreach")
# install.packages("progress")
# install.packages("glmnet")
# install.packages("ROCR")
# install.packages("future")
# install.packages("future.apply")
# install.packages("plotly")
# install.packages("progressr")
# install.packages("nestedcv")
# install.packages("future")
# BiocManager::install("kaleido")
# BiocManager::install("ChIPseeker")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("annotatr")
# BiocManager::install("wiggleplotr")
# BiocManager::install("metagene2")

## Kaleidoscope for plotly exports
# install.packages('reticulate')
# reticulate::install_miniconda()
#reticulate::conda_install('r-reticulate', 'python-kaleido')
# reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')
# reticulate::use_miniconda('r-reticulate')

# Debug source this file

### Load Libraries ----
library(caret)  # for glmnet model training
library(dplyr)
library(data.table)
#library(limma)
library(DESeq2)  # for normalization and DGE analysis
#library(sva)  # for combat-seq batch normalization
library(future)
library(future.apply)
### Glmnet libraries
library(doParallel)
library(foreach)
#library(progress)
library(glmnet)
library(caret)
library(ROCR)  # Add this line to load the ROCR package
library(future)
library(future.apply)
library(plotly)
library(limma)
library(sva)
### For removing lazy loading on nestedcv
#library(devtools)
#load_all('~/R/x86_64-pc-linux-gnu-library/4.3/nestedcvYM')
#setwd("~/R/x86_64-pc-linux-gnu-library/4.3/nestedcvYM/")
#devtools::install("~/R/x86_64-pc-linux-gnu-library/4.3/nestedcvYM/")
#library(nestedcvYM)
library(nestedcv)
library(ranger)
# For Annotations
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(annotatr)
library(progressr)
library(patchwork)

#library(smotefamily)
#library(DMwR2)

### Configurations ----
source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")

### CONFIGURATION
## General Settings
setwd("~/5hmC-Pathway-Analysis/")
# Location of raw_data_processing output

comp_source = "CANCER_HEALTHY" # "CANCER_HEALTHY" or "PMpositive_PMnegative"
comp_type = "HMR" # "HMR" or "GENE"
## Run settings 
random_conds = TRUE # negative control by randomizing conditions
use_synthetic_data_augmentation <- FALSE
part_rand_seed = seq(1,100, by=1)

plan(list(tweak(multisession, workers = 10))) #sequential, , multicore,multisession

if(comp_type == "HMR"){
  root_folder <- "./Output/Raw Data Processing/HEALTHY_PM_positive_10per_11232024_hmr_combat_healthy_noX/"
  filename_trunk <- "HEALTHY_PM_positive_DESeq2"
} else{
  root_folder <- "~/5hmC-Pathway-Analysis/Output/Raw Data Processing/HEALTHY_PM_positive_genebodies_082624_nocombat_pm_healthy_X/"
  filename_trunk <- "HEALTHY_PM_positive_DESeq2"
}

#  get current date and drop "-"
date <- Sys.Date()
date <- gsub("-", "", date)
file_version <- paste0(comp_source,"_",comp_type,"_12_29_24_",date)

output_folder <- "./Output/repeated_holdout/"
save_final_files <- TRUE
debug_progress <- TRUE
debug_data_summary <- FALSE
random_seed <- 10
set.seed(random_seed)

### Partition Settings ----
interaction_params <- c("condition","sex","age_cat")#"age_cat","sex","primary_present")#,"peritoneal_mets") # for PM: condition, primary_site 
inner_interaction_params <- c("condition","sex","age_cat")#,"peritoneal_mets")#,"peritoneal_mets")# for PM: condition, primary_site 
training_fraction <- 100 # use 100 for nested cv to evaluate inner and outer folds only
validation_fraction <- 0
holdout_fraction <- 0
if (validation_fraction > 0){
  final_train = TRUE
} else {
  final_train = FALSE
}


### DESeq2 Settings (ignored for RF filter) ----
cutoff_type = 2 # 0=padj cutoff, default; 1=lfc & pvalue cutoff; 2 = top stat genes
padj_cutoff = 0.01 # 0.1 default
pvalue_cutoff = 0.001
lfc_cutoff = 0.138 # 0.137504 ~ 10%, 0.263034 ~ 20% change, 0.32 ~ 25%, 0.415 ~ 33%, 0.58 ~ 50% change,
design_formula <- ~ condition + peritoneal_mets #for DESeq2, skipped for Random Forest
perform_deseq2_only <- 0 # 0 = DESEQ2 + Glmnet , 1 = No Repeated Glmnet
feat_occ_thresholds <- seq(0.1, 1,by = 0.01) # for filter thresh for DESeq2, skipped for Random Forest
nfolds_val = 10
standardize_bit = FALSE # should caret center and standardize before training models?
number_of_runs <- 150

### RFE & glmnet Settings ----
alpha_values <- seq(0.1 , 1, by = 0.01) # for filter thresh # 0.3 for sparse, 0.1 for HAC
num_top_genes = 50

median_percentile_cutoff <- 0 # used for pre-filtering
variance_percentile_cutoff <- 0 # used for pre-filtering
n_inner_folds_param = 10
n_outer_folds_param = 10
num_trees = 10000
max_depth = 30
min_node_size = 20
mtry_factor = 0.2
if(num_top_genes < 1000){
  min_freq_stability_plot = 0.5
} else{
  min_freq_stability_plot = 0.99
}

# Normalization settings ----
create_sig_training_counts <- FALSE
combat_norm <- FALSE
length_norm <- FALSE

#closeAllConnections()


counts_name <- paste0(root_folder, filename_trunk, "_rawcounts.csv")
normcounts_name <- paste0(root_folder, filename_trunk, "_normcounts.csv")
meta_name <- paste0(root_folder, filename_trunk, "_conditions.csv")

if(random_conds){
  file_version <- paste0(file_version,"_RAND")
  }

if(comp_source == "CANCER_HEALTHY"){
  healthy_v_diseased = TRUE # Set to True when comparing Healthy vs Diseased
} else {
  healthy_v_diseased = FALSE
}

rapid = FALSE
if(rapid){
  plan(sequential) #sequential, , multicore,multisession
  num_raw_genes = 1000
  num_top_genes = 150
  feat_occ_thresholds <- seq(0, 0.1, by = 0.1) #c(0.7,0.75,0.8,0.85,0.9,0.95,0.98,1)#
  alpha_values <- seq(0.1, 1, by = 0.1)
  part_rand_seed = seq(1,2, by=1)
  save_final_files = FALSE
}

# Define variable to post-pend to output files that incorporates key parameters
interaction_string <- paste(substr(interaction_params, 1, 4), collapse = "-")

file_version <- paste0(file_version,"_",interaction_string)#, number_of_runs, "_nfolds_", nfolds_val, "_num_top_genes_", num_top_genes, "_seeds_", length(part_rand_seed),"_params",interaction_string)
# Create output folder if it does not exist, based on file version and date
output_folder_ver <- paste0(output_folder, file_version, "_", format(Sys.Date(), "%Y%m%d"), "_", format(Sys.time(), "%H%M%S"))# if it doesn't exist, create it
if (!dir.exists(output_folder_ver)) {
  dir.create(output_folder_ver)
}

## Define color tables ----
source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
color_tables <- create_color_tables()

### LOAD AND PREPARE META AND COUNTS DATA ----
## Subselect metadata file if desired. Comment out for PM.
if(healthy_v_diseased){
  subset_conditions <- list(list("primary_site", list("CRC","ADE","healthy"))) # include these types ,"healthy"
}else{
  subset_conditions <- list(list("primary_site", list("CRC","ADE")))
}

#subset_conditions <- NULL
excluded_samples <- NULL #c("KT031", "KT043", "KT046", "KT117", "KT001", "KT003", "KT007", "KT019", "KT024", "KT047", "KT067", "KT068", "KT078", "KT135", "KT136", "KT152", "KT165", "KT168", "KT235", "KT331", "KT332")

data_prepared <- load_and_prepare_data(counts_name, 
                                       normcounts_name,
                                       meta_name, 
                                       debug_data_summary, 
                                       debug_progress, 
                                       rapid, 
                                       num_raw_genes,
                                       subset_conditions = subset_conditions,
                                       excluded_samples = NULL,
                                       random_conds = random_conds,
                                       collapse_healthy = healthy_v_diseased) # set to true when comparing cancer vs healthy

# Add stratum column based on interaction parameters
data_prepared$meta$stratum <- do.call(interaction, data_prepared$meta[interaction_params])
data_prepared$meta$inner_stratum <- do.call(interaction, data_prepared$meta[inner_interaction_params])
conditions <- define_conditions(data_prepared$meta)

### PERFORM NESTEDCV ----
if(TRUE){
  start_time <- Sys.time()
  print(paste("Started at: ",start_time))
  results_list <- NULL
  source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
  # Use future_lapply with progress monitoring
  results_list <- future.apply::future_lapply(part_rand_seed,
                                function(random_seed) {
    cat(paste("Starting on random seed:",random_seed))
    set.seed(random_seed)
  
    # Load and preprocess data
    partitions <- partition_data(data_prepared$meta, interaction_params, training_fraction, validation_fraction, holdout_fraction, random_seed, debug_progress)
    ordered_filtered_counts <- order_and_filter_data(partitions$training_meta, partitions$validation_meta, partitions$holdout_meta, data_prepared$counts_data, median_percentile_cutoff, variance_percentile_cutoff, debug_progress)
    normalized_counts <- normalize_counts(ordered_filtered_counts$training, ordered_filtered_counts$validation, ordered_filtered_counts$holdout, combat_norm, length_norm, debug_progress)
  
    # Prepare metadata and normalized counts for glmnet
    meta_training <- partitions$training_meta
    meta_validation <- partitions$validation_meta
    meta_holdout <- partitions$holdout_meta
    #class1_name <- meta_training[2, 1]
    #class2_name <- meta_training[nrow(meta_training), 1]
  
    norm_counts_data_training <- normalized_counts$training_normcounts
    norm_counts_data_validation <- if (validation_fraction == 0) NULL else normalized_counts$validation_normcounts
    norm_counts_data_holdout <- normalized_counts$holdout_normcounts
  
    conditions_vector_training <- create_conditions_vector(meta_training,"condition")
    conditions_vector_validation <- if (validation_fraction == 0) NULL else create_conditions_vector(meta_validation,"condition")
    conditions_vector_holdout <- if (validation_fraction == 0) NULL else create_conditions_vector(meta_holdout,"condition")
  
    conditions_vector_training_fac <- as.factor(make.names(conditions_vector_training))
    conditions_vector_validation_fac <- if (validation_fraction == 0) NULL else as.factor(make.names(conditions_vector_validation))
    conditions_vector_holdout_fac <- if (validation_fraction == 0) NULL else as.factor(make.names(conditions_vector_holdout))
  
    ### Calculate weights for nestedcv
    class_weights <- weight(conditions_vector_training_fac)
    ## IF WEIGHTING BY STRATUM ##
    names(class_weights) <- meta_training$condition
    # replace PMneg with X0 and PMpos with X1 in names(class_weights)
    if(healthy_v_diseased){
      names(class_weights) <- gsub("HEALTHY", "X0", names(class_weights))
      names(class_weights) <- gsub("DISEASED", "X1", names(class_weights))
    } else{
      names(class_weights) <- gsub("PMneg", "X0", names(class_weights)) 
      names(class_weights) <- gsub("PMpos", "X1", names(class_weights))
    }
    
    ## IF WEIGHTING BY STRATUM ##
    
    ## Calculate weights for RFE
    unique_weights <- unique(class_weights)
    class_labels <- c("X0", "X1")
    
    # Create the named vector for class.weights
    class_weights_vector <- setNames(unique_weights, class_labels)
    
    print("Starting Repeated Fit...")
    stratification_factor <- interaction(meta_training$stratum)
    stratification_factor_inner <- interaction(meta_training$inner_stratum)
    #stratification_factor <- if (is.matrix(conditions_vector_training_fac)) conditions_vector_training_fac[,1] else conditions_vector_training_fac
    out_folds = createFolds(stratification_factor, k = n_outer_folds_param)
    in_folds <- lapply(out_folds, function(i) {
      train_y <- stratification_factor_inner[-i]
      caret::createFolds(train_y, k = n_inner_folds_param)
    })
    
    repeated_fit <- tryCatch({
      nestcv.glmnet(
        conditions_vector_training_fac,
        t(normalized_counts$training_normcounts),
        family = "binomial",
        type.measure = "auc",
        #weights = class_weights,
        standardize = TRUE,
        outer_folds = out_folds,
        inner_folds = in_folds,
        alphaSet = alpha_values,
        cv.cores = 10,
        n_outer_folds = n_outer_folds_param,
        n_inner_folds = n_inner_folds_param,
        outer_method = "cv",
        finalCV = final_train,
        verbose = TRUE,
        min_1se = 1,
        keep = TRUE,
        outer_train_predict = TRUE,
        filterFUN = ranger_filter,
        filter_options = list(
          type = "index",
          num.trees = num_trees,
          max.depth = max_depth,
          min.node.size = min_node_size,
          #importance = "impurity_corrected"
          nfilter = num_top_genes,
          mtry = floor(mtry_factor*nrow(norm_counts_data_holdout)),
          class.weights = class_weights_vector
        )
        # filterFUN = deseq2_nestglm_filter,
        # filter_options = list(
        #   training_data = ordered_filtered_counts$training,
        #   design_formula = design_formula,
        #   cutoff_type = cutoff_type,
        #   padj_cutoff = padj_cutoff,
        #   pvalue_cutoff = pvalue_cutoff,
        #   lfc_cutoff = lfc_cutoff,
        #   num_top_genes = num_top_genes,
        #   norm_counts_data_training = normalized_counts$training_normcounts,
        #   conditions_vector_training = conditions_vector_training,
        #   alpha_values = alpha_values,
        #   number_of_runs = number_of_runs,
        #   nfolds_val = nfolds_val,
        #   standardize_bit = standardize_bit,
        #   class_weights = class_weights,
        #   meta_training = meta_training,
        #   feat_occ_thresholds = feat_occ_thresholds,
        #   random_seed = random_seed,
        #   n_outer_folds = n_outer_folds_param,
        #   n_inner_folds = n_inner_folds_param,
        #   perform_deseq2_only = perform_deseq2_only,
        #   specific_folds_inner = specific_folds
        # )
      )
    }, error = function(e) {
      message("An error occurred during repeated fitting: ", e$message)
      return(NULL)
    })
    
    
    print("Finished repeated fitting...")
    final_fit <- repeated_fit$final_fit
    
    selected_predictors = NULL
    # Try to extract the Dimnames
    try({
      selected_predictors = final_fit$glmnet.fit$beta@Dimnames[[1]]
    }, silent = TRUE)
    # if selected_predictors is NULL, try final_fit$beta@Dimnames[[1]]
    if (is.null(selected_predictors)) {
      try({
        selected_predictors = final_fit$beta@Dimnames[[1]]
      }, silent = TRUE)
    }
    
    # Check if there was an error and print the error message
    if (is.null(selected_predictors)) {
      cat("Error: Failed to extract Dimnames\n")
      pred_val <- NULL
      roc_perf_val <- NULL
      auc_val <- NULL
    }
    
    
    if (!is.null(norm_counts_data_validation) && !is.null(conditions_vector_validation)) {
      cat("Making predictions for: %d", random_seed)
      
      norm_counts_data_validation_top <- norm_counts_data_validation[selected_predictors, , drop = FALSE]
      common_predictors <- intersect(selected_predictors, rownames(norm_counts_data_validation))
      if (length(selected_predictors) == length(common_predictors)) {
        cat("Matched selected_preds and cols in norm_counts")
        probabilities_val <- predict(final_fit, newx = t(norm_counts_data_validation_top), s = repeated_fit$final_param[1], type = "response")
        pred_val <- prediction(probabilities_val, conditions_vector_validation)
        roc_perf_val <- performance(pred_val, "tpr", "fpr")
        auc_val <- as.numeric(performance(pred_val, measure = "auc")@y.values)
      } else {
        cat("Mismatched selected_preds and cols in norm_counts")
        pred_val <- NULL
        roc_perf_val <- NULL
        auc_val <- NULL
      }
    } else {
      pred_val <- NULL
      roc_perf_val <- NULL
      auc_val <- NULL
    }
    
    cat("Saving results for: %d", random_seed)
    
    result <- list(
      part_seed = random_seed,
      num_features = nrow(repeated_fit$final_coef),
      auc_tra_inner = signif(pROC::auc(innercv_roc(repeated_fit)),3),
      auc_tra_outter = signif(pROC::auc(repeated_fit$roc),3), # Outter fold auc
      filter_runs = number_of_runs,
      filter_nfolds = nfolds_val,
      #filter_alpha = best_filter_alpha,
      #filter_threshold = best_filter_threshold,
      n_inner_folds = n_inner_folds_param,
      n_outer_folds = n_outer_folds_param,
      alpha_lambda = repeated_fit$final_param,
      coefficients = repeated_fit$final_coef,
      roc_data_tra_inner = innercv_roc(repeated_fit), # Inner CV
      roc_data_tra_outter = repeated_fit$roc, # outter CV
      repeated_fit = repeated_fit,
      meta_ordered = ordered_filtered_counts$training$meta_ordered
    )
    
    # if final_train = TRUE then also add the final_fit to the result
    if (final_train) {
      result$final_fit <- final_fit
      result$roc_data_val <- roc_perf_val
      result$pred_val <- pred_val
      result$auc_val <- auc_val
    } else {
      result$final_fit <- NULL
      result$roc_data_val <- NULL
      result$pred_val <- NULL
      result$auc_val <- NULL
    }
    
    
    # Build temp results list: append the result to results_list and save it
    #results_list_temp <<- append(results_list_temp, list(result))
    
    return(result)
                                }, future.packages = c("caret", "dplyr", "limma", "DESeq2", "sva", "doParallel", "foreach", "progress", "glmnet", "ROCR", "future.apply", "BiocParallel"), future.seed = TRUE)
  end_time <- Sys.time()
  print(paste("Total Operation Time: ",end_time - start_time))
}

### Prepare results for visualization ----
# Drop failed result from results
# Assuming results_list is your list
if(TRUE){
  results_list <- lapply(results_list, function(item) {
    # Check if the item is not a function
    if (!is.function(item)) {
      return(item)
    } else {
      return(NULL)
    }
  })
  
  # Remove NULL elements from the list
  results_list <- results_list[!sapply(results_list, is.null)]
  results_df = do.call(rbind, lapply(results_list, flatten_result))
  
  # FIND MEDIAN MODEL
  if(FALSE){
    auc_tra_outters <- results_df$auc_tra_outter
    # Calculate the median auc_val
    custom_median <- function(x) {
      n <- length(x)
      if (n %% 2 == 1) {
        return(median(x))
      } else {
        sorted_x <- sort(x)
        middle_index <- n / 2
        return(sorted_x[middle_index + 1])
      }
    }
    median_auc_tra_outter <- custom_median(auc_tra_outters)
    
    # Round values to 3 significant digits
    rounded_auc_vals <- signif(auc_tra_outters, 3)
    rounded_median_auc_tra_outter <- signif(median_auc_tra_outter, 3)
    
    # Find indices of lists with the rounded median auc_val
    median_model_idx <- which(rounded_auc_vals == rounded_median_auc_tra_outter)
    # Extract lists with the median auc_val
    median_lists <- results_list[median_model_idx]
    
    # Find the number of features in each list with the median auc_val
    num_features_vals <- sapply(median_lists, function(x) x$num_features)
    
    # Find the index of the list with the smallest num_features among those with the median auc_val
    index_min_num_features_in_median <- which.min(num_features_vals)
    
    # Find the corresponding index in the original results_list
    median_model_idx <- median_model_idx[index_min_num_features_in_median]
    
    # Extract the list with the median auc_val and the smallest num_features
    median_result <- results_list[[median_model_idx]]
    
    summary(median_result$repeated_fit)
    
    median_training_innner_roc <- median_result$roc_data_tra_inner
    median_training_outter_roc <- median_result$roc_data_tra_outter
  } else{ # SKIP
    median_result <- NULL
    median_training_innner_roc <- NULL
    median_training_outter_roc <- NULL
  }
  
  ### CALCULATE COEFFICIENT FREQUENCY TABLE
  if(FALSE){
    # Pull coefficients of all models into single df
    # Initialize an empty data frame for coefficients
    coefficients_df <- data.frame()
    
    # Iterate over best result for each random seed, to build coefficients table
    for (item in results_list) {
      # Extract the best_coefficients
      best_coefficients <- rownames(item$coefficients)[-1]
      
      # Convert the coefficients to a data frame and add model_id
      coefficients_df_item <- data.frame(
        feature = rownames(item$coefficients)[-1],
        coefficient = item$coefficients$coef[-1],
        model_id = item$part_seed
      )
      
      # Bind the current coefficients to the overall coefficients_df
      coefficients_df <- rbind(coefficients_df, coefficients_df_item)
    }
    
    
    # Drop all rows where the feature is "(Intercept)" and save to new dataframe coefficients_df_features
    coefficients_df_features <- coefficients_df[coefficients_df$feature != "(Intercept)", ]
    
    # Add gene symbols to coefficients_df_features
    coefficients_df_features <- annotate_features(coefficients_df_features)
    
    # Calculate feature appearance percentage
    feature_counts <- coefficients_df_features %>%
      group_by(feature, GeneSymbol) %>%
      summarise(appearance_count = n()) %>%
      mutate(percentage = (appearance_count / length(part_rand_seed)) * 100)
  } else{ # SKIP
    coefficients_df <- NULL
    coefficients_df_features <- NULL
    feature_counts <- NULL
  }
}


### CREATE ROC PLOTS ----
if(TRUE){
  all_training_inner_roc <- sapply(results_list, function(x) x$roc_data_tra_inner) # inner
  all_training_outter_roc <- sapply(results_list, function(x) x$roc_data_tra_outter) # outter
  all_validation_roc <- if (!is.null(results_list[[1]]$roc_data_val)) {
    sapply(results_list, function(x) x$roc_data_val)
  } else {
    NULL
  }
  
  convert_pred_to_roc_data <- function(pred) {
    # Extract required data
    predictor <- pred@predictions[[1]]
    response <- pred@labels[[1]]
    
    # Ensure response is a factor with two levels
    response <- factor(response, levels = c("0", "1"))
    
    # Create ROC object
    roc_obj <- roc(response, predictor, quiet = TRUE)
    
    # Extract data from ROC object
    list(
      percent = FALSE,
      sensitivities = roc_obj$sensitivities,
      specificities = roc_obj$specificities,
      thresholds = roc_obj$thresholds,
      direction = "<",
      cases = as.numeric(roc_obj$cases),
      controls = as.numeric(roc_obj$controls),
      fun.sesp = roc_obj$fun.sesp,
      auc = as.numeric(roc_obj$auc),
      call = roc_obj$call,
      predictor = predictor,
      response = response
    )
  }
  
  # Convert all validation ROC data
  all_validation_roc_list <- lapply(results_list, function(result) {
    if (!is.null(result$pred_val)) {
      convert_pred_to_roc_data(result$pred_val)
    } else {
      NULL
    }
  })
  
  # Remove any NULL entries (if any results didn't have pred_val)
  all_validation_roc_list <- all_validation_roc_list[!sapply(all_validation_roc_list, is.null)]
  
  # Combine the list into a data frame
  all_validation_roc_converted <- as.data.frame(do.call(cbind, all_validation_roc_list))
  
  tryCatch({
    colnames(all_validation_roc_converted) <- paste0("Model", seq_along(all_validation_roc_list))
  }, error = function(e) {
    message("Validation set is empty: ", e$message)
  })
  

  source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
  training_inner_roc_plot <- plot_nested_roc_curves(all_training_inner_roc,  "")
  training_outter_roc_plot <- plot_nested_roc_curves(all_training_outter_roc, "")
  validation_roc_plot <- if (!is.null(all_validation_roc)) {
    plot_nested_roc_curves(all_validation_roc_converted, "")
  } else {
    NULL
  }

  source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")

  # Combine the plots conditionally
  if (!is.null(validation_roc_plot)) {
    combined_roc_plot <- subplot(
      training_inner_roc_plot,
      training_outter_roc_plot,
      validation_roc_plot,
      nrows = 1, titleX = TRUE, titleY = TRUE
    ) %>% layout(
      annotations = list(
        list(
          x = 0.16, y = 1.05, xref = "paper", yref = "paper",
          text = "Inner Folds CV", showarrow = FALSE,
          font = list(size = 16)
        ),
        list(
          x = 0.5, y = 1.05, xref = "paper", yref = "paper",
          text = "Outer Folds CV", showarrow = FALSE,
          font = list(size = 16)
        ),
        list(
          x = 0.84, y = 1.05, xref = "paper", yref = "paper",
          text = "Validation ROC Curves", showarrow = FALSE,
          font = list(size = 16)
        )
      )
    )
  } else {
    combined_roc_plot <- subplot(
      training_inner_roc_plot,
      training_outter_roc_plot,
      nrows = 1, titleX = TRUE, titleY = FALSE
    ) %>% layout(
      annotations = list(
        list(
          x = 0.25, y = 1.05, xref = "paper", yref = "paper",
          text = "Inner Folds CV", showarrow = FALSE,
          font = list(size = 16)
        ),
        list(
          x = 0.75, y = 1.05, xref = "paper", yref = "paper",
          text = "Outer Folds CV", showarrow = FALSE,
          font = list(size = 16)
        )
      )
    )
  }

  # Display the combined plot
  combined_roc_plot
}

### GENERATE ROC PLOTS BY GROUP ----
if(TRUE){
  source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
  
  prediction_dfs <- create_prediction_dataframe_meta(results_list)
  prediction_df_all <- do.call(rbind, prediction_dfs)
  library(pROC)
  library(dplyr)
  
  # Initialize an empty list to hold the ROC data per metadata
  roc_data_per_metadata <- list()
  
  # Define the metadata columns
  metadata_columns <- c("age_cat", "sex", "primary_site", "primary_present", "batch","race_cat")
  
  # Initialize a vector to hold metadata columns to use
  metadata_columns_to_use <- c()
  
  # Determine which metadata columns to use
  for (metadata_column in metadata_columns) {
    unique_values <- unique(prediction_df_all[[metadata_column]])
    
    # For 'primary_present', 'primary_site', 'batch', skip if healthy_v_diseased is TRUE
    if (metadata_column %in% c("primary_present", "primary_site","batch") && healthy_v_diseased) {
      next
    } else {
      metadata_columns_to_use <- c(metadata_columns_to_use, metadata_column)
    }
  }
  
  # Now loop over the filtered metadata columns to generate roc_data_per_metadata
  for (metadata_column in metadata_columns_to_use) {
    unique_values <- unique(prediction_df_all[[metadata_column]])
    roc_data_per_value <- list()
    
    # Generate ROC curves per model/run for each value
    for (value in unique_values) {
      roc_data_per_model <- list()
      for (i in seq_along(prediction_dfs)) {
        prediction_df <- prediction_dfs[[i]]
        df_value <- prediction_df[prediction_df[[metadata_column]] == value, ]
        
        if (nrow(df_value) > 1) {
          labels <- df_value$condition
          labels_numeric <- ifelse(labels %in% c("PM_positive", "DISEASED"), 1, 0)
          probabilities <- df_value$probability
          
          roc_obj <- roc(labels_numeric, probabilities, quiet = TRUE)
          roc_data_per_model[[length(roc_data_per_model) + 1]] <- roc_obj
        }
      }
      roc_data_per_value[[as.character(value)]] <- roc_data_per_model
    }
    roc_data_per_metadata[[metadata_column]] <- roc_data_per_value
  }
  
  # Data frame to store median-run DeLong results
  median_run_delong_results <- data.frame(
    metadata_column = character(),
    value1 = character(),
    value2 = character(),
    p_value = numeric(),
    statistic = numeric(),
    estimate_auc_diff = numeric(),
    conf_level = numeric(),
    conf_interval_lower = numeric(),
    conf_interval_upper = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Data frame to store distribution-based test results (e.g., Wilcoxon test)
  distribution_test_results <- data.frame(
    metadata_column = character(),
    value1 = character(),
    value2 = character(),
    test = character(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
}

## Execute ROC plot by metadata type ----
if(TRUE){
  
  source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
  
  for (metadata_column in names(roc_data_per_metadata)) {
    print("plotting ROC for metadata column:")
    print(metadata_column)
    roc_data_per_value <- roc_data_per_metadata[[metadata_column]]
    title <- paste("Average ROC Curves by", metadata_column)
    
    # Select the appropriate color table
    if (metadata_column == "condition") {
      color_table <- color_tables$condition_table
    } else if (metadata_column == "primary_site") {
      color_table <- color_tables$primary_site_table
    } else if (metadata_column == "primary_present") {
      color_table <- color_tables$primary_present_table
    } else if (metadata_column == "sex") {
      color_table <- color_tables$sex_table
    } else if (metadata_column == "batch") {
      color_table <- color_tables$batch_table
    } else if (metadata_column == "age_cat") {
      color_table <- color_tables$age_cat_table
    } else if (metadata_column == "race") {
      color_table <- color_tables$race_table
    } else if (metadata_column == "race_cat") {
      color_table <- color_tables$race_cat_table
    } else {
      next
    }
    
    # Plot average ROC per value
    p <- plot_average_roc_per_value(roc_data_per_value, title, color_table)
    print(p)
    
    # Save plots
    output_folder_roc <- paste0(output_folder_ver, "/ROC_plots_by_factor")
    dir.create(output_folder_roc, showWarnings = FALSE, recursive = TRUE)
    
    plotly::save_image(p, paste0(output_folder_roc, "/average_roc_", metadata_column, ".png"), width = 1000, height = 1000)
    plotly::save_image(p, paste0(output_folder_roc, "/average_roc_", metadata_column, ".svg"), width = 1000, height = 1000)
    
    value_names <- names(roc_data_per_value)
    
    # Compute AUC distributions for each value from the repeated runs
    auc_distributions <- lapply(roc_data_per_value, function(roc_list) {
      sapply(roc_list, auc)
    })
    
    # Identify median-run ROC for each value (based on median AUC)
    median_run_roc_per_value <- list()
    for (value_name in value_names) {
      aucs <- auc_distributions[[value_name]]
      if (length(aucs) > 0) {
        median_auc <- median(aucs)
        # Find the run closest to the median AUC
        diff_from_median <- abs(aucs - median_auc)
        median_idx <- which.min(diff_from_median)
        median_run_roc_per_value[[value_name]] <- roc_data_per_value[[value_name]][[median_idx]]
      } else {
        median_run_roc_per_value[[value_name]] <- NULL
      }
    }
    
    # Run tests only if >1 category
    if (length(value_names) > 1) {
      cat("Performing statistical tests for metadata column:", metadata_column, "\n")
      
      # Pairwise tests
      for (i in seq_len(length(value_names) - 1)) {
        for (j in seq((i+1), length(value_names))) {
          value_name_i <- value_names[i]
          value_name_j <- value_names[j]
          
          # 1) Median-run DeLong test
          median_roc_obj_i <- median_run_roc_per_value[[value_name_i]]
          median_roc_obj_j <- median_run_roc_per_value[[value_name_j]]
          
          if (!is.null(median_roc_obj_i) && !is.null(median_roc_obj_j)) {
            test_result_median <- roc.test(median_roc_obj_i, median_roc_obj_j, method = "delong")
            
            estimate_auc_diff_median <- as.numeric(test_result_median$estimate)
            conf_int_median <- test_result_median$conf.int
            conf_level_median <- attributes(conf_int_median)$conf.level
            
            median_run_delong_results <- rbind(
              median_run_delong_results,
              data.frame(
                metadata_column = metadata_column,
                value1 = value_name_i,
                value2 = value_name_j,
                p_value = test_result_median$p.value,
                statistic = if (!is.null(test_result_median$statistic)) test_result_median$statistic else NA,
                estimate_auc_diff = estimate_auc_diff_median,
                conf_level = if (!is.null(conf_level_median)) conf_level_median else NA,
                conf_interval_lower = if (!is.null(conf_int_median)) conf_int_median[1] else NA,
                conf_interval_upper = if (!is.null(conf_int_median)) conf_int_median[2] else NA,
                stringsAsFactors = FALSE
              )
            )
          }
          
          # 2) Distribution-based test (Wilcoxon test on AUC distributions)
          auc_i <- auc_distributions[[value_name_i]]
          auc_j <- auc_distributions[[value_name_j]]
          
          if (length(auc_i) > 1 && length(auc_j) > 1) {
            dist_test <- wilcox.test(auc_i, auc_j)
            distribution_test_results <- rbind(
              distribution_test_results,
              data.frame(
                metadata_column = metadata_column,
                value1 = value_name_i,
                value2 = value_name_j,
                test = "Wilcoxon rank-sum",
                p_value = dist_test$p.value,
                stringsAsFactors = FALSE
              )
            )
          }
        }
      }
    }
  }
  # Write out result data frames to CSV
  output_folder_roc <- paste0(output_folder_ver, "/ROC_plots_by_factor")
  dir.create(output_folder_roc, showWarnings = FALSE, recursive = TRUE)
  
  write.csv(median_run_delong_results, file = paste0(output_folder_roc, "/delong_results_median_run.csv"), row.names = FALSE)
}

### CREATE PREDICTION PROBABILITY BOXPLOTS: ----
if (TRUE) {
  source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
  prediction_df <- create_prediction_dataframe(results_list)
  
  if(healthy_v_diseased){
    condition_split <- "Healthy"
  } else{
    condition_split <- "PM"
  }
  
  # Define metadata columns
  metadata_columns <- c("age_cat", "race_cat","sex", "primary_site", "primary_present", "batch")
  
  generate_boxplots <- function(data, metadata_column, condition_split, color_table, color_tables) {
    library(ggplot2)
    library(ggsignif)
    library(dplyr)
    
    # Filter data
    data_filtered <- data %>%
      filter(condition %in% names(color_tables$condition_table)) %>%
      filter(!is.na(condition), !is.na(probability), !is.na(!!sym(metadata_column))) %>%
      mutate(condition = factor(condition, levels = names(color_tables$condition_table)))
    
    # Check levels
    if (length(levels(data_filtered$condition)) < 2) {
      stop("Not enough levels in 'condition' to perform comparison.")
    }
    
    # Define comparisons
    comparisons <- list(c(levels(data_filtered$condition)[1], levels(data_filtered$condition)[2]))
    
    # Calculate medians for labeling
    medians <- data_filtered %>%
      group_by(!!sym(metadata_column), condition) %>%
      summarise(median_prob = median(probability, na.rm = TRUE)) %>%
      ungroup()
    
    # Proceed with plotting
    p <- ggplot(data_filtered, aes(x = condition, y = probability, fill = condition)) +
      geom_boxplot(color = 'black', outlier.shape = NA, alpha = 0.65) +
      # Add median labels
      stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)),
                   position = position_dodge(width = 0.75), 
                   vjust = 0.5, size = 3, color = "white") +
      #geom_jitter(color = 'black',  size = 0.05, width = 0.25, alpha = 0.025) +
      facet_wrap(as.formula(paste("~", metadata_column))) +
      scale_fill_manual(values = color_tables$condition_table) +
      labs(
        title = paste0("Probability by ", condition_split, " and ", metadata_column),
        subtitle = paste0("Total Samples: n=", nrow(data_filtered)),
        y = "Probability",
        x = "Condition"
      ) +
      geom_signif(
        comparisons = comparisons,
        map_signif_level = TRUE,
        test = "wilcox.test",
        textsize = 3,
        tip_length = 0.01
      ) +
      scale_y_continuous(limits = c(0, 1)) +   # Set y-axis limits from 0 to 1
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.title.x = element_blank(),      # Remove x-axis title
        axis.text.x = element_blank(),       # Remove x-axis text
        axis.ticks.x = element_blank(),      # Remove x-axis ticks
        strip.background = element_rect(color = "black", fill = "white"),  # Facet label boxes
        strip.text = element_text(size = 12)  # Facet label text size
      )
    
    # Perform statistical tests for each level of metadata_column
    test_results <- data_filtered %>%
      group_by(!!sym(metadata_column)) %>%
      group_modify(~ {
        df <- .x
        n_levels <- length(unique(df$condition))
        if (n_levels != 2) {
          # Return NA values in expected columns
          tibble(
            p_value = NA_real_,
            statistic = NA_real_,
            method = NA_character_,
            alternative = NA_character_,
            mean_condition1 = NA_real_,
            mean_condition2 = NA_real_,
            note = "Insufficient levels in condition"
          )
        } else {
          test <- wilcox.test(probability ~ condition, data = df)
          mean_values <- df %>%
            group_by(condition) %>%
            summarise(mean_probability = mean(probability, na.rm = TRUE)) %>%
            arrange(condition)
          tibble(
            p_value = test$p.value,
            statistic = test$statistic,
            method = test$method,
            alternative = test$alternative,
            mean_condition1 = mean_values$mean_probability[1],
            mean_condition2 = mean_values$mean_probability[2],
            note = NA_character_
          )
        }
      }) %>%
      ungroup() %>%
      mutate(facet = as.character(!!sym(metadata_column))) %>%
      dplyr::select(
        facet,
        p_value,
        statistic,
        method,
        alternative,
        mean_condition1,
        mean_condition2,
        note
      )
    
    return(list(plot = p, test_results = test_results))
  }
  
  # Function to generate boxplot for 'condition' column
  generate_condition_boxplot <- function(data, condition_split, condition_table) {
    library(ggplot2)
    library(ggsignif)
    library(dplyr)
    
    # Filter data based on condition_split
    if (condition_split == "PM") {
      data_filtered <- data %>%
        filter(condition %in% c("PM_positive", "PM_negative")) %>%
        mutate(condition = factor(condition, levels = c("PM_negative","PM_positive")))
      title <- "Probability by PM Status"
    } else if (condition_split == "Healthy") {
      data_filtered <- data %>%
        filter(condition %in% c("HEALTHY", "DISEASED")) %>%
        mutate(condition = factor(condition, levels = c("HEALTHY", "DISEASED")))
      title <- "Probability by Health Status"
    } else {
      stop("Invalid condition_split. Choose either 'PM' or 'Healthy'.")
    }
    
    # Define comparisons for significance testing
    comparisons <- list(c(levels(data_filtered$condition)[1], levels(data_filtered$condition)[2]))
    
    # Calculate medians for labeling
    medians <- data_filtered %>%
      group_by(condition) %>%
      summarise(median_prob = median(probability, na.rm = TRUE)) %>%
      ungroup()
    
    # Calculate total sample size for subtitle
    total_samples <- nrow(data_filtered)
    
    # Generate boxplot without faceting
    p <- ggplot(data_filtered, aes(x = condition, y = probability, fill = condition)) +
      geom_boxplot(color = 'black', outlier.shape = NA, alpha = 0.65) +
      # Add median labels
      stat_summary(fun = median, geom = "text", aes(label = round(..y.., 2)),
                   position = position_dodge(width = 0.75), 
                   vjust = 0.5, size = 3, color = "white") +
      #geom_jitter(color = 'black', size = 0.05, width = 0.25, alpha = 0.025) +
      scale_fill_manual(values = condition_table) +
      labs(
        title = title,
        subtitle = paste0("Total Samples: n=", total_samples),
        y = "Probability",
        x = "Condition"
      ) +
      geom_signif(comparisons = comparisons, 
                  map_signif_level = TRUE, 
                  test = "wilcox.test",
                  textsize = 3,
                  tip_length = 0.01) +
      scale_y_continuous(limits = c(0, 1)) +   # Set y-axis limits from 0 to 1
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),     # Remove x-axis text
        axis.title.x = element_blank(),    # Remove x-axis title
        axis.ticks.x = element_blank(),    # Remove x-axis ticks
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        legend.position = "none"
      )
    
    # Perform the statistical test
    test <- wilcox.test(probability ~ condition, data = data_filtered)
    mean_values <- data_filtered %>%
      group_by(condition) %>%
      summarise(mean_probability = mean(probability, na.rm = TRUE))
    
    test_results <- data.frame(
      metadata_column = "condition",
      facet = NA,
      p_value = test$p.value,
      statistic = test$statistic,
      method = test$method,
      alternative = test$alternative,
      mean_condition1 = mean_values$mean_probability[mean_values$condition == levels(data_filtered$condition)[1]],
      mean_condition2 = mean_values$mean_probability[mean_values$condition == levels(data_filtered$condition)[2]]
    )
    
    return(list(plot = p, test_results = test_results))
  }
  
  # Initialize lists to store plots and test results
  boxplot_list <- list()
  test_results_list <- list()
  
  # Generate boxplot and test results for 'condition' first
  result_condition <- generate_condition_boxplot(
    data = prediction_df,
    condition_split = condition_split,
    condition_table = color_tables$condition_table
  )
  p_condition <- result_condition$plot
  test_results_condition <- result_condition$test_results
  boxplot_list[["condition"]] <- p_condition
  test_results_list[["condition"]] <- test_results_condition
  
  # Loop through each metadata column and generate boxplots and test results
  for (metadata_column in metadata_columns) {
    result <- generate_boxplots(
      data = prediction_df, 
      metadata_column = metadata_column, 
      condition_split = condition_split, 
      color_table = color_tables[[paste0(metadata_column, "_table")]],
      color_tables = color_tables
    )
    p <- result$plot
    test_results <- result$test_results
    boxplot_list[[metadata_column]] <- p
    test_results_list[[metadata_column]] <- test_results
  }
  
  # Ensure the condition plot is displayed first
  boxplot_list <- boxplot_list[c("condition", metadata_columns)]
  
  # Combine all test results into a single data frame
  test_results_df <- bind_rows(test_results_list)
  # if output_folder_ver /prob_boxes/ does not exist, create it
  # Check if the folder exists; if not, create it
  output_dir <- file.path(output_folder_ver, "prob_boxes")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Write test results to CSV
  write.csv(test_results_df, file = paste0(output_folder_ver, "/prob_boxes/statistical_test_results_", condition_split, ".csv"), row.names = FALSE)
  
  # Arrange all boxplots into subplots using patchwork
  combined_boxplots <- wrap_plots(boxplot_list, ncol = 2) + 
    plot_layout(guides = 'collect') & 
    theme(legend.position = 'bottom')
  
  # Display the combined boxplots
  print(combined_boxplots)
  
  # Save the combined boxplots as PNG and SVG
  ggsave(
    filename = paste0(output_folder_ver, "/prob_boxes/combined_boxplots_", condition_split, ".png"),
    plot = combined_boxplots, 
    width = 8, height = 12, dpi = 300
  )
  ggsave(
    filename = paste0(output_folder_ver, "/prob_boxes/combined_boxplots_", condition_split, ".svg"),
    plot = combined_boxplots, 
    width = 8, height = 12
  )
} # Version 1


if(TRUE){
  source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
  prediction_df <- create_prediction_dataframe(results_list)
  
  # If conditions include DISEASED and HEALTHY
  # set probability to 1-probabilit
  if("HEALTHY" %in% unique(prediction_df$condition) & "DISEASED" %in% unique(prediction_df$condition)){
    prediction_df$probability <- 1 - prediction_df$probability
  }

  source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
  combined_prob_boxplot_pm <- generate_boxplots_ggplot_pm(prediction_df, color_tables)
  
  source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
  combined_prob_boxplot_healthy <- generate_boxplots_ggplot_healthy(prediction_df, color_tables)
  
  # -------------------------------
  # Calculate Statistics for PM Boxplots
  # -------------------------------
  
  # Plot 1: Peritoneal Metastasis Status (Grouped by 'condition')
  pm_stats_p1 <- prediction_df %>%
    group_by(condition) %>%
    summarise(
      mean_probability = mean(probability, na.rm = TRUE),
      median_probability = median(probability, na.rm = TRUE),
      sd_probability = sd(probability, na.rm = TRUE)
    ) %>%
    mutate(plot = "Peritoneal Metastasis Status")
  
  # Plot 2: Peritoneal Metastasis Status and Presence of Primary (Grouped by 'condition' and 'primary_present')
  pm_stats_p2 <- prediction_df %>%
    group_by(condition, primary_present) %>%
    summarise(
      mean_probability = mean(probability, na.rm = TRUE),
      median_probability = median(probability, na.rm = TRUE),
      sd_probability = sd(probability, na.rm = TRUE)
    ) %>%
    mutate(plot = "Peritoneal Metastasis Status and Presence of Primary")
  
  # Plot 3: Peritoneal Metastasis Status and Primary Site (Grouped by 'condition' and 'primary_site')
  pm_stats_p3 <- prediction_df %>%
    group_by(condition, primary_site) %>%
    summarise(
      mean_probability = mean(probability, na.rm = TRUE),
      median_probability = median(probability, na.rm = TRUE),
      sd_probability = sd(probability, na.rm = TRUE)
    ) %>%
    mutate(plot = "Peritoneal Metastasis Status and Primary Site")
  
  # Combine all PM statistics into one dataframe
  pm_stats <- bind_rows(pm_stats_p1, pm_stats_p2, pm_stats_p3)
  
  # -------------------------------
  # Calculate Statistics for Healthy Boxplots
  # -------------------------------
  
  # Plot 1: Condition (Grouped by 'condition')
  healthy_stats_p1 <- prediction_df %>%
    group_by(condition) %>%
    summarise(
      mean_probability = mean(probability, na.rm = TRUE),
      median_probability = median(probability, na.rm = TRUE),
      sd_probability = sd(probability, na.rm = TRUE)
    ) %>%
    mutate(plot = "Condition")
  
  # Plot 2: Peritoneal Metastasis Status (Grouped by 'condition' and 'peritoneal_mets')
  healthy_stats_p2 <- prediction_df %>%
    group_by(condition, peritoneal_mets) %>%
    summarise(
      mean_probability = mean(probability, na.rm = TRUE),
      median_probability = median(probability, na.rm = TRUE),
      sd_probability = sd(probability, na.rm = TRUE)
    ) %>%
    mutate(plot = "Peritoneal Metastasis Status")
  
  # Plot 3: Primary Tumor Presence (Grouped by 'condition' and 'primary_present')
  healthy_stats_p3 <- prediction_df %>%
    group_by(condition, primary_present) %>%
    summarise(
      mean_probability = mean(probability, na.rm = TRUE),
      median_probability = median(probability, na.rm = TRUE),
      sd_probability = sd(probability, na.rm = TRUE)
    ) %>%
    mutate(plot = "Primary Tumor Presence")
  
  # Plot 4: Primary Site (Grouped by 'condition' and 'primary_site')
  healthy_stats_p4 <- prediction_df %>%
    group_by(condition, primary_site) %>%
    summarise(
      mean_probability = mean(probability, na.rm = TRUE),
      median_probability = median(probability, na.rm = TRUE),
      sd_probability = sd(probability, na.rm = TRUE)
    ) %>%
    mutate(plot = "Primary Site")
  
  # Plot 5: Age Category (Grouped by 'condition' and 'age_cat')
  healthy_stats_p5 <- prediction_df %>%
    group_by(condition, age_cat) %>%
    summarise(
      mean_probability = mean(probability, na.rm = TRUE),
      median_probability = median(probability, na.rm = TRUE),
      sd_probability = sd(probability, na.rm = TRUE)
    ) %>%
    mutate(plot = "Age Category")
  
  # Plot 6: Sex (Grouped by 'condition' and 'sex')
  healthy_stats_p6 <- prediction_df %>%
    group_by(condition, sex) %>%
    summarise(
      mean_probability = mean(probability, na.rm = TRUE),
      median_probability = median(probability, na.rm = TRUE),
      sd_probability = sd(probability, na.rm = TRUE)
    ) %>%
    mutate(plot = "Sex")
  
  # Plot 7: Batch (Grouped by 'condition' and 'batch')
  healthy_stats_p7 <- prediction_df %>%
    group_by(condition, batch) %>%
    summarise(
      mean_probability = mean(probability, na.rm = TRUE),
      median_probability = median(probability, na.rm = TRUE),
      sd_probability = sd(probability, na.rm = TRUE)
    ) %>%
    mutate(plot = "Batch")
  
  # Combine all Healthy statistics into one dataframe
  healthy_stats <- bind_rows(
    healthy_stats_p1, 
    healthy_stats_p2, 
    healthy_stats_p3, 
    healthy_stats_p4, 
    healthy_stats_p5, 
    healthy_stats_p6, 
    healthy_stats_p7
  )
  ### Added Code Ends Here ###
} # version 2, both are good but different

### CREATE CONFUSION MATRIX OF OUTER FOLDS: ----
if(TRUE){source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
library(caret)
library(ggplot2)
library(boot)

# Function to calculate performance metrics
process_result <- function(result) {
  if (is.null(result) || is.null(result$roc_data_tra_outter)) {
    return(NULL)
  }
  
  roc_data <- result$roc_data_tra_outter
  
  # Extract predictions and actual values
  predictions <- roc_data$predictor
  actual <- roc_data$response
  
  # Convert probabilities to binary predictions
  binary_preds <- factor(ifelse(predictions > 0, "X1", "X0"), levels = levels(actual))
  # Create confusion matrix
  cm <- confusionMatrix(binary_preds, actual,positive = "X1")
  
  # Extract metrics with custom names
  metrics <- c(
    TP = cm$table[2, 2],
    FP = cm$table[2, 1],
    TN = cm$table[1, 1],
    FN = cm$table[1, 2],
    Sensitivity = cm$byClass["Sensitivity"],
    Specificity = cm$byClass["Specificity"],
    Accuracy = cm$overall["Accuracy"],
    Precision = cm$byClass["Precision"],
    F1 = cm$byClass["F1"],
    AUC = as.numeric(roc_data$auc),
    NPV = cm$byClass["Neg Pred Value"]
  )
  
  # Return the metrics as a named vector without repeating names
  names(metrics) <- c("TP", "FP", "TN", "FN", "Sensitivity", "Specificity", "Accuracy", "Precision", "F1", "AUC","NPV")
  
  return(metrics)
}

# Process all results
all_metrics <- lapply(results_list, process_result)
all_metrics <- do.call(rbind, all_metrics[!sapply(all_metrics, is.null)])

# Function to calculate confidence intervals
calculate_ci <- function(data, indices) {
  d <- data[indices,]
  return(colMeans(d))
}

# Calculate bootstrap confidence intervals
set.seed(123)  # for reproducibility
boot_results <- boot(all_metrics, calculate_ci, R = 1000)

# Calculate mean metrics and confidence intervals
mean_metrics <- colMeans(all_metrics)
ci_metrics <- t(sapply(1:ncol(all_metrics), function(i) {
  ci <- boot.ci(boot_results, type = "perc", index = i)
  c(ci$percent[4], ci$percent[5])
}))

# Combine results
cm_results_df <- data.frame(
  Metric = names(mean_metrics),
  Mean = mean_metrics,
  CI_Lower = ci_metrics[,1],
  CI_Upper = ci_metrics[,2]
)

# Function to create average confusion matrix
create_average_confusion_matrix <- function(results_list) {
  all_cms <- lapply(results_list, function(result) {
    if (is.null(result) || is.null(result$roc_data_tra_outter)) return(NULL)
    predictions <- result$roc_data_tra_outter$predictor
    actual <- result$roc_data_tra_outter$response
    binary_preds <- factor(ifelse(predictions > 0, "X1", "X0"), levels = levels(actual))
    table(binary_preds, actual)
  })
  all_cms <- do.call(abind::abind, c(all_cms[!sapply(all_cms, is.null)], along = 3))
  avg_cm <- apply(all_cms, c(1, 2), mean)
  return(avg_cm)
}
avg_cm <- create_average_confusion_matrix(results_list)

# Print results
print(cm_results_df)

source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
library(caret)
library(ggplot2)
library(boot)
library(gridExtra)
create_enhanced_confusion_matrix_plot <- function(cm, filename, conditions, cm_results_df) {
  library(ggplot2)
  library(gridExtra)
  
  # Define the mapping from X0 and X1 to actual condition names
  condition_mapping <- c(X0 = as.character(conditions$first_condition),
                         X1 = as.character(conditions$second_condition))
  
  # Convert the matrix to a data frame
  cm_df <- as.data.frame(as.table(cm))
  names(cm_df) <- c("Predicted", "Actual", "Freq")
  
  # Replace X0 and X1 with actual condition names
  cm_df$Predicted <- factor(condition_mapping[cm_df$Predicted], 
                            levels = condition_mapping)
  cm_df$Actual <- factor(condition_mapping[cm_df$Actual], 
                         levels = condition_mapping)
  
  # Calculate percentages
  total_samples <- sum(cm_df$Freq)
  cm_df$Percentage <- cm_df$Freq / total_samples * 100
  
  # Extract metrics with CI
  get_metric_with_ci <- function(metric) {
    mean_val <- cm_results_df$Mean[cm_results_df$Metric == metric] * 100
    ci_lower <- cm_results_df$CI_Lower[cm_results_df$Metric == metric] * 100
    ci_upper <- cm_results_df$CI_Upper[cm_results_df$Metric == metric] * 100
    sprintf("%.1f%%\n(95%%CI: %.1f%% - %.1f%%)", mean_val, ci_lower, ci_upper)
  }
  
  sensitivity <- get_metric_with_ci("Sensitivity")
  specificity <- get_metric_with_ci("Specificity")
  accuracy <- get_metric_with_ci("Accuracy")
  precision <- get_metric_with_ci("Precision")
  npv <- get_metric_with_ci("NPV")  # Approximation for binary classification
  
  # Create a 3x3 grid for the enhanced confusion matrix
  grid_data <- expand.grid(x = 1:3, y = 1:3)
  grid_data$value <- NA
  grid_data$label <- NA
  grid_data$border <- TRUE
  
  # Remove borders for Pred. Metrics row and Class Metrics column
  grid_data$border[grid_data$x == 3 | grid_data$y == 3] <- FALSE
  
  # Fill in the confusion matrix values
  grid_data$value[grid_data$x == 1 & grid_data$y == 1] <- cm_df$Freq[cm_df$Predicted == condition_mapping["X1"] & cm_df$Actual == condition_mapping["X1"]]
  grid_data$value[grid_data$x == 2 & grid_data$y == 1] <- cm_df$Freq[cm_df$Predicted == condition_mapping["X0"] & cm_df$Actual == condition_mapping["X1"]]
  grid_data$value[grid_data$x == 1 & grid_data$y == 2] <- cm_df$Freq[cm_df$Predicted == condition_mapping["X1"] & cm_df$Actual == condition_mapping["X0"]]
  grid_data$value[grid_data$x == 2 & grid_data$y == 2] <- cm_df$Freq[cm_df$Predicted == condition_mapping["X0"] & cm_df$Actual == condition_mapping["X0"]]
  
  # Calculate percentages for confusion matrix cells
  grid_data$percentage <- grid_data$value / total_samples * 100
  
  # Calculate min and max percentages for color scale
  min_percent <- floor(min(grid_data$percentage[!is.na(grid_data$percentage)]) / 10) * 10
  max_percent <- ceiling(max(grid_data$percentage[!is.na(grid_data$percentage)]) / 10) * 10
  
  # Add labels with percentages
  grid_data$label[grid_data$x == 1 & grid_data$y == 1] <- sprintf("TP\n%.0f\n(%.1f%%)", grid_data$value[grid_data$x == 1 & grid_data$y == 1], grid_data$percentage[grid_data$x == 1 & grid_data$y == 1])
  grid_data$label[grid_data$x == 2 & grid_data$y == 1] <- sprintf("FN\n%.0f\n(%.1f%%)", grid_data$value[grid_data$x == 2 & grid_data$y == 1], grid_data$percentage[grid_data$x == 2 & grid_data$y == 1])
  grid_data$label[grid_data$x == 1 & grid_data$y == 2] <- sprintf("FP\n%.0f\n(%.1f%%)", grid_data$value[grid_data$x == 1 & grid_data$y == 2], grid_data$percentage[grid_data$x == 1 & grid_data$y == 2])
  grid_data$label[grid_data$x == 2 & grid_data$y == 2] <- sprintf("TN\n%.0f\n(%.1f%%)", grid_data$value[grid_data$x == 2 & grid_data$y == 2], grid_data$percentage[grid_data$x == 2 & grid_data$y == 2])
  
  # Add performance metrics with CI
  grid_data$label[grid_data$x == 3 & grid_data$y == 1] <- sprintf("Sensitivity\n%s", sensitivity)
  grid_data$label[grid_data$x == 3 & grid_data$y == 2] <- sprintf("Specificity\n%s", specificity)
  grid_data$label[grid_data$x == 1 & grid_data$y == 3] <- sprintf("Precision\n%s", precision)
  grid_data$label[grid_data$x == 2 & grid_data$y == 3] <- sprintf("NPV\n%s", npv)
  grid_data$label[grid_data$x == 3 & grid_data$y == 3] <- sprintf("Accuracy\n%s", accuracy)
  
  # Create the plot
  plot <- ggplot(grid_data, aes(x = x, y = y, fill = percentage)) +
    geom_tile(aes(width = ifelse(border, 0.95, 1), height = ifelse(border, 0.95, 1)),
              color = ifelse(grid_data$border, "black", NA)) +
    geom_text(aes(label = label), size = 3, fontface = "bold") +
    scale_fill_gradient(low = "white", high = "steelblue", na.value = "white", 
                        name = "Percentage", 
                        labels = scales::percent_format(scale = 1),
                        limits = c(min_percent, max_percent)) +
    theme_void() +
    theme(aspect.ratio = 1) +
    scale_x_continuous(breaks = 1:3, labels = c(condition_mapping["X1"], condition_mapping["X0"], "Class Metrics"), position = "top") +
    scale_y_continuous(breaks = 1:3, labels = c(condition_mapping["X1"], condition_mapping["X0"], "Pred. Metrics"), trans = "reverse") +
    labs(x = "Actual Class", y = "Predicted Class", title = "Enhanced Confusion Matrix") +
    theme(axis.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0),
          axis.title.x = element_text(margin = margin(b = 10)),
          axis.title.y = element_text(margin = margin(r = 10)),
          plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 10)),
          legend.position = "right")
  
  # Return the plot object
  return(plot)
}

# Create and save the enhanced confusion matrix plot
enhanced_cm_plot <- create_enhanced_confusion_matrix_plot(avg_cm, "enhanced_confusion_matrix.png", conditions, cm_results_df)

# Display the plot
print(enhanced_cm_plot)
}

### CREATE COEFFICIENT STABILITY PLOT ----
if(TRUE){
  source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
  ### Incorporate feature stability
  tryCatch({
    # Function to extract feature coefficients from outer folds.a
    # cv_coef(results_list[[1]]$repeated_fit)
    all_model_stability <- create_stability_plot(results_list, min_freq_stability_plot*length(part_rand_seed)*(n_outer_folds_param),length(part_rand_seed)*(n_outer_folds_param))
    # output first list element
    all_model_stability[[1]]
    stability_df_grouped <- all_model_stability[[2]]
    # sort by frequency_percentage
    stability_df_grouped <- stability_df_grouped[order(-stability_df_grouped$frequency_percentage),]
    # reset index
    rownames(stability_df_grouped) <- NULL
    
    # No thresh plot
    all_model_stability_no_thresh <- create_stability_plot(results_list, 0*length(part_rand_seed)*(n_outer_folds_param),length(part_rand_seed)*(n_outer_folds_param))
    
    ### Annotate enhancers:
    # Load enhancer lookup file
    enhancer_lookup_df <- fread("~/reference_genomes/GeneHancer_AnnotSV_gene_association_scores_v5.20.txt")
    enhancer_lookup_df <- enhancer_lookup_df[is_elite == 1]
    
    # Create a lookup table for GH symbols
    gh_lookup <- enhancer_lookup_df[startsWith(GHid, "GH") & nchar(GHid) >= 6, .(GHid, symbol)]
    setkey(gh_lookup, GHid)
    
    data <- all_model_stability_no_thresh[[2]][, c("feature", setdiff(names(all_model_stability_no_thresh[[2]]), "feature"))]
    # Function to process symbols efficiently
    process_symbols <- function(symbols, gh_lookup) {
      symbols_list <- strsplit(symbols, ", ")[[1]]
      gh_symbols <- symbols_list[startsWith(symbols_list, "GH") & nchar(symbols_list) >= 6]
      non_gh_symbols <- symbols_list[!(symbols_list %in% gh_symbols)]
      
      if (length(gh_symbols) > 0) {
        matched_symbols <- gh_lookup[gh_symbols, on = "GHid"]$symbol
        return(c(matched_symbols, non_gh_symbols))
      } else {
        return(symbols_list)
      }
    }
    
    setDT(data)
    
    annotated_df_enh <- data[, .(
      symbols_gh_full = list(process_symbols(all_symbols, gh_lookup))
    ), by = names(data)]
    
    annotated_df_enh <- annotated_df_enh[, .(
      feature = rep(feature, sapply(symbols_gh_full, length)),
      symbols_gh_full = unlist(symbols_gh_full)
    ), by = setdiff(names(annotated_df_enh), c("feature", "symbols_gh_full"))]
    
    
  }, error = function(e) {
    message("An error occurred during stability plot creation: ", e$message)
    all_model_stability <- NULL
    stability_df_grouped <- NULL
  })
}

### CREATE FINAL SIG GENE TABLES ----
if(TRUE){
  source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
  sig_gene_table <- process_csv(NULL,annotated_df_enh)
  #"~/5hmC-Pathway-Analysis/Output/repeated_holdout/NESTEDCV_v3_num_runs_paramscond-prim_20241014_114901/annotated_df_enh.csv"
}

### GENERATE HEATMAPS ----
if(TRUE){
  source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
  tryCatch({
    combined_heatmap <- generate_heatmaps(training_data = data_prepared$normcounts_data,
                                        meta_training = data_prepared$meta,
                                        annotated_features = stability_df_grouped,
                                        color_tables= color_tables)

    combined_heatmap
  }, error = function(e) {
    message("An error occurred during heatmap generation: ", e$message)
    combined_heatmap <- NULL
  })
}

### Generate enrichment boxplots by gene ----
source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")

# SEt to true if you want to include healthy in plot
if(TRUE){
  if(healthy_v_diseased != TRUE){
  subset_conditions_gene <- list(list("primary_site", list("CRC","ADE","healthy"))) # include these types ,"healthy"
  
  #meta_name <- "~/5hmC-Pathway-Analysis/Output/Raw Data Processing/HEALTHY_PM_positive_10per_080424_hmr_nocombat_pm_healthy_X/HEALTHY_PM_positive_DESeq2_conditions.csv"
  #normcounts_name <- "~/5hmC-Pathway-Analysis/Output/Raw Data Processing/HEALTHY_PM_positive_10per_080424_hmr_nocombat_pm_healthy_X/HEALTHY_PM_positive_DESeq2_normcounts.csv"
  #counts_name <- "~/5hmC-Pathway-Analysis/Output/Raw Data Processing/HEALTHY_PM_positive_10per_080424_hmr_nocombat_pm_healthy_X/HEALTHY_PM_positive_DESeq2_rawcounts.csv"
  
  data_prepared_gene <- load_and_prepare_data(counts_name, 
                                         normcounts_name,
                                         meta_name, 
                                         debug_data_summary, 
                                         debug_progress, 
                                         rapid, 
                                         num_raw_genes,
                                         subset_conditions = subset_conditions_gene,
                                         excluded_samples = NULL,
                                         random_conds = random_conds,
                                         collapse_healthy = FALSE) # set to true when comparing cancer vs healthy
  
  
  # Add stratum column based on interaction parameters
  data_prepared_gene$meta$stratum <- do.call(interaction, data_prepared_gene$meta[interaction_params])
  data_prepared_gene$meta$inner_stratum <- do.call(interaction, data_prepared_gene$meta[inner_interaction_params])
  conditions <- define_conditions(data_prepared_gene$meta)
  }
}

if(TRUE){
source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
tryCatch({
  boxplots <- generate_boxplots_by_gene(training_data = data_prepared_gene$normcounts_data,
                                meta_training = data_prepared_gene$meta,
                                annotated_features = stability_df_grouped,
                                color_tables=color_tables)
  
  # View the statistical results
  print(boxplots$stats_results)
  
  # Print the negative mean plot
  print(boxplots$negative)
  
  # Print the positive mean plot
  print(boxplots$positive)
  
  # Print the faceted plot for genes with negative mean
  print(boxplots$faceted_negative)
  
  # Print the faceted plot for genes with positive mean
  print(boxplots$faceted_positive)
  
}, error = function(e) {
  message("An error occurred during boxplot generation: ", e$message)
  print(e$call)
  print(str(e))
})
}

### PLOT ENRICHMENT ON INDIVIDUAL TRANSCRIPTS ----
if(FALSE){
  library("wiggleplotr")
  library("dplyr")
  library("GenomicRanges")
  library("GenomicFeatures")
  require("org.Hs.eg.db")
  require("TxDb.Hsapiens.UCSC.hg19.knownGene")
  # install metagene2 with bioconductor
  library("metagene2")
  
  plot_transcript_coverage <- function(stability_df_grouped, tx_choice = 1, rescale_introns = FALSE) {
    # transcript_ids <- stability_df_grouped$annot.tx_id where locationType != Enhancer or Intergenic
    transcript_ids <- stability_df_grouped$annot.tx_id[stability_df_grouped$locationType != "Enhancer" & stability_df_grouped$locationType != "Intergenic"]
    # Extract start and end coordinates where locationType is neither Enhancer nor Intergenic
    starts <- stability_df_grouped$start[stability_df_grouped$locationType != "Enhancer" & stability_df_grouped$locationType != "Intergenic"]
    ends <- stability_df_grouped$end[stability_df_grouped$locationType != "Enhancer" & stability_df_grouped$locationType != "Intergenic"]
    # Construct regions
    regions <- data.frame(start = starts, end = ends)
    # subtract 1000 from start and add 1000 to end
    regions$start <- regions$start - 5000
    regions$end <- regions$end + 5000
    # Load TxDb object for hg19
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  
    # Filter transcript_ids based on locationType
    transcript_ids <- stability_df_grouped$annot.tx_id[stability_df_grouped$locationType != "Enhancer" & stability_df_grouped$locationType != "Intergenic"]
  
    # Extract start and end coordinates
    starts <- stability_df_grouped$start[stability_df_grouped$locationType != "Enhancer" & stability_df_grouped$locationType != "Intergenic"]
    ends <- stability_df_grouped$end[stability_df_grouped$locationType != "Enhancer" & stability_df_grouped$locationType != "Intergenic"]
  
    # Construct regions and adjust coordinates
    regions <- data.frame(start = starts, end = ends)
    # find region center and add/subtract 5000 for new start and end
    regions$center <- floor((regions$start + regions$end) / 2)
    regions$start <- regions$center - 5000
    regions$end <- regions$center + 5000
  
    # Load TxDb object for hg19
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  
    # Extract exons, cds, and metadata for the chosen transcript
    exons <- exonsBy(txdb, by = "tx", use.names = TRUE)[transcript_ids[[tx_choice]]]
    cdss <- cdsBy(txdb, by = "tx", use.names = TRUE)[transcript_ids[[tx_choice]]]
  
    # Prepare metadata as a data frame
    metadata_df <- data.frame(
      transcript_id = stability_df_grouped$annot.tx_id[stability_df_grouped$locationType != "Enhancer" & stability_df_grouped$locationType != "Intergenic"][[tx_choice]],
      gene_id = stability_df_grouped$annot.gene_id[stability_df_grouped$locationType != "Enhancer" & stability_df_grouped$locationType != "Intergenic"][[tx_choice]],
      gene_name = stability_df_grouped$annot.symbol[stability_df_grouped$locationType != "Enhancer" & stability_df_grouped$locationType != "Intergenic"][[tx_choice]],
      strand = "*"
    )
  
    # Plot transcripts
    plotTranscripts(exons, cdss, metadata_df, rescale_introns = rescale_introns)
  
    # Prepare meta_tracks
    meta_tracks <- data_prepared$meta
    meta_tracks$scaling_factor <- 1
    meta_tracks$sample_id <- rownames(meta_tracks)
    meta_tracks <- meta_tracks[, c(ncol(meta_tracks), 1:(ncol(meta_tracks)-1))]
    meta_tracks$sample_id <- as.factor(meta_tracks$sample_id)
    rownames(meta_tracks) <- NULL
  
    meta_tracks <- meta_tracks %>%
      dplyr::mutate(bigWig = paste0("/home/turagapm/peritoneal_processing/trimmed_data_bw_9_18_2023_bigwig_CPM/", sample_id, "_bowtie2.bw"))
  
    # Prepare track_data
    track_data <- dplyr::mutate(meta_tracks, track_id = "PM_status", colour_group = condition)
  
    # Plot coverage
    plotCoverage(exons, cdss, metadata_df, track_data,
                 heights = c(2, 1), fill_palette = getGenotypePalette(), coverage_type = "line", rescale_introns = rescale_introns, region_coords = unlist(regions[tx_choice, ]))
  }
  
  plot_transcript_coverage(stability_df_grouped, tx_choice = 1)
}

### PLOTTNG METAGENE PLOTS FOR SPECIFIC REGIONS ----
if(FALSE){
  regions <- NULL
  regions <- data.frame(
    seqnames = stability_df_grouped$chr,
    start = stability_df_grouped$start,
    end = stability_df_grouped$end,
    coefs = stability_df_grouped$mean,
    strand = "*",
    name = stability_df_grouped$annot.symbol,
    frequency_percentage = stability_df_grouped$frequency_percentage,
    score = 0
  )
  
  # Ensure the column names are correctly set
  colnames(regions) <- c("seqnames","start", "end", "coefs", "strand", "name", "frequency_percentage","score")
  
  # Split the regions data frame into two based on the 'coefs' column where coefs >=0 and where frequency_percentage > 25
  regions_pos <- regions[regions$coefs >= 0 & regions$frequency_percentage >= 0, ]
  regions_neg <- regions[regions$coefs < 0 & regions$frequency_percentage >= 0, ]
  
  
  # Function to convert a data frame to a GRanges object
  df_to_granges <- function(df) {
    GRanges(
      seqnames = df$seqnames,  # Replace with appropriate chromosome info if available
      ranges = IRanges(start = df$start, end = df$end),
      strand = Rle(df$strand),
      name = df$name,
      score = df$score
    )
  }
  
  # # Convert the data frames to GRanges objects
  gr_pos <- df_to_granges(regions_pos)
  gr_neg <- df_to_granges(regions_neg)
  #
  # # Combine the GRanges objects into a GRangesList
  gr_list <- GRangesList(PMpos_enriched = gr_pos, PMneg_enriched = gr_neg)
  # lapply(gr_list, length)
  
  #gr_list <- df_to_granges(regions)
  #gr_list <- df_to_granges(regions_neg)
  
  ## Prepare metatracks
  meta_tracks <- data_prepared$meta
  # drop rows where condition == PMneg
  #meta_tracks <- meta_tracks[meta_tracks$condition != "PMpos", ]
  meta_tracks$scaling_factor <- 1
  meta_tracks$sample_id <- rownames(meta_tracks)
  meta_tracks <- meta_tracks[, c(ncol(meta_tracks), 1:(ncol(meta_tracks)-1))]
  meta_tracks$sample_id <- as.factor(meta_tracks$sample_id)
  rownames(meta_tracks) <- NULL
  
  meta_tracks <- meta_tracks %>%
    dplyr::mutate(bam_fn = paste0("/home/turagapm/peritoneal_processing/trimmed_data_bam_9_18_2023/", sample_id, "_bowtie2.bam"))
  
  # Create design data frame
  metagene_design <- data.frame(sample_name = meta_tracks$sample_id,
                                Samples = meta_tracks$bam_fn,
                                PMpos = ifelse(meta_tracks$condition == "PM_positive", 1, 0),
                                PMneg = ifelse(meta_tracks$condition == "PM_negative", 1, 0)
                                )
  # set rownames to sample_name
  rownames(metagene_design) <- metagene_design$sample_name
  # drop sample_name column
  metagene_design <- metagene_design[, -1]
  head(metagene_design)
  
  mg_bam_files <- meta_tracks$bam_fn
  
  mg <- metagene2::metagene2$new(regions = gr_list,
                      bam_files = mg_bam_files,
                      assay='chipseq',
                      paired_end=TRUE,
                      cores=45,
                      bin_count=200,
                      region_mod="separate",
                      padding_size = 100
                      #extend_reads=180
                      )
  
  # We can then plot coverage over those regions across all bam files.
  metagene_plots_nest <- mg$produce_metagene(title = "Demo metagene plot",design=metagene_design, normalization = "RPM",alpha=1,facet_by=~region)
  metagene2_heatmap(mg)
}

### BOXPLOTS OF ENRCIHMENT ----
# if not healthy_v_diseased
#if(healthy_v_diseased != TRUE){
if(TRUE){
  ### For plotting healthy with PMpos and PM neg in heatmap
  meta_name_all <- meta_name
  counts_name_all <- counts_name
  counts_data_all <- read.csv(counts_name_all,row.names = 1)
  #columns in counts_data_all
  meta_all <-  read.csv(meta_name_all,row.names = 1)
  
  # drop all columns in meta where condition is not equal to healthy or PM_positive
  meta_all <- meta_all[meta_all$condition %in% c("PM_positive", "PM_negative","HEALTHY"),] #, ,"HEALTHY""PM_negative"
  
  # drop columns from counts_data_all that are not in meta_all
  counts_data_all <- counts_data_all[, colnames(counts_data_all) %in% rownames(meta_all)]
  
  # Create DESeq2 object without specifying a design
  dds_all <- DESeqDataSetFromMatrix(countData = counts_data_all, 
                                    colData = meta_all, 
                                    design = ~ 1)  # Using '~ 1' as a placeholder design
  
  # Estimate size factors and normalize
  dds_all <- estimateSizeFactors(dds_all)
  normalized_counts_all <- counts(dds_all, normalized=TRUE)
  normalized_counts_all_tb <- data.frame(normalized_counts_all) %>% 
    rownames_to_column(var="gene") %>% 
    as_tibble()
  
} else {
  meta_all <- data_prepared$meta
  if (TRUE) {
    normalized_counts_all_tb <- data.frame(data_prepared$normcounts_data) %>% 
      rownames_to_column(var="gene") %>% 
      as_tibble()
  }
}

if (TRUE) {
  # Prepare the data
  annotation <- meta_all %>%
    dplyr::select(condition, primary_present, ovr_histopath, chemo_6weeks, primary_site)
  
  # Convert to data.frame to avoid tibble rownames warning
  normalized_counts_all_tb_box <- as.data.frame(normalized_counts_all_tb)
  rownames(normalized_counts_all_tb_box) <- normalized_counts_all_tb_box$gene
  normalized_counts_all_tb_box$gene <- NULL
  
  log_counts_data_mat <- data.matrix(normalized_counts_all_tb_box, rownames.force = NA)
  log_counts_data_mat <- t(scale(t(log_counts_data_mat), scale = FALSE, center = FALSE))
  
  # Prepare data for boxplot
  boxplot_data <- as.data.frame(log_counts_data_mat)
  boxplot_data <- boxplot_data[rownames(boxplot_data) %in% stability_df_grouped$feature, ]
  boxplot_data$gene <- rownames(boxplot_data)
  boxplot_data_long <- tidyr::pivot_longer(boxplot_data, -gene, names_to = "sample", values_to = "normalized_count")
  
  # Add condition information
  boxplot_data_long$condition <- annotation$condition[match(boxplot_data_long$sample, rownames(annotation))]
  
  # Convert the condition column to a factor
  boxplot_data_long$condition <- factor(boxplot_data_long$condition)
  
  # Determine the unique conditions
  unique_conditions <- levels(boxplot_data_long$condition)
  
  # Determine the baseline condition (assuming it's the last one in alphabetical order)
  baseline_condition <- sort(unique_conditions, decreasing = TRUE)[2]
  alternative_condition <- setdiff(unique_conditions, baseline_condition)[1]
  
  # Add gene group information based on conditions
  gene_groups <- ifelse(stability_df_grouped$mean < 0, 
                        paste0(baseline_condition, "_enriched"),
                        paste0(alternative_condition, "_enriched"))
  names(gene_groups) <- stability_df_grouped$feature
  boxplot_data_long$gene_group <- gene_groups[boxplot_data_long$gene]
  
  # Calculate fold change values with different baselines for each gene group
  boxplot_data_long <- boxplot_data_long %>%
    group_by(gene, gene_group) %>%
    mutate(
      baseline_mean = case_when(
        gene_group == paste0(baseline_condition, "_enriched") ~ mean(normalized_count[condition == alternative_condition], na.rm = TRUE),
        gene_group == paste0(alternative_condition, "_enriched") ~ mean(normalized_count[condition == baseline_condition], na.rm = TRUE),
        TRUE ~ NA_real_
      ),
      fold_change = (normalized_count / baseline_mean - 1) * 100  # Convert to percentage
    ) %>%
    ungroup()
  
  # Remove rows with non-finite fold change values
  boxplot_data_long <- boxplot_data_long %>% filter(is.finite(fold_change))
  
  # Define colors for conditions as per your specifications
  condition_colors <- c(
    'DISEASED' = color_tables$condition_table[["DISEASED"]],      # Red
    'PM_positive' = color_tables$condition_table[["PM_positive"]],   # Red
    'PM_negative' = color_tables$condition_table[["PM_negative"]],   # Blue
    'HEALTHY' = color_tables$condition_table[["HEALTHY"]]        # Dark grey
  )
  
  
  # Ensure all conditions have assigned colors
  missing_conditions <- setdiff(unique_conditions, names(condition_colors))
  if (length(missing_conditions) > 0) {
    warning("The following conditions are missing colors: ", paste(missing_conditions, collapse = ", "))
  }
  
  # Set y-axis limits
  y_min <- -65  # -65% change
  y_max <- 90
  
  # Update the test function to return the full test object
  test_meaningful_difference <- function(x, y) {
    wilcox_result <- wilcox.test(x, y, alternative = "two.sided")
    return(wilcox_result)
  }
  
  # Create a full test function to collect additional statistics
  test_meaningful_difference_full <- function(x, y) {
    wilcox_result <- wilcox.test(x, y, alternative = "two.sided")
    return(list(
      p.value = wilcox_result$p.value,
      statistic = wilcox_result$statistic,
      method = wilcox_result$method
    ))
  }
  
  # Calculate the number of regions (genes) in each gene group
  n_regions <- boxplot_data_long %>%
    group_by(gene_group) %>%
    summarise(n_regions = n_distinct(gene), .groups = 'drop')
  
  # Create a named vector for easier lookup
  n_regions_vector <- setNames(n_regions$n_regions, n_regions$gene_group)
  
  # Create a custom labeller function to include (n=n_regions) in subtitles
  gene_group_labeller <- function(value) {
    n <- n_regions_vector[value]
    return(paste(value, " (n=", n, " regions)", sep = ""))
  }
  
  # Define comparisons
  comparisons <- combn(levels(as.factor(boxplot_data_long$condition)), 2, simplify = FALSE)
  
  # Calculate mean values for each condition and gene group
  mean_values <- boxplot_data_long %>%
    group_by(condition, gene_group) %>%
    summarise(mean_fold_change = mean(fold_change), .groups = 'drop')
  
  # Create a dataframe to store statistical test results
  stat_results_list <- list()
  
  # Loop over each gene group and perform statistical tests
  for (current_gene_group in unique(boxplot_data_long$gene_group)) {
    # Subset data to the current gene group
    data_subset <- boxplot_data_long %>% filter(gene_group == current_gene_group)
    
    # Get all pairwise comparisons of conditions
    comparisons <- combn(levels(data_subset$condition), 2, simplify = FALSE)
    
    # Perform tests for each comparison
    for (comp in comparisons) {
      condition1 <- comp[1]
      condition2 <- comp[2]
      
      # Get data for the two conditions
      data_condition1 <- data_subset %>% filter(condition == condition1)
      data_condition2 <- data_subset %>% filter(condition == condition2)
      
      # Get fold change values
      fold_change1 <- data_condition1$fold_change
      fold_change2 <- data_condition2$fold_change
      
      # Perform the test
      test_result <- test_meaningful_difference_full(fold_change1, fold_change2)
      
      # Calculate mean fold change for each condition
      mean_fold_change1 <- mean(fold_change1, na.rm = TRUE)
      mean_fold_change2 <- mean(fold_change2, na.rm = TRUE)
      
      # Calculate sample sizes
      n_condition1 <- length(fold_change1)
      n_condition2 <- length(fold_change2)
      
      # Collect the results
      stat_result <- data.frame(
        gene_group = current_gene_group,
        condition1 = condition1,
        condition2 = condition2,
        test_used = test_result$method,
        p_value = test_result$p.value,
        statistic = as.numeric(test_result$statistic),
        mean_fold_change_condition1 = mean_fold_change1,
        mean_fold_change_condition2 = mean_fold_change2,
        n_condition1 = n_condition1,
        n_condition2 = n_condition2
      )
      
      # Append to the list
      stat_results_list[[length(stat_results_list) + 1]] <- stat_result
    }
  }
  
  # Combine all results into a single dataframe
  stat_sign_gene_grp <- do.call(rbind, stat_results_list)
  
  # Print the statistical results dataframe
  print(stat_sign_gene_grp)
  
  # Create boxplot with scatter and mean values
  boxplot_enriched <- ggplot(boxplot_data_long, aes(x = condition, y = fold_change, fill = condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    facet_wrap(~gene_group, scales = "free_x", ncol = 1, labeller = as_labeller(gene_group_labeller)) +
    theme_minimal(base_size = 20) +
    labs(title = "Fold Change of 5hmC",
         x = "Condition",
         y = "Fold Change (%)") +
    geom_signif(
      comparisons = comparisons,
      step_increase = 0.02,
      map_signif_level = FALSE,
      test = test_meaningful_difference,
      y_position = seq(60, 60, length.out = length(comparisons)),
      tip_length = 0.0025,
    ) +
    #geom_point(aes(y = fold_change), position = position_jitterdodge(jitter.width = 0.25), alpha = 0.2, size = 0.5) +
    scale_fill_manual(values = condition_colors) +
    theme(
      legend.position = "right",
      strip.text = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA)) +
    scale_y_continuous(
      limits = c(y_min, y_max),
      breaks = seq(y_min - y_min %% 50, y_max, by = 25),  # Major breaks at multiples of 50 including 0
      labels = function(x) paste0(x, "%")
    ) +
    # Add mean values as text
    geom_text(data = mean_values, 
              aes(x = condition, y = y_max, label = sprintf("%.1f%%", mean_fold_change)),
              vjust = 5, size = 3, fontface = "bold")
  
  # Display the plot
  print(boxplot_enriched)
  
  # Save boxplot
}


### PLOT INDIVIDUAL ROC PLOTS FOR EACH SIG GENE ----
# Assuming 'data' is your full dataset
if(TRUE){
  source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
  plots <- plot_roc_curves_per_feature_nestedcv_logreg(data_prepared$normcounts_data, 
                                                       data_prepared$meta, 
                                                       stability_df_grouped, 
                                                       "ROC Curves for Individual Features (Nested CV)",
                                                       n_outer_folds = 10)
  
  # To view the negative coefficients plot (in Purple-Blue colors)
  plots$negative_plot
  
  # To view the positive coefficients plot (in Yellow-Orange-Red colors)
  plots$positive_plot
}

### GRIDSEARCH ----
if(FALSE){
  # Load required libraries
  library(future)
  library(future.apply)
  library(glmnet)
  library(caret)
  library(ROCR)
  library(ggplot2)
  library(reshape2)
  library(parallel)
  library(dplyr)
  library(gridExtra)
  library(pROC)
  library(tidyr)
  
  # Set up parallel processing
  plan(list(tweak(multisession, workers = 10)))
  
  # Define linked parameter vectors
  # Example: each index corresponds to a specific set of linked parameters
  n_outer_folds_values <- c(10)
  n_inner_folds_values <- c(10)
  training_fraction_values <- c(100)
  
  # Check that lengths match if more than one value is provided
  param_lengths <- c(length(n_outer_folds_values), length(n_inner_folds_values), length(training_fraction_values))
  if (any(param_lengths > 1) && (length(unique(param_lengths[param_lengths > 1])) != 1)) {
    stop("When providing multiple values for n_outer_folds, n_inner_folds, and training_fraction, their lengths must match.")
  }
  
  # Other parameters
  alpha_min_values <- c(0.1)
  num_top_genes_values <- c(10,15,25,50,100,150,250,500,1000,2000,5000,10000,20000)
  num_trees_values <- c(10000)
  max_depth_values <- c(30)
  min_node_size_values <- c(20)
  mtry_factor_values <- c(0.2)
  seed_values <- 1:30
  
  # Create parameter grid by expanding over alpha_min, num_top_genes, and seeds,
  # while using the linked triplets of (n_outer_folds, n_inner_folds, training_fraction)
  # Note: For each triple (outer, inner, fraction), we combine with all other param combinations.
  
  param_grid_list <- list()
  for (j in seq_along(n_outer_folds_values)) {
    for (a in alpha_min_values) {
      for (g in num_top_genes_values) {
        for (t in num_trees_values) {
          for (d in max_depth_values) {
            for (mns in min_node_size_values) {
              for (mf in mtry_factor_values) {
                for (s in seed_values) {
                  param_grid_list[[length(param_grid_list) + 1]] <- data.frame(
                    alpha_min = a,
                    num_top_genes = g,
                    n_inner_folds = n_inner_folds_values[j],
                    n_outer_folds = n_outer_folds_values[j],
                    training_fraction = training_fraction_values[j],
                    num_trees = t,
                    max_depth = d,
                    min_node_size = mns,
                    mtry_factor = mf,
                    seed = s,
                    stringsAsFactors = FALSE
                  )
                }
              }
            }
          }
        }
      }
    }
  }
  
  param_grid <- do.call(rbind, param_grid_list)
  
  ### PERFORM GRID SEARCH WITH NESTED CV
  start_time <- Sys.time()
  print(paste("Started at: ", start_time))
  
  results_list <- future_lapply(1:nrow(param_grid), function(i) {
    current_params <- param_grid[i, ]
    cat(paste("Starting on parameter set:", i, "\n"))
    
    # Set the random seed for reproducibility
    set.seed(current_params$seed)
    
    # PARTITION DATA with linked training_fraction
    partitions <- partition_data(
      data_prepared$meta,
      interaction_params,
      current_params$training_fraction,
      validation_fraction,
      holdout_fraction,
      current_params$seed,
      debug_progress
    )
    
    # ORDER AND FILTER DATA
    ordered_filtered_counts <- order_and_filter_data(
      partitions$training_meta,
      partitions$validation_meta,
      partitions$holdout_meta,
      data_prepared$counts_data,
      median_percentile_cutoff,
      variance_percentile_cutoff,
      debug_progress
    )
    
    # NORMALIZE COUNTS
    normalized_counts <- normalize_counts(
      ordered_filtered_counts$training,
      ordered_filtered_counts$validation,
      ordered_filtered_counts$holdout,
      combat_norm,
      length_norm,
      debug_progress
    )
    
    meta_training <- partitions$training_meta
    meta_validation <- partitions$validation_meta
    meta_holdout <- partitions$holdout_meta
    
    norm_counts_data_training <- normalized_counts$training_normcounts
    norm_counts_data_validation <- if (validation_fraction == 0) NULL else normalized_counts$validation_normcounts
    norm_counts_data_holdout <- normalized_counts$holdout_normcounts
    
    conditions_vector_training <- create_conditions_vector(meta_training, "condition")
    conditions_vector_validation <- if (validation_fraction == 0) NULL else create_conditions_vector(meta_validation, "condition")
    conditions_vector_holdout <- if (validation_fraction == 0) NULL else create_conditions_vector(meta_holdout, "condition")
    
    conditions_vector_training_fac <- as.factor(make.names(conditions_vector_training))
    conditions_vector_validation_fac <- if (validation_fraction == 0) NULL else as.factor(make.names(conditions_vector_validation))
    conditions_vector_holdout_fac <- if (validation_fraction == 0) NULL else as.factor(make.names(conditions_vector_holdout))
    
    # Calculate weights
    class_weights <- weight(conditions_vector_training_fac)
    names(class_weights) <- meta_training$condition
    names(class_weights) <- gsub("PMneg", "X0", names(class_weights))
    names(class_weights) <- gsub("PMpos", "X1", names(class_weights))
    
    unique_weights <- unique(class_weights)
    class_labels <- c("X0", "X1")
    class_weights_vector <- setNames(unique_weights, class_labels)
    
    print("Starting Nested CV...")
    stratification_factor <- interaction(meta_training$stratum)
    stratification_factor_inner <- interaction(meta_training$inner_stratum)
    out_folds <- createFolds(stratification_factor, k = current_params$n_outer_folds)
    in_folds <- lapply(out_folds, function(i) {
      train_y <- stratification_factor_inner[-i]
      caret::createFolds(train_y, k = current_params$n_inner_folds)
    })
    
    repeated_fit <- tryCatch({
      nestcv.glmnet(
        conditions_vector_training_fac,
        t(normalized_counts$training_normcounts),
        family = "binomial",
        type.measure = "auc",
        standardize = TRUE,
        outer_folds = out_folds,
        inner_folds = in_folds,
        alphaSet = seq(current_params$alpha_min, 1, by = 0.05),
        cv.cores = 10,
        n_outer_folds = current_params$n_outer_folds,
        n_inner_folds = current_params$n_inner_folds,
        outer_method = "cv",
        finalCV = FALSE,
        verbose = TRUE,
        min_1se = 1,
        keep = TRUE,
        outer_train_predict = TRUE,
        filterFUN = ranger_filter,
        filter_options = list(
          type = "index",
          num.trees = current_params$num_trees,
          max.depth = current_params$max_depth,
          min.node.size = current_params$min_node_size,
          nfilter = current_params$num_top_genes,
          mtry = floor(current_params$mtry_factor * nrow(norm_counts_data_holdout)),
          class.weights = class_weights_vector
        )
      )
    }, error = function(e) {
      message("An error occurred during nested CV: ", e$message)
      return(NULL)
    })
    
    #print("Finished nested CV...")
    final_fit <- NULL
    #final_fit <- repeated_fit$final_fit
    
    # If final fit is NULL, return NULL
    # if (is.null(final_fit)) {
    #   return(NULL)
    # }
    
    selected_predictors <- NULL
    try({
      selected_predictors <- final_fit$glmnet.fit$beta@Dimnames[[1]]
    }, silent = TRUE)
    if (is.null(selected_predictors)) {
      try({
        selected_predictors <- final_fit$beta@Dimnames[[1]]
      }, silent = TRUE)
    }
    
    if (is.null(selected_predictors)) {
      cat("Error: Failed to extract Dimnames\n")
      pred_val <- NULL
      roc_perf_val <- NULL
      auc_val <- NULL
    }
    
    if (!is.null(norm_counts_data_validation) && !is.null(conditions_vector_validation)) {
      cat("Making predictions for parameter set:", i, "\n")
      
      norm_counts_data_validation_top <- norm_counts_data_validation[selected_predictors, , drop = FALSE]
      common_predictors <- intersect(selected_predictors, rownames(norm_counts_data_validation))
      if (length(selected_predictors) == length(common_predictors)) {
        cat("Matched selected_preds and cols in norm_counts")
        probabilities_val <- predict(final_fit, newx = t(norm_counts_data_validation_top), s = repeated_fit$final_param[1], type = "response")
        pred_val <- prediction(probabilities_val, conditions_vector_validation)
        roc_perf_val <- performance(pred_val, "tpr", "fpr")
        auc_val <- as.numeric(performance(pred_val, measure = "auc")@y.values)
      } else {
        cat("Mismatched selected_preds and cols in norm_counts")
        pred_val <- NULL
        roc_perf_val <- NULL
        auc_val <- NULL
      }
    } else {
      pred_val <- NULL
      roc_perf_val <- NULL
      auc_val <- NULL
    }
    
    cat("Saving results for parameter set:", i, "\n")
    
    # Compute sample counts
    total_training_samples <- nrow(meta_training)
    test_fold_sizes <- sapply(out_folds, length)
    avg_test_size <- mean(test_fold_sizes)
    
    ### NOTE: NEED TO UPDATE TO SUBTRACT TEST SIZE FROM TRAINING SIZE
    ### FOR PAPER, CHANGED ILLUSTRATOR DIRECTLY.
    
    result <- list(
      param_set = i,
      seed = current_params$seed,
      alpha_min = current_params$alpha_min,
      num_top_genes = current_params$num_top_genes,
      n_inner_folds = current_params$n_inner_folds,
      n_outer_folds = current_params$n_outer_folds,
      training_fraction = current_params$training_fraction,
      num_trees = current_params$num_trees,
      max_depth = current_params$max_depth,
      min_node_size = current_params$min_node_size,
      mtry_factor = current_params$mtry_factor,
      num_features = nrow(repeated_fit$final_coef),
      auc_tra_inner = signif(pROC::auc(innercv_roc(repeated_fit)), 3),
      auc_tra_outer = signif(pROC::auc(repeated_fit$roc), 3),
      auc_val = auc_val,
      alpha_lambda = repeated_fit$final_param,
      coefficients = repeated_fit$final_coef,
      roc_data_tra_inner = innercv_roc(repeated_fit),
      roc_data_tra_outer = repeated_fit$roc,
      roc_data_val = roc_perf_val,
      pred_val = pred_val,
      repeated_fit = repeated_fit,
      final_fit = final_fit,
      meta_ordered = ordered_filtered_counts$training$meta_ordered,
      total_training_samples = total_training_samples,
      avg_test_size = avg_test_size
    )
    
    return(result)
  }, future.packages = c("caret", "dplyr", "limma", "DESeq2", "sva", "doParallel", "foreach", "progress", "glmnet", "ROCR", "future.apply", "BiocParallel"), future.seed = TRUE)
  
  end_time <- Sys.time()
  print(paste("Total Operation Time: ", end_time - start_time))
  
  # Remove NULL elements
  results_list <- results_list[!sapply(results_list, is.null)]
  
  results_df <- do.call(rbind, lapply(results_list, function(x) {
    data.frame(
      param_set = x$param_set,
      seed = x$seed,
      alpha_min = x$alpha_min,
      num_top_genes = x$num_top_genes,
      n_inner_folds = x$n_inner_folds,
      n_outer_folds = x$n_outer_folds,
      training_fraction = x$training_fraction,
      num_trees = x$num_trees,
      max_depth = x$max_depth,
      min_node_size = x$min_node_size,
      mtry_factor = x$mtry_factor,
      num_features = x$num_features,
      auc_tra_inner = x$auc_tra_inner,
      auc_tra_outer = x$auc_tra_outer,
      auc_val = ifelse(is.null(x$auc_val), NA, x$auc_val),
      total_training_samples = ifelse(is.null(x$total_training_samples), NA, x$total_training_samples),
      avg_test_size = ifelse(is.null(x$avg_test_size), NA, x$avg_test_size),
      stringsAsFactors = FALSE
    )
  }))
  
  results_summary <- results_df %>%
    dplyr::group_by(alpha_min, num_top_genes, n_inner_folds, n_outer_folds, training_fraction, num_trees, max_depth, min_node_size, mtry_factor) %>%
    summarize(
      avg_auc_tra_inner = mean(auc_tra_inner, na.rm = TRUE),
      avg_auc_tra_outer = mean(auc_tra_outer, na.rm = TRUE),
      sd_auc_tra_inner = sd(auc_tra_inner, na.rm = TRUE),
      sd_auc_tra_outer = sd(auc_tra_outer, na.rm = TRUE),
      avg_num_features = mean(num_features, na.rm = TRUE),
      sd_num_features = sd(num_features, na.rm = TRUE),
      avg_training_samples = mean(total_training_samples, na.rm = TRUE),
      avg_test_samples = mean(avg_test_size, na.rm = TRUE),
      num_seeds = n(),
      .groups = 'drop'
    )
  
  # Save results_summary
  output_folder <- file.path(output_folder_ver, "summary_statistics")
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  write.csv(results_summary, file = file.path(output_folder, paste0("results_summary_", file_version, ".csv")), row.names = FALSE)
  
  # Find best params
  best_params <- results_summary %>%
    dplyr::arrange(desc(avg_auc_tra_outer)) %>%
    dplyr::slice(1)
  
  print("Best performing parameters:")
  print(best_params)
  
  # Identify varying parameters
  param_columns <- c("alpha_min", "num_top_genes", "n_inner_folds", "n_outer_folds", 
                     "training_fraction", "num_trees", "max_depth", "min_node_size", "mtry_factor")
  
  varying_params <- sapply(results_summary[, param_columns], function(x) length(unique(x)) > 1)
  varying_param_names <- names(varying_params)[varying_params]
  
  print("Varying parameters:")
  print(varying_param_names)
  
  if (length(varying_param_names) == 0) {
    print("No varying parameters found.")
  } else {
    # Create individual plots for varying parameters vs avg_auc_tra_outer
    plot_param <- function(data, x_var, y_var = "avg_auc_tra_outer") {
      ggplot(data, aes_string(x = x_var, y = y_var)) +
        geom_point() +
        theme_minimal() +
        labs(title = paste(x_var, "vs", y_var),
             x = x_var,
             y = y_var)
    }
    
    plots <- lapply(varying_param_names, function(param) {
      plot_param(results_summary, param)
    })
    
    # Save individual plots
    output_folder <- file.path(output_folder_ver, "parameter_plots")
    if (!dir.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
    
    for (i in seq_along(plots)) {
      plot_name <- varying_param_names[i]
      ggsave(filename = file.path(output_folder, paste0("plot_", plot_name, "_", file_version, ".png")),
             plot = plots[[i]], width = 8, height = 6, dpi = 300)
      ggsave(filename = file.path(output_folder, paste0("plot_", plot_name, "_", file_version, ".svg")),
             plot = plots[[i]], width = 8, height = 6)
    }
    
    # If at least two varying parameters, create heatmap of top two correlated
    if (length(varying_param_names) >= 2) {
      correlations <- sapply(results_summary[, varying_param_names], function(x) cor(x, results_summary$avg_auc_tra_outer, use = "complete.obs"))
      correlation_df <- data.frame(
        Parameter = varying_param_names,
        Correlation = correlations
      )
      
      print("Correlations with avg_auc_tra_outer:")
      print(correlation_df)
      
      # top two correlated parameters
      top_params <- names(sort(abs(correlations), decreasing = TRUE))[1:2]
      
      heatmap_data <- results_summary %>%
        group_by(!!sym(top_params[1]), !!sym(top_params[2])) %>%
        summarize(mean_auc = mean(avg_auc_tra_outer, na.rm = TRUE), .groups = 'drop')
      
      heatmap_plot <- ggplot(heatmap_data, aes_string(x = top_params[1], y = top_params[2], fill = "mean_auc")) +
        geom_tile() +
        scale_fill_gradient(low = "white", high = "red") +
        theme_minimal() +
        labs(title = paste("Heatmap of Mean AUC for", top_params[1], "and", top_params[2]),
             x = top_params[1],
             y = top_params[2],
             fill = "Mean AUC")
      
      ggsave(filename = file.path(output_folder, paste0("heatmap_", file_version, ".png")),
             plot = heatmap_plot, width = 8, height = 6, dpi = 300)
      ggsave(filename = file.path(output_folder, paste0("heatmap_", file_version, ".svg")),
             plot = heatmap_plot, width = 8, height = 6)
    }
  }
  
  # Summary stats
  summary_stats <- results_summary %>%
    summarize(
      mean_avg_auc = mean(avg_auc_tra_outer, na.rm = TRUE),
      median_avg_auc = median(avg_auc_tra_outer, na.rm = TRUE),
      min_avg_auc = min(avg_auc_tra_outer, na.rm = TRUE),
      max_avg_auc = max(avg_auc_tra_outer, na.rm = TRUE),
      sd_avg_auc = sd(avg_auc_tra_outer, na.rm = TRUE)
    )
  
  print("Summary statistics for average AUC:")
  print(summary_stats)
  
  # Save summary stats
  write.csv(summary_stats, file = file.path(output_folder, paste0("summary_stats_", file_version, ".csv")), row.names = FALSE)
  
  # Top 5 parameter combos
  top_5_params <- results_summary %>%
    dplyr::arrange(desc(avg_auc_tra_outer)) %>%
    dplyr::select(all_of(param_columns), avg_auc_tra_outer) %>%
    dplyr::slice_head(n = 5)
  
  print("Top 5 parameter combinations:")
  print(top_5_params)
  
  write.csv(top_5_params, file = file.path(output_folder, paste0("top_5_params_", file_version, ".csv")), row.names = FALSE)
  
  # Plot AUC vs features
  auc_vs_features_plot <- ggplot(results_summary, aes(x = avg_num_features, y = avg_auc_tra_outer)) +
    geom_point() +
    geom_text(
      aes(label = paste0("α=", alpha_min, ", Features=", num_top_genes)),
      hjust = 0.5,
      vjust = -1,
      size = 3,
      check_overlap = TRUE
    ) +
    theme_minimal() +
    labs(
      title = "Average AUC vs. Average Number of Features",
      x = "Average Number of Features",
      y = "Average AUC (Outer CV)"
    ) +
    scale_x_continuous(trans = "log2")
  
  ggsave(
    filename = file.path(output_folder, paste0("auc_vs_features_", file_version, ".png")),
    plot = auc_vs_features_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(output_folder, paste0("auc_vs_features_", file_version, ".svg")),
    plot = auc_vs_features_plot,
    width = 8,
    height = 6
  )
  
  # Create the requested final plot:
  # primary y axis = AUC
  # secondary y axis = number of samples
  # x axis = n_outer_folds
  max_samples <- max(results_summary$avg_training_samples, results_summary$avg_test_samples, na.rm = TRUE)
  scale_factor <- max_samples
  
  sample_data <- results_summary %>%
    dplyr::select(n_outer_folds, training_fraction, avg_auc_tra_outer, avg_training_samples, avg_test_samples) %>%
    pivot_longer(cols = c("avg_training_samples", "avg_test_samples"), 
                 names_to = "sample_type", values_to = "sample_count")
  
  # set avg_training_samples to avg_training_samples - avg_test_samples
  sample_data$sample_count <- ifelse(sample_data$sample_type == "avg_training_samples", 
                                     sample_data$sample_count - sample_data$sample_count[sample_data$sample_type == "avg_test_samples"],
                                     sample_data$sample_count)
  
  final_plot <- ggplot(sample_data, aes(x = factor(training_fraction ), group = sample_type)) +
    geom_col(aes(y = sample_count / scale_factor, fill = sample_type),
             position = position_dodge(width = 0.9), width = 0.8) +
    geom_point(data = distinct(sample_data[, c("training_fraction", "avg_auc_tra_outer")]),
               aes(x = factor(training_fraction), y = avg_auc_tra_outer), 
               color = "black", size = 3, shape = 19, inherit.aes = FALSE) +
    theme_minimal() +
    labs(
      title = "AUC vs. Training/Test Sample Counts by training_fraction",
      x = "N Outer Folds",
      y = "AUC (Primary Axis)"
    ) +
    scale_fill_manual(values = c("avg_training_samples" = "blue", "avg_test_samples" = "red"),
                      labels = c("Training Samples", "Test Samples")) +
    scale_y_continuous(
      sec.axis = sec_axis(~ . * scale_factor, name = "Number of Samples (Secondary Axis)")
    ) +
    # Add labels to bars
    geom_text(aes(y = sample_count / scale_factor, label = round(sample_count)),
              position = position_dodge(width = 0.9), vjust = -0.5, size = 3, color = "black") +
    # Add labels to points
    geom_text(data = distinct(sample_data[, c("training_fraction", "avg_auc_tra_outer")]),
              aes(x = factor(training_fraction), y = avg_auc_tra_outer, label = round(avg_auc_tra_outer, 3)),
              vjust = -1.0, size = 3, color = "black", inherit.aes = FALSE)
  
  # Save the final plot
  ggsave(filename = file.path(output_folder, paste0("AUC_vs_Samples_", file_version, ".png")),
         plot = final_plot, width = 10, height = 7, dpi = 300)
  ggsave(filename = file.path(output_folder, paste0("AUC_vs_Samples_", file_version, ".svg")),
         plot = final_plot, width = 10, height = 7)
  
  # Step 1: Generate a reshaped data frame with training_fraction as columns
  sample_info <- results_summary %>%
    dplyr::select(training_fraction, avg_training_samples, avg_test_samples) %>%
    distinct() %>%
    arrange(training_fraction) %>%
    # Reshape to long format
    pivot_longer(
      cols = c(avg_training_samples, avg_test_samples),
      names_to = "Sample_Type",
      values_to = "Count"
    ) %>%
    # Clean up Sample_Type names for better readability
    mutate(Sample_Type = recode(Sample_Type, 
                                avg_training_samples = "Training Samples",
                                avg_test_samples = "Test Samples")) %>%
    # Reshape back to wide format with training_fraction as columns
    pivot_wider(
      names_from = training_fraction,
      values_from = Count
    )
  
  # Step 2: Create a table grob from the reshaped sample_info
  sample_table <- tableGrob(
    sample_info,
    rows = NULL,  # Remove row names
    theme = ttheme_minimal(base_size = 14)  # Use a minimal theme for clarity
  )
  
  # Step 3: Create a boxplot of AUC vs Training Fraction with Median Labels
  # Calculate median AUC for each training_fraction
  median_info <- results_df %>%
    group_by(training_fraction) %>%
    summarize(median_auc = median(auc_tra_outer, na.rm = TRUE))
  
  # Create the boxplot
  boxplot_auc <- ggplot(results_df, aes(x = factor(training_fraction), y = auc_tra_outer)) +
    geom_boxplot(color = "black", fill= NA,size = 1.2) +  # Added black border for clarity
    theme_minimal() +
    labs(
      title = "AUC Distribution by Training Fraction",
      x = "Training Fraction",
      y = "AUC"
    ) +
    # Increase font size
    theme(
      text = element_text(size = 16),
      plot.title = element_text(hjust = 0.5)  # Center the title
    ) +
    # Add median labels
    geom_text(
      data = median_info,
      aes(x = factor(training_fraction), y = median_auc, label = round(median_auc, 2)),
      vjust = -0.5,  # Position above the median point
      color = "black",
      size = 5
    ) +
    # Adjust y-axis limits to accommodate labels
    ylim(NA, max(results_df$auc_tra_outer, na.rm = TRUE) * 1.05)
  
  # Step 4: Arrange the boxplot and the table vertically with aligned columns
  final_plot <- grid.arrange(
    boxplot_auc,
    sample_table,
    ncol = 1,         # Stack vertically
    heights = c(4, 1)  # Allocate more space to the boxplot
  )
  
  # Save the final plot
  ggsave(filename = file.path(output_folder, paste0("AUC_boxplot_with_table_", file_version, ".png")),
         plot = final_plot, width = 4.5, height = 7, dpi = 300)
  
  ggsave(filename = file.path(output_folder, paste0("AUC_boxplot_with_table_", file_version, ".svg")),
         plot = final_plot, width = 5, height = 7)
  
  source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
  # First, retrieve all color tables
  color_tables <- create_color_tables()
  
  # For shape, we'll manually map alpha_min (e.g. 0.1 -> 16, 0.3 -> 17).
  # For color, we use the newly added num_top_genes_table in color_tables.
  auc_vs_features_plot <- ggplot(
    results_summary, 
    aes(
      x = avg_num_features, 
      y = avg_auc_tra_outer,
      color = as.factor(num_top_genes),
      shape = as.factor(alpha_min)
    )
  ) +
    geom_point(size = 3) +
    
    # Use log2 scaling on x-axis; set custom breaks at multiples of 50 up to 20000
    scale_x_continuous(
      trans = "log2",
      breaks = c(12,25,50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 20000),
      labels = c("12","25","50", "100", "200", "400", "800", "1.6k", "3.2k", "6.4k", "12.8k", "20k")
    ) +
    
    # Map color and shape manually
    scale_color_manual(
      name = "Num Top Genes",
      values = color_tables$num_top_genes_table
    ) +
    scale_shape_manual(
      name = "Alpha Min",
      values = c("0.1" = 16, "0.3" = 17)
    ) +
    
    # A minimal, “no-background” theme:
    theme_classic(base_size = 14) +
    theme(
      panel.background  = element_blank(), 
      plot.background   = element_blank(),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank()
    ) +
    
    labs(
      title = "Average AUC vs. Average Number of Features",
      x     = "Average Number of Features (log2)",
      y     = "Average AUC (Outer CV)"
    )
  
  ggsave(
    filename = file.path(output_folder, 
                         paste0("auc_vs_features_colorByNumTopGenes_shapeByAlphaMin_", file_version, ".png")),
    plot   = auc_vs_features_plot,
    width  = 8,
    height = 6,
    dpi    = 300
  )
  
  ggsave(
    filename = file.path(output_folder, 
                         paste0("auc_vs_features_colorByNumTopGenes_shapeByAlphaMin_", file_version, ".svg")),
    plot   = auc_vs_features_plot,
    width  = 8,
    height = 6
  )
  
  # 1) Compute the average number of features for each num_top_genes
  df_with_avgs <- results_df %>%
    group_by(num_top_genes) %>%
    mutate(avg_num_features = mean(num_features, na.rm = TRUE)) %>%
    ungroup()
  
  # 2) Compute the median AUC for each avg_num_features
  df_medians <- df_with_avgs %>%
    group_by(avg_num_features) %>%
    summarize(median_auc = median(auc_tra_outer, na.rm = TRUE), .groups = "drop")
  
  # 3) Create the boxplot of AUC vs. average number of features
  boxplot_auc_features <- ggplot(df_with_avgs, aes(x = factor(round(avg_num_features, 2)), y = auc_tra_outer)) +
    geom_boxplot(outlier.shape = NA, color = "black", fill = NA) +
    # Optionally add jittered points if you want to see individual data points:
    # geom_jitter(width = 0.2, alpha = 0.5)
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.line = element_line(color = "black"),
      panel.grid.major = element_line(color = "gray80"),
      panel.grid.minor = element_line(color = "gray90"),
      text = element_text(color = "black", size = 14),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(color = "black", hjust = 0.5)
    ) +
    labs(
      title = "Distribution of Outer AUC by Average Number of Features",
      x = "Average Number of Features (rounded)",
      y = "AUC (Outer CV)"
    ) +
    # 4) Add text labels for median values above each box
    geom_text(
      data = df_medians,
      aes(
        x = factor(round(avg_num_features, 2)),
        y = median_auc,
        label = round(median_auc, 2)
      ),
      vjust = -0.5,
      color = "black",
      fontface = "bold",
      size = 4
    )
  
  # Print or save the plot
  print(boxplot_auc_features)
  
  # 4. Save the plot to file
  ggsave(
    filename = file.path(
      output_folder,
      paste0("auc_boxplot_by_avgNumFeatures_", file_version, ".png")
    ),
    plot   = boxplot_auc_features,
    width  = 8,
    height = 6,
    dpi    = 300
  )
  
  ggsave(
    filename = file.path(
      output_folder,
      paste0("auc_boxplot_by_avgNumFeatures_", file_version, ".svg")
    ),
    plot   = boxplot_auc_features,
    width  = 8,
    height = 6
  )
  
  write.csv(results_df, file = file.path(output_folder, paste0("resultsdf_", file_version, ".csv")), row.names = FALSE)
}

### SAVE FILES ----
if(TRUE){
  ### Save key results and files to output folder
  
  ## Save enrichment boxplots by gene
  tryCatch({
    # Define the path for the new directory
    output_folder <- file.path(output_folder_ver, "boxplot_enrichment_by_gene")
    
    # Create the directory if it does not exist
    if (!dir.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
    
    # Save boxplots by gene
    ggsave(filename = file.path(output_folder, paste0("boxplot_negative_", file_version, ".png")), 
           plot = boxplots$negative, width = 20, height = 10)
    ggsave(filename = file.path(output_folder, paste0("boxplot_positive_", file_version, ".png")), 
           plot = boxplots$positive, width = 20, height = 10)
    
    # Save as SVGs
    ggsave(filename = file.path(output_folder, paste0("boxplot_negative_", file_version, ".svg")), 
           plot = boxplots$negative, width = 20, height = 10)
    ggsave(filename = file.path(output_folder, paste0("boxplot_positive_", file_version, ".svg")), 
           plot = boxplots$positive, width = 20, height = 10)
    
    # Save faceted boxplots by gene
    ggsave(filename = file.path(output_folder, paste0("boxplot_negative_faceted_", file_version, ".png")), 
           plot = boxplots$faceted_negative, width = 15, height = 15)
    ggsave(filename = file.path(output_folder, paste0("boxplot_positive_faceted_", file_version, ".png")), 
           plot = boxplots$faceted_positive, width = 15, height = 15)
    
    # Save faceted boxplots as SVGs
    ggsave(filename = file.path(output_folder, paste0("boxplot_negative_faceted_", file_version, ".svg")), 
           plot = boxplots$faceted_negative, width = 15, height = 15)
    ggsave(filename = file.path(output_folder, paste0("boxplot_positive_faceted_", file_version, ".svg")), 
           plot = boxplots$faceted_positive, width = 15, height = 15)
    
    # Save the statistical results to CSV
    write.csv(boxplots$stats_results, 
              file = file.path(output_folder, paste0("statistical_results_", file_version, ".csv")), 
              row.names = FALSE)
    
  }, error = function(e) {
    message("An error occurred, skipping: ", e$message)
  })
  
  # Save key results and files to output folder
  
  # Create a subfolder for results_df
  output_folder <- file.path(output_folder_ver, "results_df")
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  # Save results_df_out
  is_list_column <- sapply(results_df, is.list)
  list_columns <- names(results_df)[is_list_column]
  results_df_out <- results_df %>%
    unnest_wider(all_of(list_columns))
  write.csv(results_df_out, file = file.path(output_folder, paste0("resultsdf_", file_version, ".csv")), row.names = FALSE)
  
  # Create a subfolder for RDS files
  output_folder <- file.path(output_folder_ver, "RDS_files")
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  # Save RDS files
  saveRDS(median_result$repeated_fit, file = file.path(output_folder, paste0("median_model_", file_version, ".rds")))
  saveRDS(training_inner_roc_plot, file = file.path(output_folder, paste0("training_inner_roc_plot", file_version, ".rds")))
  saveRDS(training_outter_roc_plot, file = file.path(output_folder, paste0("training_outter_roc_plot", file_version, ".rds")))
  saveRDS(validation_roc_plot, file = file.path(output_folder, paste0("validation_roc_plot", file_version, ".rds")))
  
  # Create a subfolder for model stability plots
  output_folder <- file.path(output_folder_ver, "model_stability")
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  # Save model stability plots
  plotly::save_image(all_model_stability[[1]], file = file.path(output_folder, paste0("all_model_stability_plot_", file_version, ".svg")), width=800, height=700)
  plotly::save_image(all_model_stability[[1]], file = file.path(output_folder, paste0("all_model_stability_plot_", file_version, ".jpg")), width=800, height=700)
  
  # Create a subfolder for coefficients
  output_folder <- file.path(output_folder_ver, "coefficients")
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  # Save coefficients CSV
  write.csv(all_model_stability_no_thresh[2], file = file.path(output_folder, paste0("coefficients_", file_version, ".csv")), row.names = FALSE)
  
  # Create a subfolder for annotated data
  output_folder <- file.path(output_folder_ver, "annotated_data")
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  # Save annotated data CSVs
  write.csv(annotated_df_enh, file = file.path(output_folder, "annotated_df_enh.csv"), row.names = FALSE)
  write.csv(sig_gene_table, file = file.path(output_folder, "sig_gene_table.csv"), row.names = FALSE)
  
  ### Save all ROC plots as SVG and JPG files
  # Create a subfolder for ROC curves
  output_folder <- file.path(output_folder_ver, "roc_curves")
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  # Save ROC plots
  plotly::save_image(training_inner_roc_plot, file = file.path(output_folder, paste0("training_inner_roc_curves_", file_version, ".svg")))
  plotly::save_image(training_outter_roc_plot, file = file.path(output_folder, paste0("training_outter_roc_curves_", file_version, ".svg")))
  plotly::save_image(training_inner_roc_plot, file = file.path(output_folder, paste0("training_inner_roc_curves_", file_version, ".jpg")))
  plotly::save_image(training_outter_roc_plot, file = file.path(output_folder, paste0("training_outter_roc_curves_", file_version, ".jpg")))
  
  # Save validation ROC plot if it exists
  if (!is.null(validation_roc_plot)) {
    plotly::save_image(validation_roc_plot, file = file.path(output_folder, paste0("validation_roc_curves_", file_version, ".svg")))
    plotly::save_image(validation_roc_plot, file = file.path(output_folder, paste0("validation_roc_curves_", file_version, ".jpg")))
  }
  
  ## Save individual gene ROC plots
  output_folder <- file.path(output_folder_ver, "individual_gene_roc_plots")
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  tryCatch({
    plotly::save_image(plots$positive_plot, file = file.path(output_folder, paste0("positive_plot_", file_version, ".svg")), width = 1200, height = 640)
    plotly::save_image(plots$negative_plot, file = file.path(output_folder, paste0("negative_plot_", file_version, ".svg")), width = 1200, height = 640)
    plotly::save_image(plots$positive_plot, file = file.path(output_folder, paste0("positive_plot_", file_version, ".jpg")), width = 1200, height = 640, dpi = 300)
    plotly::save_image(plots$negative_plot, file = file.path(output_folder, paste0("negative_plot_", file_version, ".jpg")), width = 1200, height = 640, dpi = 300)
  }, error = function(e) {
    message("Individual gene ROC plots do not exist, skipping...: ", e$message)
  })
  
  # Save the configuration file in the main output folder
  config_content <- capture.output({
    cat("Configuration and Parameters\n")
    cat("Counts file: ", counts_name, "\n")
    cat("Normalized counts file: ", normcounts_name, "\n")
    # combat_seq on or off?
    cat("Combat_seq normalization: ", combat_norm, "\n")
    cat("Metadata file: ", meta_name, "\n")
    cat("Training fraction: ", training_fraction, "\n")
    cat("Validation fraction: ", validation_fraction, "\n")
    cat("Holdout fraction: ", holdout_fraction, "\n")
    cat("median_percentile_cutoff: ", median_percentile_cutoff, "\n")
    cat("variance_percentile_cutoff: ", variance_percentile_cutoff, "\n")
    cat("Create significant training counts: ", create_sig_training_counts, "\n")
    cat("Combat normalization: ", combat_norm, "\n")
    cat("Interaction parameters: ", paste(interaction_params, collapse = ", "), "\n")
    cat("Design formula: ", as.character(design_formula), "\n")
    cat("Cutoff type: ", cutoff_type, "\n")
    cat("Number of top genes: ", num_top_genes, "\n")
    cat("Number of folds for validation: ", nfolds_val, "\n")
    cat("Feature occurrence thresholds: ", paste(feat_occ_thresholds, collapse = ", "), "\n")
    
    cat("Standardize: ", standardize_bit, "\n")
    cat("Number of runs: ", number_of_runs, "\n")
    cat("Partition random seed: ", paste(part_rand_seed, collapse = ", "), "\n")
    cat("Number of inner folds: ", n_inner_folds_param, "\n")
    cat("Number of outer folds: ", n_outer_folds_param, "\n")
    cat("min.node.size: ", min_node_size, "\n")
    cat("num.trees: ", num_trees, "\n")
    cat("max.depth: ", max_depth, "\n")
    cat("Rapid mode: ", rapid, "\n")
    cat("Save final files: ", save_final_files, "\n")
    cat("Debug progress: ", debug_progress, "\n")
    cat("Debug data summary: ", debug_data_summary, "\n")
    cat("Random seed: ", random_seed, "\n")
    cat("Output folder: ", output_folder, "\n")
    cat("File version: ", file_version, "\n")
    cat("Output folder version: ", output_folder_ver, "\n")
    cat("Start time: ", as.character(start_time), "\n")
    cat("End time: ", as.character(end_time), "\n")
    # include glmnet settings: alp
    ## glmnet Settings
    cat("### GLMNET SETTINGS: ###")
    cat("Alpha values: ", paste(alpha_values, collapse = ", "), "\n")
    cat("n_inner_folds_param: ", n_inner_folds_param, "\n")
    cat("n_outer_folds_param: ", n_outer_folds_param, "\n")
    cat("num_trees: ", num_trees, "\n")
    cat("max_depth: ", max_depth, "\n")
    cat("min_node_size: ", min_node_size, "\n")
    cat("mtry_factor: ", mtry_factor, "\n")
    
  })
  
  # Write the captured output to the file
  writeLines(config_content, con = file.path(output_folder_ver, "config.txt"))
  
  # Save heatmaps
  output_folder <- file.path(output_folder_ver, "heatmaps")
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  svg(file.path(output_folder, paste0("combined_heatmap_", file_version, ".svg")), width = 9,height=7)
  draw(combined_heatmap)
  dev.off()
  
  png(file.path(output_folder, paste0("combined_heatmap_", file_version, ".png")), width = 9,height=7, units = "in", res = 300)
  draw(combined_heatmap)
  dev.off()
  
  # Save probability boxplots
  output_folder <- file.path(output_folder_ver, "probability_boxplots")
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  ggsave(file.path(output_folder, paste0("combined_prob_boxplot_pm", file_version, ".jpg")), plot = combined_prob_boxplot_pm, width = 4, height = 10, units = "in", dpi = 300)
  ggsave(file.path(output_folder, paste0("combined_prob_boxplot_pm", file_version, ".svg")), plot = combined_prob_boxplot_pm, width = 4, height = 10, units = "in")
  
  ggsave(file.path(output_folder, paste0("combined_prob_boxplot_healthy", file_version, ".jpg")), plot = combined_prob_boxplot_healthy, width = 6, height = 10, units = "in", dpi = 300)
  ggsave(file.path(output_folder, paste0("combined_prob_boxplot_healthy", file_version, ".svg")), plot = combined_prob_boxplot_healthy, width = 6, height = 10, units = "in")
  
  # save healthy_stats and pm_stats as csvs
  write.csv(healthy_stats, file = file.path(output_folder, paste0("healthy_stats_", file_version, ".csv")), row.names = FALSE)
  write.csv(pm_stats, file = file.path(output_folder, paste0("pm_stats_", file_version, ".csv")), row.names = FALSE)
  
  # Save confusion matrix plots
  tryCatch({
    # Define the path for the new directory
    output_folder <- file.path(output_folder_ver, "confusion_matrix")
    
    # Create the directory if it does not exist
    if (!dir.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
    # Save enhanced confusion matrix to PNG and SVG
    ggsave(file.path(output_folder, paste0("enhanced_confusion_matrix_", file_version, ".png")), plot = enhanced_cm_plot , width = 10, height = 10, units = "in", dpi = 300)
    ggsave(file.path(output_folder, paste0("enhanced_confusion_matrix_", file_version, ".svg")), plot = enhanced_cm_plot , width = 10, height = 10, units = "in")
    # Save confusion matrix results as CSV
    write.csv(cm_results_df, file = file.path(output_folder, paste0("cm_results_df_", file_version, ".csv")), row.names = FALSE)
  }, error = function(e) {
    message("Enhanced confusion matrix plot does not exist, skipping...: ", e$message)
  })
  
  # Save boxplot enrichment by group
  output_folder <- file.path(output_folder_ver, "boxplot_enrichment_by_group")
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  tryCatch({
    # Save enrichment by gene group plot
    ggsave(file.path(output_folder, paste0("sig_enriched_plots_", paste(file_version, collapse = "_"), ".png")),
           boxplot_enriched, width = 7.5, height = 12.5, dpi = 300)
    
    ggsave(file.path(output_folder, paste0("sig_enriched_plots_", paste(file_version, collapse = "_"), ".svg")),
           boxplot_enriched, width = 7.5, height = 12.5)
  }, error = function(e) {
    message("sig_enriched_plots do not exist, skipping...: ", e$message)
  })
  
  tryCatch({
    # Save the statistical results dataframe as a CSV file
    write.csv(stat_sign_gene_grp, file = file.path(output_folder, paste0("stat_sign_gene_grp_", paste(file_version, collapse = "_"), ".csv")), row.names = FALSE)
  }, error = function(e) {
    message("Failed to save stat_sign_gene_grp CSV file: ", e$message)
  })
  
  ### Save metagene regions
  output_folder <- file.path(output_folder_ver, "metagene_plots")
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  tryCatch({
    ggsave(file.path(output_folder, paste0("metagene_plot_", file_version, ".svg")), plot = metagene_plots_nest, width = 12, height = 6, units = "in")
  }, error = function(e) {
    message("Metagene plots do not exist, skipping...: ", e$message)
  })
  
  # Save stability_df_grouped as CSV
  output_folder <- file.path(output_folder_ver, "stability_df_grouped")
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  # Group annotated_df_enh by chr, start, and end, keeping all columns
  annotated_df_enh_grouped_temp <- annotated_df_enh %>%
    group_by(chr, start, end) %>%
    dplyr::slice(1) %>%
    ungroup()
  # Save the grouped stability dataframe
  write.csv(annotated_df_enh_grouped_temp, file = file.path(output_folder, paste0("stability_df_grouped_", file_version, ".csv")), row.names = FALSE)

  
}


### Plot number of features
if(FALSE){
  # Load required libraries
  library(ggplot2)
  
  # 1. Read CSVs into R dataframes
  df_hacc <- read.csv(
    "~/5hmC-Pathway-Analysis/Output/repeated_holdout/PMpositive_PMnegative_HMR_HAC_12_21_24_20241221_cond-sex-age__20241221_225341/results_df/resultsdf_PMpositive_PMnegative_HMR_HAC_12_21_24_20241221_cond-sex-age_.csv",
    header = TRUE
  )
  
  df_sparse_pm <- read.csv(
    "~/5hmC-Pathway-Analysis/Output/repeated_holdout/PMpositive_PMnegative_HMR_SPARSE_12_21_24_20241221_cond-sex-age__20241221_225305/results_df/resultsdf_PMpositive_PMnegative_HMR_SPARSE_12_21_24_20241221_cond-sex-age_.csv",
    header = TRUE
  )
  
  df_sparse_cancer <- read.csv(
    "~/5hmC-Pathway-Analysis/Output/repeated_holdout/CANCER_HEALTHY_HMR_12_29_24_20241230_cond-sex-age__20241230_003307/results_df/resultsdf_CANCER_HEALTHY_HMR_12_29_24_20241230_cond-sex-age_.csv",
    header = TRUE
  )
  
  # 2. Add a label to distinguish the dataframes
  df_hacc$method_label        <- "high-accuracy (PM)"
  df_sparse_pm$method_label   <- "sparse (PM)"
  df_sparse_cancer$method_label <- "sparse (cancer)"
  
  # 3. Combine into a single dataframe
  df_all <- rbind(df_hacc, df_sparse_pm, df_sparse_cancer)
  
  # 4. Set the factor levels for `method_label` to define the order
  df_all$method_label <- factor(
    df_all$method_label,
    levels = c("high-accuracy (PM)", "sparse (PM)", "sparse (cancer)")
  )
  
  # 5. Compute medians for labeling
  medians <- aggregate(num_features ~ method_label, data = df_all, FUN = median)
  
  # 6. Create the boxplot
  boxplot_enriched <- ggplot(df_all, aes(x = method_label, y = num_features)) +
    geom_boxplot(
      fill = "transparent",  # No fill for the boxes
      color = "black",       # Black stroke for the boxes and whiskers
      size = 0.8             # Set stroke thickness (adjust as needed)
    ) +
    # Label medians with one decimal
    geom_text(
      data = medians,
      aes(x = method_label, y = num_features, label = sprintf("%.0f", round(num_features, 1))),
      inherit.aes = FALSE,
      color = "black",
      vjust = -0.5,
      size = 3.4
    ) +
    facet_wrap(~method_label, scales = "free", ncol = 3) +
    labs(
      x = NULL,              # or "Method" if you prefer
      y = "Mean model features"
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.line.x      = element_line(color = "black"),
      axis.line.y      = element_line(color = "black"),
      axis.text        = element_text(color = "black"),
      axis.title       = element_text(color = "black"),
      strip.background = element_blank(),   # remove facet strip background
      strip.text       = element_text(color = "black", size = 14),
      panel.border     = element_blank(),
      panel.background = element_blank(),
      plot.background  = element_blank()
    )
  
  # ----------------
  # Saving the plots
  # ----------------
  
  # Define an output folder and file version
  output_folder_ver <- "~/5hmC-Pathway-Analysis/Output/plots" 
  file_version      <- c("num_features", "by_method", "3_methods")
  
  # Create the subfolder for saving
  output_folder <- file.path(output_folder_ver, "boxplot_enrichment_by_group")
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  # Save the plot in both PNG and SVG formats
  tryCatch({
    ggsave(
      file.path(output_folder, paste0("sig_enriched_plots_", paste(file_version, collapse = "_"), ".png")),
      boxplot_enriched,
      width = 6, 
      height = 4.5, 
      dpi = 300
    )
    
    ggsave(
      file.path(output_folder, paste0("sig_enriched_plots_", paste(file_version, collapse = "_"), ".svg")),
      boxplot_enriched,
      width = 6, 
      height = 4.5
    )
  }, error = function(e) {
    message("Error saving plots: ", e$message)
  })
  
  
  
}

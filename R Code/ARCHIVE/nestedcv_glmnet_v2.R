###dataset_group_randomizer.R###
# Purpose: Split a DESeq2-ready file into training, validation, and holdout sets.

# ### INSTALL LIBRARIES IF NECESSARY
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
#install.packages('reticulate')
#reticulate::install_miniconda()
#reticulate::conda_install('r-reticulate', 'python-kaleido')
#reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')
#reticulate::use_miniconda('r-reticulate')

# Debug source this file

### Load Libraries
library(caret)  # for glmnet model training
library(dplyr)
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

source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")


#closeAllConnections()
### CONFIGURATION
## General Settings
setwd("~/5hmC-Pathway-Analysis/")
# Location of raw_data_processing output
root_folder <- "./Output/Raw Data Processing/PM_negative_PM_positive_10per_072624_hmr_nocombat_CRC_ADE/"
filename_trunk <- "PM_negative_PM_positive_DESeq2"
counts_name <- paste0(root_folder, filename_trunk, "_rawcounts.csv")
normcounts_name <- paste0(root_folder, filename_trunk, "_normcounts.csv")
meta_name <- paste0(root_folder, filename_trunk, "_conditions.csv")
file_version <- "hmr_nocombat_07262024_CRCADE_seed10"
output_folder <- "./Output/repeated_holdout/"
save_final_files <- FALSE
debug_progress <- TRUE
debug_data_summary <- FALSE
random_seed <- 10
set.seed(random_seed)


## Partition Settings
training_fraction <- 100
validation_fraction <- 0
holdout_fraction <- 0
median_percentile_cutoff <- 0.25
variance_percentile_cutoff <- 0.25
create_sig_training_counts <- FALSE
combat_norm <- FALSE
length_norm <- FALSE
interaction_params <- c("condition","primary_present")#, "primary_present")#,"primary_site")
inner_interaction_params <- c("condition","primary_present")#, "primary_present")

## DESeq2 Settings
cutoff_type = 2 # 0=padj cutoff, default; 1=lfc & pvalue cutoff; 2 = top stat genes
padj_cutoff = 0.01 # 0.1 default
pvalue_cutoff = 0.001
lfc_cutoff = 0.138 # 0.137504 ~ 10%, 0.263034 ~ 20% change, 0.32 ~ 25%, 0.415 ~ 33%, 0.58 ~ 50% change,
num_top_genes = 25
design_formula <- ~ condition + peritoneal_mets + primary_present #condition
perform_deseq2_only <- 0 # 0 = DESEQ2 + Glmnet , 1 = No Repeated Glmnet

## glmnet Settings
nfolds_val = 5
rep_val = 1
feat_occ_thresholds <- seq(0.7, 1,by = 0.05) # for filter thresh
alpha_values <- seq(0,1, by = 0.01) # for filter thresh
standardize_bit = FALSE # should caret center and standardize before training models?
number_of_runs <- 15
part_rand_seed = seq(1,100, by=1)
#plan(sequential)
plan(list(tweak(multisession, workers = 15))) #sequential, , multicore,multisession
n_inner_folds_param = 10
n_outer_folds_param = 10

rapid = FALSE
if(rapid){
  plan(sequential) #sequential, , multicore,multisession
  num_raw_genes = 1000
  num_top_genes = 500
  feat_occ_thresholds <- seq(0, 0.1, by = 0.1) #c(0.7,0.75,0.8,0.85,0.9,0.95,0.98,1)#
  alpha_values <- seq(0.1, 0.2, by = 0.1)
  standardize_bit = FALSE # should caret center and standardize before training models?
  number_of_runs <- 10
  part_rand_seed = seq(1,2, by=1)
  save_final_files = FALSE
}

#### Main Script
## Subselect metadata file if desired
#subset_conditions <- list(list("primary_site", "HGA"))
subset_conditions <- NULL

### LOAD AND PREPARE META AND COUNTS DATA
data_prepared <- load_and_prepare_data(counts_name, 
                                       normcounts_name, 
                                       meta_name, 
                                       debug_data_summary, 
                                       debug_progress, 
                                       rapid, 
                                       num_raw_genes,
                                       subset_conditions = subset_conditions)

# Add stratum column based on interaction parameters
data_prepared$meta$stratum <- do.call(interaction, data_prepared$meta[interaction_params])
data_prepared$meta$inner_stratum <- do.call(interaction, data_prepared$meta[inner_interaction_params])
conditions <- define_conditions(data_prepared$meta)

### PERFORM NESTEDCV
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
  names(class_weights) <- gsub("PMneg", "X0", names(class_weights))
  names(class_weights) <- gsub("PMpos", "X1", names(class_weights))
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
      cv.cores = 5,
      n_outer_folds = n_outer_folds_param,
      n_inner_folds = n_inner_folds_param,
      outer_method = "cv",
      finalCV = TRUE, 
      verbose = TRUE,
      min_1se = 1,
      keep = TRUE,
      outer_train_predict = TRUE,
      filterFUN = ranger_filter,
      filter_options = list(
        type = "index",
        num.trees = 50000,
        max.depth = 20,
        min.node.size = 20,
        #importance = "impurity_corrected"
        nfilter = num_top_genes,
        mtry = floor(0.2*nrow(norm_counts_data_holdout)),
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
  
  # If the final fit is NULL, return NULL
  if (is.null(final_fit)) {
    return(NULL)
  }
  
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
    auc_val = auc_val,
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
    roc_data_val = roc_perf_val,
    pred_val = pred_val,
    repeated_fit = repeated_fit,
    final_fit = final_fit,
    meta_ordered = ordered_filtered_counts$training$meta_ordered
  )
  
  
  # Build temp results list: append the result to results_list and save it
  #results_list_temp <<- append(results_list_temp, list(result))
  
  return(result)
}, future.packages = c("caret", "dplyr", "limma", "DESeq2", "sva", "doParallel", "foreach", "progress", "glmnet", "ROCR", "future.apply", "BiocParallel"), future.seed = TRUE)
end_time <- Sys.time()
print(paste("Total Operation Time: ",end_time - start_time))

### Prepare results for visualization
# Drop failed result from results
# Assuming results_list is your list
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
} else{
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
} else{
  coefficients_df <- NULL
  coefficients_df_features <- NULL
  feature_counts <- NULL

}

### CREATE ROC PLOTS
if(TRUE){
  all_training_inner_roc <- sapply(results_list, function(x) x$roc_data_tra_inner) # inner
  all_training_outter_roc <- sapply(results_list, function(x) x$roc_data_tra_outter) # outter
  all_validation_roc <- if (!is.null(results_list[[1]]$roc_data_val)) {
    sapply(results_list, function(x) x$roc_data_val)
  } else {
    NULL
  }
  
  source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
  training_inner_roc_plot <- plot_nested_roc_curves(all_training_inner_roc,  "")
  training_outter_roc_plot <- plot_nested_roc_curves(all_training_outter_roc, "")
  validation_roc_plot <- if (!is.null(all_validation_roc)) {
    plot_rocr_roc_curves(all_validation_roc, "Holdout ROC")
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

### CREATE PREDICTION PROBABILITY BOXPLOTS:
if(FALSE){
    
  source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
  prediction_df <- create_prediction_dataframe(results_list)
  
  
  source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
  generate_boxplots_ggplot2(prediction_df)
  
  ## Shapiro-Wilk Test for normality
  #shapiro_test_A <- shapiro.test(prediction_df$probability[prediction_df$condition == "PMpos"])
  #shapiro_test_B <- shapiro.test(prediction_df$probability[prediction_df$condition == "PMneg"])
  #print(shapiro_test_A)
  #print(shapiro_test_B)
}

### CREATE COEFFICIENT STABILITY PLOT
if(TRUE){
  source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
  ### Incorporate feature stability
  tryCatch({
    all_model_stability <- create_stability_plot(results_list, 0.5*length(part_rand_seed)*n_outer_folds_param,length(part_rand_seed)*n_outer_folds_param)
    # output first list element
    all_model_stability[[1]]
    stability_df_grouped <- all_model_stability[[2]]
  }, error = function(e) {
    message("An error occurred during stability plot creation: ", e$message)
    all_model_stability <- NULL
    stability_df_grouped <- NULL
  })
}

### GENERATE HEATMAPS
if(TRUE){
  source("~/5hmC-Pathway-Analysis/R Code/repeated_glmnet_funcs_v1.R")
  tryCatch({
    combined_heatmap <- generate_heatmaps(training_data = data_prepared$counts_data, 
                                        meta_training = data_prepared$meta, 
                                        annotated_features = stability_df_grouped)
  
    combined_heatmap
  }, error = function(e) {
    message("An error occurred during heatmap generation: ", e$message)
    combined_heatmap <- NULL
  })
}

### SAVE FILES
if(save_final_files){
  ### Save key results and files to output folder
  # Define variable to post-pend to output files that incorporates key parameters
  interaction_string <- paste(substr(interaction_params, 1, 4), collapse = "-")
  
  file_version <- paste0("NESTEDCV_v2_num_runs_", number_of_runs, "_nfolds_", nfolds_val, "_num_top_genes_", num_top_genes, "_seeds_", length(part_rand_seed),"_params",interaction_string)
  # Create output folder if it does not exist, based on file version and date
  output_folder_ver <- paste0(output_folder, file_version, "_", format(Sys.Date(), "%Y%m%d"), "_", format(Sys.time(), "%H%M%S"))# if it doesn't exist, create it
  if (!dir.exists(output_folder_ver)) {
    dir.create(output_folder_ver)
  }
  # Save key results and files to output folder
  write.csv(results_df, file = paste0(output_folder_ver, "/resultsdf_", file_version, ".csv"), row.names = FALSE)
  write.csv(coefficients_df, file = paste0(output_folder_ver, "/coefficients_", file_version, ".csv"), row.names = FALSE)
  saveRDS(median_result$repeated_fit, file = paste0(output_folder_ver, "/median_model_", file_version, ".rds"))
  saveRDS(training_inner_roc_plot, file = paste0(output_folder_ver, "/training_inner_roc_plot", file_version, ".rds"))
  saveRDS(training_outter_roc_plot, file = paste0(output_folder_ver, "/training_outter_roc_plot", file_version, ".rds"))
  saveRDS(validation_roc_plot, file = paste0(output_folder_ver, "/validation_roc_plot", file_version, ".rds"))
  htmlwidgets::saveWidget(training_inner_roc_plot, paste0(output_folder_ver, "/training_inner_roc_curves_", file_version, ".html"))
  htmlwidgets::saveWidget(training_outter_roc_plot, paste0(output_folder_ver, "/training_outter_roc_curves_", file_version, ".html"))
  htmlwidgets::saveWidget(combined_roc_plot, paste0(output_folder_ver, "/combined_roc_curves_", file_version, ".html"))
  
  # save all_model_stability as html
  htmlwidgets::saveWidget(all_model_stability[[1]], paste0(output_folder_ver, "/all_model_stability_plot_", file_version, ".html"))
  #save stability_df_grouped
  write.csv(stability_df_grouped, file = paste0(output_folder_ver, "/stability_df_grouped_", file_version, ".csv"), row.names = FALSE)
  
  
  ### Save all plots as SVG files
  # Save all plots as SVG files using kaleido
  plotly::save_image(training_inner_roc_plot, file = paste0(output_folder_ver, "/training_inner_roc_curves_", file_version, ".svg"))
  plotly::save_image(all_training_outter_roc, file = paste0(output_folder_ver, "/training_outter_roc_curves_", file_version, ".svg"))
  plotly::save_image(training_inner_roc_plot, file = paste0(output_folder_ver, "/training_inner_roc_curves_", file_version, ".jpg"))
  plotly::save_image(all_training_outter_roc, file = paste0(output_folder_ver, "/training_outter_roc_curves_", file_version, ".jpg"))
  # Save plots only if they are not NULL
  if (!is.null(validation_roc_plot)) {
    plotly::save_image(validation_roc_plot, file = paste0(output_folder_ver, "/validation_roc_curves_", file_version, ".svg"))
    plotly::save_image(validation_roc_plot, file = paste0(output_folder_ver, "/validation_roc_curves_", file_version, ".jpg"))
    htmlwidgets::saveWidget(validation_roc_plot, paste0(output_folder_ver, "/validation_roc_curves_", file_version, ".html"))
  }
  #plotly::save_image(box_plot, file = paste0(output_folder_ver, "/box_plot", file_version, ".jpg"))
  #plotly::save_image(bar_plot, file = paste0(output_folder_ver, "/bar_plot", file_version, ".jpg"))
  #plotly::save_image(box_plot, file = paste0(output_folder_ver, "/box_plot", file_version, ".svg"))
  #plotly::save_image(bar_plot, file = paste0(output_folder_ver, "/bar_plot", file_version, ".svg"))
  #htmlwidgets::saveWidget(box_plot, paste0(box_plot, "/combined_roc_curves_", file_version, ".html"))
  #htmlwidgets::saveWidget(bar_plot, paste0(bar_plot, "/combined_roc_curves_", file_version, ".html"))
  
  #plotly::save_image(bar_plot, file = paste0(output_folder_ver, "/bar_plot_", file_version, ".svg"))
  #plotly::save_image(box_plot, file = paste0(output_folder_ver, "/box_plot_", file_version, ".svg"))
  
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
    cat("Number of repeats for validation: ", rep_val, "\n")
    cat("Feature occurrence thresholds: ", paste(feat_occ_thresholds, collapse = ", "), "\n")
    cat("Alpha values: ", paste(alpha_values, collapse = ", "), "\n")
    cat("Standardize: ", standardize_bit, "\n")
    cat("Number of runs: ", number_of_runs, "\n")
    cat("Partition random seed: ", paste(part_rand_seed, collapse = ", "), "\n")
    cat("Number of inner folds: ", n_inner_folds_param, "\n")
    cat("Number of outer folds: ", n_outer_folds_param, "\n")
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
  })
  
  # Write the captured output to the file
  writeLines(config_content, con = paste0(output_folder_ver, "/config.txt"))
}

### SAve heatmaps and probability boxplots
if(save_final_files){
  svg(paste0(output_folder_ver, "/combined_heatmap_", file_version, ".svg"))
  draw(combined_heatmap)
  dev.off()
  
  jpeg(paste0(output_folder_ver, "/combined_heatmap_", file_version, ".jpg"))
  draw(combined_heatmap)
  dev.off()
  ggsave(paste0(output_folder_ver, "/probability_boxplots_", file_version, ".svg"), width = 10, height = 10, units = "in")
  ggsave(paste0(output_folder_ver, "/probability_boxplots_", file_version, ".jpg"), width = 10, height = 10, units = "in")
}



### TO DO:
## Heatmap of most important genes
## Try adding validation set, confirm performance holds

## Look up enhancers for distal intergenic regions in annotations

## Fix subtitles for ROC plots, 
# Update labeling, colors
## Plot prediction score for outer folds (training vs validation) per different groups (PMpos, PMneg, Primary_Present, Primary_absent, CRC, HGA, Healthy)
## Plot heatmap for resulting genes, on median model training vs validation set
# Plot 5hmc Levels across each gene (PM+ vs PM-) using metagene2 or https://showteeth.github.io/ggcoverage/


### Plotting individual genes:
#  Install wiggle plotr
library("wiggleplotr")
library("dplyr")
library("GenomicRanges")
library("GenomicFeatures")
require("org.Hs.eg.db")
require("TxDb.Hsapiens.UCSC.hg19.knownGene")
# install metagene2 with bioconductor
library("metagene2")

### FUNCTION
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


### PLOTTNG METAGENE PLOTS FOR SPECIFIC REGIONS
## Construct regions
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
regions_pos <- regions[regions$coefs >= 0 & regions$frequency_percentage >= 50, ]
regions_neg <- regions[regions$coefs < 0 & regions$frequency_percentage >= 50, ]


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
                              PMpos = ifelse(meta_tracks$condition == "PMpos", 1, 0),
                              PMneg = ifelse(meta_tracks$condition == "PMneg", 1, 0)
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
mg$produce_metagene(title = "Demo metagene plot",design=metagene_design, normalization = "RPM",alpha=1,facet_by=~region)

metagene2_heatmap(mg)

### Plotting metagene plot for all gene bodies, and enhancers:
# Import health control meta

meta_healthy <- read.csv("~/5hmC-Pathway-Analysis/Raw Input/Working Inputs/CRC_HGA_HEALTHY_05042024.csv", row.names = 1)

all_gene_bodies = GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
# Take random 1000 sample
set.seed(123)  # Set seed for reproducibility
sampled_gene_bodies <- all_gene_bodies[sample(1:length(all_gene_bodies), 1000)]

all_TSS = GenomicFeatures::promoters(TxDb.Hsapiens.UCSC.hg19.knownGene)

meta_tracks <- meta_healthy
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
                              CRC_HGA = ifelse(meta_tracks$condition == "CRC_HGA", 1, 0),
                              HEALTHY = ifelse(meta_tracks$condition == "HEALTHY", 1, 0)
)
# set rownames to sample_name
rownames(metagene_design) <- metagene_design$sample_name
# drop sample_name column
metagene_design <- metagene_design[, -1]
head(metagene_design)

mg_bam_files <- meta_tracks$bam_fn

mg <- metagene2::metagene2$new(regions = sampled_gene_bodies, 
                               bam_files = mg_bam_files, 
                               assay='chipseq',
                               paired_end=TRUE,
                               cores=80,
                               padding_size = 1000,
                               force_seqlevels = TRUE
                               #extend_reads=180
)

# We can then plot coverage over those regions across all bam files.
mg$produce_metagene(title = "Demo metagene plot",design=metagene_design, facet_by = NULL, normalization = "RPM",alpha=1, bin_count=1000) #facet_by=~region, 


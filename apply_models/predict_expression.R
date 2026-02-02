#!/usr/bin/Rscript

#' Predict Gene Expression Using BMI-Dynamic Models
#' 
#' This script predicts gene expression in biobank cohorts using pre-trained
#' BMI-dynamic models. It accepts genotype data, BMI values, and outputs
#' predicted expression (GREX) for genes on a specified chromosome and part.
#' 
#' The part system splits genes into chunks for memory management in large cohorts.
#'
#' @author [Your Name]
#' @date [Date]

# ============================================================================
# LOAD REQUIRED LIBRARIES
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(optparse)
})

# ============================================================================
# PARSE COMMAND LINE ARGUMENTS
# ============================================================================

option_list <- list(
  make_option(c("--tissue"), type="character", default=NULL,
              help="Tissue name (e.g., Adipose-Subcutaneous, Brain-Cortex)",
              metavar="CHARACTER"),
  
  make_option(c("--chr"), type="integer", default=NULL,
              help="Chromosome number (1-22)",
              metavar="INTEGER"),
  
  make_option(c("--part"), type="integer", default=NULL,
              help="Part number (1-10, splits genes for memory management)",
              metavar="INTEGER"),
  
  make_option(c("--genotype-file"), type="character", default=NULL,
              help="Path to genotype file (.raw format with GTEx column names)",
              metavar="PATH"),
  
  make_option(c("--bmi-file"), type="character", default=NULL,
              help="Path to BMI file (columns: SUBJID, BMI)",
              metavar="PATH"),
  
  make_option(c("--output-dir"), type="character", default=NULL,
              help="Output directory for results",
              metavar="PATH"),
  
  make_option(c("--model-dir"), type="character", default="models",
              help="Directory containing model files [default: %default]",
              metavar="PATH"),
  
  make_option(c("--n-parts"), type="integer", default=10,
              help="Total number of parts to split genes into [default: %default]",
              metavar="INTEGER")
)

opt_parser <- OptionParser(
  option_list=option_list,
  description="Predict gene expression using BMI-dynamic models",
  epilogue="Example:\n  Rscript predict_expression.R --tissue Adipose-Subcutaneous --chr 1 --part 1 --genotype-file biobank_chr1.raw.gz --bmi-file biobank_bmi.txt --output-dir results/"
)

opt <- parse_args(opt_parser)

# ============================================================================
# VALIDATE ARGUMENTS
# ============================================================================

validate_arguments <- function(opt) {
  errors <- c()
  
  # Check required arguments
  if (is.null(opt$tissue)) errors <- c(errors, "  --tissue is required")
  if (is.null(opt$chr)) errors <- c(errors, "  --chr is required")
  if (is.null(opt$part)) errors <- c(errors, "  --part is required")
  if (is.null(opt$`genotype-file`)) errors <- c(errors, "  --genotype-file is required")
  if (is.null(opt$`bmi-file`)) errors <- c(errors, "  --bmi-file is required")
  if (is.null(opt$`output-dir`)) errors <- c(errors, "  --output-dir is required")
  
  if (length(errors) > 0) {
    cat("Error: Missing required arguments:\n")
    cat(paste(errors, collapse="\n"), "\n\n")
    print_help(opt_parser)
    quit(status=1)
  }
  
  # Validate chromosome
  if (opt$chr < 1 || opt$chr > 22) {
    cat("Error: Chromosome must be between 1 and 22\n")
    quit(status=1)
  }
  
  # Validate part
  if (opt$part < 1 || opt$part > opt$`n-parts`) {
    cat(sprintf("Error: Part must be between 1 and %d\n", opt$`n-parts`))
    quit(status=1)
  }
  
  # Check input files exist
  if (!file.exists(opt$`genotype-file`)) {
    cat(sprintf("Error: Genotype file not found: %s\n", opt$`genotype-file`))
    quit(status=1)
  }
  
  if (!file.exists(opt$`bmi-file`)) {
    cat(sprintf("Error: BMI file not found: %s\n", opt$`bmi-file`))
    quit(status=1)
  }
  
  # Check model directory exists
  if (!dir.exists(opt$`model-dir`)) {
    cat(sprintf("Error: Model directory not found: %s\n", opt$`model-dir`))
    quit(status=1)
  }
  
  # Check model file exists
  model_file <- sprintf("%s/tissues/%s_r2_0.01_p_0.05_betas.txt.gz", 
                       opt$`model-dir`, opt$tissue)
  if (!file.exists(model_file)) {
    cat(sprintf("Error: Model file not found: %s\n", model_file))
    cat(sprintf("  Available tissues can be found in: %s/tissues/\n", opt$`model-dir`))
    quit(status=1)
  }
  
  # Create output directory if needed
  if (!dir.exists(opt$`output-dir`)) {
    dir.create(opt$`output-dir`, recursive=TRUE, showWarnings=FALSE)
    cat(sprintf("Created output directory: %s\n", opt$`output-dir`))
  }
}

validate_arguments(opt)

# ============================================================================
# SETUP
# ============================================================================

tissue_name <- opt$tissue
chr_here <- opt$chr
part_here <- opt$part
n_parts <- opt$`n-parts`
genotype_file <- opt$`genotype-file`
bmi_file <- opt$`bmi-file`
output_dir <- opt$`output-dir`
model_dir <- opt$`model-dir`

cat("\n")
cat("================================================================================\n")
cat(" Predicting Expression Using BMI-Dynamic Models\n")
cat("================================================================================\n")
cat(sprintf("  Tissue:           %s\n", tissue_name))
cat(sprintf("  Chromosome:       %d\n", chr_here))
cat(sprintf("  Part:             %d of %d\n", part_here, n_parts))
cat(sprintf("  Genotype file:    %s\n", genotype_file))
cat(sprintf("  BMI file:         %s\n", bmi_file))
cat(sprintf("  Model directory:  %s\n", model_dir))
cat(sprintf("  Output directory: %s\n", output_dir))
cat("================================================================================\n\n")

# ============================================================================
# LOAD MODEL COEFFICIENTS
# ============================================================================

cat("Loading model coefficients...\n")
model_file <- sprintf("%s/tissues/%s_r2_0.01_p_0.05_betas.txt.gz", model_dir, tissue_name)

all_models <- tryCatch({
  as.data.frame(fread(model_file, header=TRUE, sep="\t", stringsAsFactors=FALSE))
}, error = function(e) {
  cat(sprintf("Error reading model file: %s\n", e$message))
  quit(status=1)
})

# Filter to this chromosome
models_here <- all_models %>% filter(Chr == chr_here)

if (nrow(models_here) == 0) {
  cat(sprintf("Warning: No models found for chromosome %d in tissue %s\n", chr_here, tissue_name))
  cat("  This may be expected if no genes on this chromosome passed QC.\n")
  quit(status=0)
}

cat(sprintf("  Loaded models for chromosome %d\n", chr_here))

# ============================================================================
# SPLIT GENES INTO PARTS
# ============================================================================

# Get all genes for this chromosome
all_genes_chr <- unique(models_here$Gene)
n_genes_total <- length(all_genes_chr)

# Split genes into equal parts
gene_parts <- split(all_genes_chr, cut(seq_along(all_genes_chr), n_parts, labels = FALSE))
gene_part_here <- gene_parts[[part_here]]

if (is.null(gene_part_here) || length(gene_part_here) == 0) {
  cat(sprintf("Warning: Part %d has no genes for chromosome %d\n", part_here, chr_here))
  quit(status=0)
}

cat(sprintf("  Total genes on chromosome %d: %d\n", chr_here, n_genes_total))
cat(sprintf("  Genes in part %d: %d\n", part_here, length(gene_part_here)))

# Filter models to genes in this part
models_here <- models_here %>% filter(Gene %in% gene_part_here)

# ============================================================================
# LOAD BMI DATA
# ============================================================================

cat("Loading BMI data...\n")
BMI_df <- tryCatch({
  as.data.frame(fread(bmi_file, header=TRUE, sep="\t", stringsAsFactors=FALSE))
}, error = function(e) {
  cat(sprintf("Error reading BMI file: %s\n", e$message))
  quit(status=1)
})

# Check required columns
if (!("SUBJID" %in% colnames(BMI_df)) || !("BMI" %in% colnames(BMI_df))) {
  cat("Error: BMI file must contain columns 'SUBJID' and 'BMI'\n")
  cat(sprintf("  Found columns: %s\n", paste(colnames(BMI_df), collapse=", ")))
  quit(status=1)
}

# Keep only required columns
BMI_df <- BMI_df %>% select(SUBJID, BMI)

cat(sprintf("  Loaded BMI data for %d individuals\n", nrow(BMI_df)))

# ============================================================================
# LOAD GENOTYPE DATA
# ============================================================================

cat("Loading genotype data...\n")
geno_raw <- tryCatch({
  as.data.frame(fread(genotype_file, header=TRUE, sep="\t", stringsAsFactors=FALSE))
}, error = function(e) {
  cat(sprintf("Error reading genotype file: %s\n", e$message))
  quit(status=1)
})

# Check for IID column
if (!("IID" %in% colnames(geno_raw))) {
  cat("Error: Genotype file must contain column 'IID'\n")
  cat(sprintf("  Found columns: %s\n", paste(colnames(geno_raw)[1:min(5, ncol(geno_raw))], collapse=", ")))
  quit(status=1)
}

# Set IID as rownames and rename to SUBJID for consistency
geno_raw <- geno_raw %>% 
  rename(SUBJID = IID) %>%
  column_to_rownames("SUBJID")

cat(sprintf("  Loaded genotype data: %d individuals, %d SNPs\n", 
            nrow(geno_raw), ncol(geno_raw)))

# Check overlap with BMI data
n_overlap <- length(intersect(rownames(geno_raw), BMI_df$SUBJID))
cat(sprintf("  Overlap with BMI data: %d individuals\n", n_overlap))

if (n_overlap == 0) {
  cat("Error: No overlapping individuals between genotype and BMI files\n")
  cat("  Check that SUBJID in BMI file matches IID in genotype file\n")
  quit(status=1)
}

# ============================================================================
# PREDICT EXPRESSION FOR EACH GENE IN THIS PART
# ============================================================================

cat(sprintf("\nPredicting expression for %d genes in part %d...\n", length(gene_part_here), part_here))

n_genes <- length(gene_part_here)

# Initialize results
all_genes <- data.frame(SUBJID = rownames(geno_raw))
impute_stats <- NULL

# Progress tracking
progress_interval <- max(1, floor(n_genes / 10))

for (i in seq_along(gene_part_here)) {
  gene_here <- gene_part_here[i]
  
  # Show progress
  if (i %% progress_interval == 0) {
    cat(sprintf("  Progress: %d/%d genes (%.1f%%)\n", i, n_genes, 100*i/n_genes))
  }
  
  # Get predictors for this gene
  predictors_here <- models_here %>% filter(Gene == gene_here)
  
  # Get SNPs (exclude BMI from predictor list)
  SNPs_here <- setdiff(predictors_here$name, 'BMI')
  SNPs_overlap <- intersect(SNPs_here, colnames(geno_raw))
  
  # Prepare genotype data for this gene
  if (length(SNPs_overlap) > 0) {
    # Subset genotype matrix
    if (length(SNPs_overlap) == 1) {
      geno_here <- data.frame(geno_raw[, SNPs_overlap])
      colnames(geno_here) <- SNPs_overlap
    } else {
      geno_here <- geno_raw[, SNPs_overlap, drop=FALSE]
    }
    
    # Merge with BMI
    geno_merge <- geno_here %>% rownames_to_column('SUBJID')
    pred_df <- inner_join(BMI_df, geno_merge, by='SUBJID')
    
  } else {
    # No SNPs available - only BMI
    pred_df <- BMI_df
  }
  
  # ====================================================================
  # Create interaction terms (SNP*BMI)
  # ====================================================================
  
  all_predictors <- unique(predictors_here$Predictor)
  interactions <- all_predictors[grepl('\\*', all_predictors)]
  
  if (length(interactions) > 0) {
    # Extract SNPs that are in interactions
    interacting_snps <- gsub('\\*BMI', '', interactions)
    interacting_snps_here <- intersect(interacting_snps, colnames(pred_df))
    
    # Create interaction terms
    if (length(interacting_snps_here) > 0) {
      for (snp_int in interacting_snps_here) {
        int_name <- paste0(snp_int, '*BMI')
        # Multiply SNP dosage by BMI
        pred_df[[int_name]] <- pred_df[[snp_int]] * pred_df[['BMI']]
      }
    }
  }
  
  # ====================================================================
  # Calculate BMI-aware predicted expression (BGREX)
  # ====================================================================
  
  # Set SUBJID as rownames for matrix operations
  pred_matrix <- pred_df %>% 
    column_to_rownames('SUBJID') %>%
    as.matrix()
  
  # Transpose to get predictors as rows
  predictor_df_t <- t(pred_matrix)
  
  # Get available predictors
  pred_used <- intersect(rownames(predictor_df_t), predictors_here$Predictor)
  
  if (length(pred_used) > 0) {
    # Subset to available predictors
    predictor_subset <- predictor_df_t[pred_used, , drop=FALSE]
    class(predictor_subset) <- "numeric"
    
    # Calculate weighted sum (BGREX)
    if (length(pred_used) == 1) {
      # Single predictor
      Beta_here <- predictors_here[predictors_here$Predictor == pred_used, 'Beta']
      BGREX_vec <- predictor_subset * Beta_here
      save_thisgene <- data.frame(
        SUBJID = colnames(predictor_df_t),
        value = as.numeric(BGREX_vec)
      )
      
    } else {
      # Multiple predictors
      BGREX_mat <- matrix(0, nrow=nrow(predictor_subset), ncol=ncol(predictor_subset))
      
      for (row_idx in 1:nrow(predictor_subset)) {
        pred_name <- rownames(predictor_subset)[row_idx]
        Beta_here <- predictors_here[predictors_here$Predictor == pred_name, 'Beta']
        BGREX_mat[row_idx, ] <- predictor_subset[row_idx, ] * Beta_here
      }
      
      # Sum across predictors
      BGREX_sum <- colSums(BGREX_mat)
      save_thisgene <- data.frame(
        SUBJID = names(BGREX_sum),
        value = as.numeric(BGREX_sum)
      )
    }
    
    # Record imputation statistics
    impute_here <- data.frame(
      gene = gene_here,
      N_pred_all = nrow(predictors_here),
      N_pred_used = length(pred_used),
      pred_used = paste(pred_used, collapse=",")
    )
    
  } else {
    # No predictors available
    save_thisgene <- data.frame(
      SUBJID = colnames(predictor_df_t),
      value = NA
    )
    
    impute_here <- data.frame(
      gene = gene_here,
      N_pred_all = nrow(predictors_here),
      N_pred_used = 0,
      pred_used = NA
    )
  }
  
  # Rename value column to gene name
  colnames(save_thisgene)[2] <- gene_here
  
  # Add to results
  all_genes <- left_join(all_genes, save_thisgene, by='SUBJID')
  impute_stats <- rbind(impute_stats, impute_here)
}

cat(sprintf("  Completed: %d/%d genes\n", n_genes, n_genes))

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("\nSaving results...\n")

# Output filenames include part number
bgrex_file <- sprintf("%s/%s_%d_%d_BGREX.txt.gz", output_dir, tissue_name, chr_here, part_here)
stats_file <- sprintf("%s/%s_%d_%d_imputestats.txt.gz", output_dir, tissue_name, chr_here, part_here)

# Save predicted expression
write.table(all_genes, 
            file = gzfile(bgrex_file),
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

cat(sprintf("  Saved predicted expression: %s\n", bgrex_file))

# Save imputation statistics
write.table(impute_stats,
            file = gzfile(stats_file),
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

cat(sprintf("  Saved imputation statistics: %s\n", stats_file))

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\nSummary:\n")
cat("================================================================================\n")
cat(sprintf("  Chromosome:             %d\n", chr_here))
cat(sprintf("  Part:                   %d of %d\n", part_here, n_parts))
cat(sprintf("  Genes in this part:     %d\n", n_genes))
cat(sprintf("  Genes predicted:        %d (%.1f%%)\n", 
            sum(impute_stats$N_pred_used > 0),
            100 * sum(impute_stats$N_pred_used > 0) / n_genes))
cat(sprintf("  Genes not predicted:    %d (%.1f%%)\n",
            sum(impute_stats$N_pred_used == 0),
            100 * sum(impute_stats$N_pred_used == 0) / n_genes))
cat(sprintf("  Mean predictors used:   %.1f / %.1f (%.1f%%)\n",
            mean(impute_stats$N_pred_used),
            mean(impute_stats$N_pred_all),
            100 * mean(impute_stats$N_pred_used / impute_stats$N_pred_all)))
cat("================================================================================\n")
cat("\nDone!\n\n")
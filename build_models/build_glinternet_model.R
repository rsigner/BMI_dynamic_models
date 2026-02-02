#!/usr/bin/Rscript

#' Build BMI-Dynamic Gene Expression Prediction Models
#' 
#' This script trains glinternet models to predict gene expression from genotype
#' data, incorporating BMI as a dynamic variable that modulates genetic effects.
#' 
#' The script processes genes in parts for memory management and uses 4-fold
#' cross-validation with one held-out fold for validation.
#'
#' @author [Your Name]
#' @date 2026

# ============================================================================
# LOAD REQUIRED LIBRARIES
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(glinternet)
  library(optparse)
})

# ============================================================================
# PARSE COMMAND LINE ARGUMENTS
# ============================================================================

option_list <- list(
  make_option(c("--tissue"), type="character", default=NULL,
              help="Tissue name (e.g., Brain-Cortex, Adipose-Subcutaneous)",
              metavar="CHARACTER"),
  
  make_option(c("--chr"), type="integer", default=NULL,
              help="Chromosome number (1-22)",
              metavar="INTEGER"),
  
  make_option(c("--part"), type="integer", default=NULL,
              help="Part number (1-30)",
              metavar="INTEGER"),
  
  make_option(c("--gene-list-file"), type="character", default=NULL,
              help="Path to file listing genes to process (one Ensembl ID per line)",
              metavar="PATH"),
  
  make_option(c("--fold-file"), type="character", default=NULL,
              help="Path to CV fold assignments (from create_folds.py)",
              metavar="PATH"),
  
  make_option(c("--phenotype-file"), type="character", default=NULL,
              help="Path to phenotype file (columns: SUBJID, BMI)",
              metavar="PATH"),
  
  make_option(c("--expression-file"), type="character", default=NULL,
              help="Path to expression residuals file (rows: samples, cols: genes)",
              metavar="PATH"),
  
  make_option(c("--genotype-file"), type="character", default=NULL,
              help="Path to genotype dosage file for this chromosome",
              metavar="PATH"),
  
  make_option(c("--feature-file"), type="character", default=NULL,
              help="Path to gene feature/annotation file",
              metavar="PATH"),
  
  make_option(c("--output-dir"), type="character", default=NULL,
              help="Output directory for model files",
              metavar="PATH"),
  
  make_option(c("--n-parts"), type="integer", default=30,
              help="Total number of parts to split genes into [default: %default]",
              metavar="INTEGER"),
  
  make_option(c("--cis-window"), type="integer", default=1000000,
              help="Cis-regulatory window size in bp [default: %default = 1Mb]",
              metavar="INTEGER"),
  
  make_option(c("--rsq-threshold"), type="numeric", default=0.01,
              help="Minimum R² threshold for saving models [default: %default]",
              metavar="NUMERIC"),
  
  make_option(c("--n-folds-cv"), type="integer", default=3,
              help="Number of internal CV folds for glinternet [default: %default]",
              metavar="INTEGER"),
  
  make_option(c("--random-seed"), type="integer", default=1702,
              help="Random seed for reproducibility [default: %default]",
              metavar="INTEGER")
)

opt_parser <- OptionParser(
  option_list=option_list,
  description="Build BMI-dynamic gene expression prediction models using glinternet",
  epilogue="Example:\n  Rscript build_glinternet_model.R --tissue Brain-Cortex --chr 1 --part 1 --gene-list-file genes_chr1.txt --fold-file folds.txt.gz --phenotype-file phenotypes.txt.gz --expression-file expression.txt.gz --genotype-file genotypes_chr1.txt.gz --feature-file features.txt.gz --output-dir output/"
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
  if (is.null(opt$`gene-list-file`)) errors <- c(errors, "  --gene-list-file is required")
  if (is.null(opt$`fold-file`)) errors <- c(errors, "  --fold-file is required")
  if (is.null(opt$`phenotype-file`)) errors <- c(errors, "  --phenotype-file is required")
  if (is.null(opt$`expression-file`)) errors <- c(errors, "  --expression-file is required")
  if (is.null(opt$`genotype-file`)) errors <- c(errors, "  --genotype-file is required")
  if (is.null(opt$`feature-file`)) errors <- c(errors, "  --feature-file is required")
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
  
  # Check files exist
  required_files <- c(
    opt$`gene-list-file`,
    opt$`fold-file`,
    opt$`phenotype-file`,
    opt$`expression-file`,
    opt$`genotype-file`,
    opt$`feature-file`
  )
  
  for (file in required_files) {
    if (!file.exists(file)) {
      cat(sprintf("Error: File not found: %s\n", file))
      quit(status=1)
    }
  }
  
  # Create output directory
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
gene_list_file <- opt$`gene-list-file`
fold_file <- opt$`fold-file`
phenotype_file <- opt$`phenotype-file`
expression_file <- opt$`expression-file`
genotype_file <- opt$`genotype-file`
feature_file <- opt$`feature-file`
output_dir <- opt$`output-dir`
n_parts <- opt$`n-parts`
cis_window <- opt$`cis-window`
rsq_threshold <- opt$`rsq-threshold`
n_folds_cv <- opt$`n-folds-cv`
random_seed <- opt$`random-seed`

set.seed(random_seed)

cat("\n")
cat("================================================================================\n")
cat(" Building BMI-Dynamic Gene Expression Prediction Models\n")
cat("================================================================================\n")
cat(sprintf("  Tissue:          %s\n", tissue_name))
cat(sprintf("  Chromosome:      %d\n", chr_here))
cat(sprintf("  Part:            %d of %d\n", part_here, n_parts))
cat(sprintf("  Gene list:       %s\n", gene_list_file))
cat(sprintf("  CV folds:        %s\n", fold_file))
cat(sprintf("  Phenotype:       %s\n", phenotype_file))
cat(sprintf("  Expression:      %s\n", expression_file))
cat(sprintf("  Genotype:        %s\n", genotype_file))
cat(sprintf("  Features:        %s\n", feature_file))
cat(sprintf("  Output dir:      %s\n", output_dir))
cat(sprintf("  Cis window:      ±%d bp (%.1f Mb)\n", cis_window, cis_window/1e6))
cat(sprintf("  R² threshold:    %.2f\n", rsq_threshold))
cat(sprintf("  Internal CV:     %d folds\n", n_folds_cv))
cat(sprintf("  Random seed:     %d\n", random_seed))
cat("================================================================================\n\n")

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading data...\n")

# Load gene features
cat("  Loading gene features...\n")
featuredt <- as.data.frame(fread(feature_file, header=T, sep="\t", stringsAsFactors=F))
cat(sprintf("    Loaded: %d genes\n", nrow(featuredt)))

# Load phenotype data
cat("  Loading phenotype data...\n")
phenodt_raw <- as.data.frame(fread(phenotype_file, header=T, sep="\t", stringsAsFactors=F))
pheno_here <- phenodt_raw[, c('SUBJID', 'BMI')]
cat(sprintf("    Loaded: %d individuals\n", nrow(pheno_here)))

# Load expression residuals
cat("  Loading expression residuals...\n")
residualsdf_raw <- as.data.frame(fread(expression_file, header=T, sep="\t", stringsAsFactors=F))
residuals <- residualsdf_raw %>% column_to_rownames("ID")
response_here <- residuals[pheno_here$SUBJID, ]
cat(sprintf("    Loaded: %d samples, %d genes\n", nrow(response_here), ncol(response_here)))

# Load genotype data
cat("  Loading genotype data...\n")
genotypedt <- as.data.frame(fread(genotype_file, header=T, sep="\t", stringsAsFactors=F))
genotypedt_mat <- genotypedt %>% column_to_rownames('ID')

# Remove SNPs with any NAs
snp_NA <- rownames(genotypedt_mat[rowSums(is.na(genotypedt_mat)) > 0, ])
genotypedt_mat_noNA <- genotypedt_mat[!(rownames(genotypedt_mat) %in% snp_NA), ]
cat(sprintf("    Loaded: %d samples, %d SNPs\n", nrow(genotypedt_mat_noNA), ncol(genotypedt_mat_noNA)))

if (length(snp_NA) > 0) {
  cat(sprintf("    Removed %d SNPs with missing values\n", length(snp_NA)))
}

# Load CV fold assignments
cat("  Loading CV fold assignments...\n")
fold_dat <- as.data.frame(fread(fold_file, header=T, sep="\t", stringsAsFactors=F))
cat(sprintf("    Loaded: %d individuals, %d folds\n", nrow(fold_dat), max(fold_dat$Fold_in_test)))

# Subset to training folds (hold out fold 4 for validation)
fold_CV <- fold_dat %>% filter(Fold_in_test != 4)
cat(sprintf("    Training folds (1-3): %d individuals\n", nrow(fold_CV)))
cat(sprintf("    Validation fold (4):  %d individuals\n", sum(fold_dat$Fold_in_test == 4)))

# Subset data to training samples
response_final <- response_here[fold_CV$SUBJID, ]
geno_final <- genotypedt_mat_noNA[fold_CV$SUBJID, ]

# Prepare phenotype for merging
samp2subj <- pheno_here

# ============================================================================
# PREPARE SNP POSITION DATA
# ============================================================================

cat("\nPreparing SNP position data...\n")

# Extract chromosome and position from SNP names (GTEx format: chr1_12345_A_G_b38)
snp_names <- colnames(geno_final)
snp_info <- data.frame(SNP = snp_names)

# Parse SNP names to get positions
snp_info$Chr <- sapply(strsplit(snp_names, "_"), function(x) gsub("chr", "", x[1]))
snp_info$pos <- sapply(strsplit(snp_names, "_"), function(x) x[2])

snp_clean_df <- snp_info %>% filter(Chr == chr_here)
cat(sprintf("  SNPs on chromosome %d: %d\n", chr_here, nrow(snp_clean_df)))

# ============================================================================
# SPLIT GENES INTO PARTS
# ============================================================================

cat("\nDetermining genes to process...\n")

# Read gene list
genes_all <- read_lines(gene_list_file)
cat(sprintf("  Total genes in list: %d\n", length(genes_all)))

# Split into parts
gene_parts <- split(genes_all, cut(seq_along(genes_all), n_parts, labels = FALSE))
gene_part_here <- gene_parts[[part_here]]

if (is.null(gene_part_here) || length(gene_part_here) == 0) {
  cat(sprintf("Warning: Part %d has no genes\n", part_here))
  quit(status=0)
}

cat(sprintf("  Genes in part %d: %d\n", part_here, length(gene_part_here)))

# ============================================================================
# PROCESS GENES
# ============================================================================

cat("\n")
cat("================================================================================\n")
cat(sprintf(" Processing %d genes in part %d\n", length(gene_part_here), part_here))
cat("================================================================================\n\n")

n_genes <- length(gene_part_here)
n_success <- 0
n_failed <- 0
progress_interval <- max(1, floor(n_genes / 10))

for (i in seq_along(gene_part_here)) {
  gene <- gene_part_here[i]
  
  # Show progress
  if (i %% progress_interval == 0) {
    cat(sprintf("Progress: %d/%d genes (%.1f%%)\n", i, n_genes, 100*i/n_genes))
  }
  
  cat(sprintf("  [%d/%d] Processing %s...\n", i, n_genes, gene))
  
  # Check if gene exists in expression data
  if (!(gene %in% colnames(response_final))) {
    cat(sprintf("    Warning: Gene %s not found in expression data, skipping\n", gene))
    n_failed <- n_failed + 1
    next
  }
  
  # Get expression for this gene
  exp_dat <- response_final %>% select(all_of(gene))
  exp_dat <- exp_dat[fold_CV$SUBJID, , drop=FALSE]
  
  # Get gene features
  features_here <- featuredt %>% filter(Name == gene)
  
  if (nrow(features_here) == 0) {
    cat(sprintf("    Warning: Gene %s not found in feature file, skipping\n", gene))
    n_failed <- n_failed + 1
    next
  }
  
  # Define cis-window
  start_site <- as.numeric(features_here$start)
  cis_upper <- start_site + cis_window
  cis_lower <- max(0, start_site - cis_window)
  
  cat(sprintf("    Cis window: %d - %d (±%.1f Mb)\n", cis_lower, cis_upper, cis_window/1e6))
  
  # Identify SNPs in cis-window
  snp_cols <- snp_clean_df %>%
    filter(as.numeric(pos) >= cis_lower) %>%
    filter(as.numeric(pos) <= cis_upper) %>%
    .$SNP
  
  snp_shared <- intersect(snp_cols, colnames(geno_final))
  
  if (length(snp_shared) == 0) {
    cat(sprintf("    No SNPs in cis-window, skipping\n"))
    
    # Save failed model summary
    model_summary <- data.frame(
      Gene = gene,
      Tissue = tissue_name,
      Rsq = NA,
      p_Rsq = NA,
      CV_Lambda = NA,
      CV_MStEr = NA,
      CV_MSE = NA,
      R2_rep = NA,
      p_rep = NA,
      n_predictors = 0
    )
    
    out_summary <- paste0(output_dir, "/", tissue_name, "_", gene, "_model_summary.txt.gz")
    write.table(model_summary, file = gzfile(out_summary), sep = "\t", 
                col.names=T, row.names = F, quote=F)
    
    n_failed <- n_failed + 1
    next
  }
  
  cat(sprintf("    SNPs in window: %d\n", length(snp_shared)))
  
  # Prepare genotype data
  geno_here <- geno_final %>%
    as.data.frame %>%
    select(all_of(snp_shared)) %>%
    rownames_to_column('SUBJID')
  
  # Merge with BMI
  predictors <- left_join(geno_here, samp2subj[, c('BMI', 'SUBJID')], by='SUBJID')
  predictors_final <- predictors[, c('BMI', snp_shared, 'SUBJID')]
  predictors_final <- predictors_final %>% column_to_rownames('SUBJID')
  
  # Run glinternet
  cat("    Running glinternet...\n")
  
  # All predictors are continuous (for dosage data)
  num_levels <- c(1, rep(1, length(snp_shared)))
  
  tryCatch({
    fit <- glinternet.cv(
      X = predictors_final,
      Y = exp_dat,
      nFolds = n_folds_cv,
      interactionCandidates = c(1),  # BMI interactions
      numLevels = num_levels
    )
    
    # Get best lambda
    best_lambda <- fit$lambdaHat
    idx_best <- which(fit$lambda == best_lambda)
    min_cvErr <- fit$cvErr[idx_best]
    min_cvErrStd <- fit$cvErrStd[idx_best]
    
    # Calculate training R²
    y_actual <- exp_dat
    y_predicted <- fit$fitted
    
    if (var(y_predicted) == 0) {
      R_squared <- 0
      pval_lm <- NA
    } else {
      res <- summary(lm(y_actual ~ y_predicted))
      R_squared <- res$r.squared
      pval_lm <- res$coefficients[2, 4]
    }
    
    cat(sprintf("    Training R² = %.4f, P = %.2e\n", R_squared, pval_lm))
    
    # Check if model meets threshold
    if (R_squared > rsq_threshold) {
      cat("    Model meets R² threshold, validating on held-out fold...\n")
      
      # Validate on held-out fold 4
      fold_oos <- fold_dat %>% filter(Fold_in_test == 4)
      geno_rep <- genotypedt_mat_noNA[fold_oos$SUBJID, ]
      geno_rep <- geno_rep %>%
        as.data.frame %>%
        select(all_of(snp_shared)) %>%
        rownames_to_column('SUBJID')
      
      # Merge with BMI
      predictors_rep <- left_join(geno_rep, samp2subj[, c('BMI', 'SUBJID')], by='SUBJID')
      predictors_rep_final <- predictors_rep[, c('BMI', snp_shared, 'SUBJID')]
      predictors_rep_final <- predictors_rep_final %>% column_to_rownames('SUBJID')
      
      # Predict on validation fold
      predicted_11th <- predict(fit, predictors_rep_final)
      
      # Get observed values
      response_rep <- response_here[fold_oos$SUBJID, gene]
      
      # Calculate validation R²
      res_rep <- summary(lm(response_rep ~ predicted_11th))
      rsq_rep <- res_rep$r.squared
      
      if (rsq_rep > 0) {
        pval_rep <- res_rep$coefficients[2, 4]
      } else {
        pval_rep <- NA
      }
      
      cat(sprintf("    Validation R² = %.4f, P = %.2e\n", rsq_rep, pval_rep))
      
      # Extract coefficients
      coef_all <- coef(fit, lambdaIndex=idx_best)
      
      # Main effects
      idx_main <- coef_all$mainEffects$cont
      snp_names_main <- colnames(predictors_final)[idx_main]
      coef_main <- coef_all$mainEffectsCoef$cont
      names(coef_main) <- snp_names_main
      main_df <- enframe(coef_main) %>% unnest(cols=2)
      main_df$predictor_type <- 'Main Effect'
      main_df$Predictor <- main_df$name
      
      # Interactions
      idx_int <- coef_all$interactions$contcont[, 2]
      snp_names_int <- colnames(predictors_final)[idx_int]
      coef_int <- coef_all$interactionsCoef$contcont
      names(coef_int) <- snp_names_int
      int_df <- enframe(coef_int) %>% unnest(cols=2)
      int_df$predictor_type <- 'Interaction Effect'
      int_df$Predictor <- paste(int_df$name, 'BMI', sep='*')
      
      # Combine
      all_pred <- as.data.frame(rbind(main_df, int_df))
      betas <- all_pred %>% dplyr::rename(Beta=value)
      
      n_predictors <- length(unique(all_pred$name))
      cat(sprintf("    Predictors in model: %d\n", n_predictors))
      
      # Save betas
      out_beta <- paste0(output_dir, "/", tissue_name, "_", gene, "_betas.txt.gz")
      write.table(betas, file = gzfile(out_beta), sep = "\t",
                  col.names=T, row.names = F, quote=F)
      
      # Save model summary
      model_summary <- data.frame(
        Gene = gene,
        Tissue = tissue_name,
        Rsq = R_squared,
        p_Rsq = pval_lm,
        CV_Lambda = best_lambda,
        CV_MStEr = min_cvErrStd,
        CV_MSE = min_cvErr,
        R2_rep = rsq_rep,
        p_rep = pval_rep,
        n_predictors = n_predictors
      )
      
      out_summary <- paste0(output_dir, "/", tissue_name, "_", gene, "_model_summary.txt.gz")
      write.table(model_summary, file = gzfile(out_summary), sep = "\t",
                  col.names=T, row.names = F, quote=F)
      
      # Save model fit
      out_fit <- paste0(output_dir, "/", tissue_name, "_", gene, "_fit.rds")
      saveRDS(fit, file = out_fit)
      
      cat("    SUCCESS: Model saved\n")
      n_success <- n_success + 1
      
    } else {
      cat(sprintf("    Model R² (%.4f) below threshold (%.2f), not saving\n", R_squared, rsq_threshold))
      
      # Save failed model summary
      model_summary <- data.frame(
        Gene = gene,
        Tissue = tissue_name,
        Rsq = R_squared,
        p_Rsq = pval_lm,
        CV_Lambda = best_lambda,
        CV_MStEr = min_cvErrStd,
        CV_MSE = min_cvErr,
        R2_rep = NA,
        p_rep = NA,
        n_predictors = NA
      )
      
      out_summary <- paste0(output_dir, "/", tissue_name, "_", gene, "_model_summary.txt.gz")
      write.table(model_summary, file = gzfile(out_summary), sep = "\t",
                  col.names=T, row.names = F, quote=F)
      
      n_failed <- n_failed + 1
    }
    
  }, error = function(e) {
    cat(sprintf("    Error fitting model: %s\n", e$message))
    
    # Save failed model summary
    model_summary <- data.frame(
      Gene = gene,
      Tissue = tissue_name,
      Rsq = NA,
      p_Rsq = NA,
      CV_Lambda = NA,
      CV_MStEr = NA,
      CV_MSE = NA,
      R2_rep = NA,
      p_rep = NA,
      n_predictors = NA
    )
    
    out_summary <- paste0(output_dir, "/", tissue_name, "_", gene, "_model_summary.txt.gz")
    write.table(model_summary, file = gzfile(out_summary), sep = "\t",
                col.names=T, row.names = F, quote=F)
    
    n_failed <- n_failed + 1
  })
  
  cat("\n")
}

# ============================================================================
# SUMMARY
# ============================================================================

cat("================================================================================\n")
cat(" Summary\n")
cat("================================================================================\n")
cat(sprintf("  Chromosome:           %d\n", chr_here))
cat(sprintf("  Part:                 %d of %d\n", part_here, n_parts))
cat(sprintf("  Total genes:          %d\n", n_genes))
cat(sprintf("  Successful models:    %d (%.1f%%)\n", n_success, 100*n_success/n_genes))
cat(sprintf("  Failed models:        %d (%.1f%%)\n", n_failed, 100*n_failed/n_genes))
cat("================================================================================\n")
cat("\nDone!\n\n")
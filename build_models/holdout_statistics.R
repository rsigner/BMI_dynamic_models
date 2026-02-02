#!/usr/bin/Rscript

#' Calculate Holdout Statistics and Combine Model Results
#' 
#' This script combines glinternet model results across all genes, calculates
#' validation statistics on the held-out fold, filters for significant models,
#' and creates validation plots.
#'
#' @author Rebecca Signer
#' @date 01/31/2026

# ============================================================================
# LOAD REQUIRED LIBRARIES
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(optparse)
})

# ============================================================================
# PARSE COMMAND LINE ARGUMENTS
# ============================================================================

option_list <- list(
  make_option(c("--tissue"), type="character", default=NULL,
              help="Tissue name (e.g., Brain-Cortex, Adipose-Subcutaneous)",
              metavar="CHARACTER"),
  
  make_option(c("--gene-list-file"), type="character", default=NULL,
              help="Path to file listing all genes tested (one Ensembl ID per line)",
              metavar="PATH"),
  
  make_option(c("--model-dir"), type="character", default=NULL,
              help="Directory containing individual gene model files",
              metavar="PATH"),
  
  make_option(c("--feature-file"), type="character", default=NULL,
              help="Path to gene feature/annotation file",
              metavar="PATH"),
  
  make_option(c("--output-dir"), type="character", default=NULL,
              help="Output directory for combined results",
              metavar="PATH"),
  
  make_option(c("--plot-dir"), type="character", default=NULL,
              help="Output directory for plots [default: same as output-dir]",
              metavar="PATH"),
  
  make_option(c("--rsq-threshold"), type="numeric", default=0.01,
              help="Minimum R² threshold for significance [default: %default]",
              metavar="NUMERIC"),
  
  make_option(c("--pval-threshold"), type="numeric", default=0.05,
              help="Maximum p-value threshold for significance [default: %default]",
              metavar="NUMERIC")
)

opt_parser <- OptionParser(
  option_list=option_list,
  description="Calculate holdout statistics and combine glinternet model results",
  epilogue="Example:\n  Rscript holdout_statistics.R --tissue Brain-Cortex --gene-list-file genes_all.txt --model-dir models/ --feature-file features.txt.gz --output-dir results/ --plot-dir plots/"
)

opt <- parse_args(opt_parser)

# ============================================================================
# VALIDATE ARGUMENTS
# ============================================================================

validate_arguments <- function(opt) {
  errors <- c()
  
  # Check required arguments
  if (is.null(opt$tissue)) errors <- c(errors, "  --tissue is required")
  if (is.null(opt$`gene-list-file`)) errors <- c(errors, "  --gene-list-file is required")
  if (is.null(opt$`model-dir`)) errors <- c(errors, "  --model-dir is required")
  if (is.null(opt$`feature-file`)) errors <- c(errors, "  --feature-file is required")
  if (is.null(opt$`output-dir`)) errors <- c(errors, "  --output-dir is required")
  
  if (length(errors) > 0) {
    cat("Error: Missing required arguments:\n")
    cat(paste(errors, collapse="\n"), "\n\n")
    print_help(opt_parser)
    quit(status=1)
  }
  
  # Check files/directories exist
  if (!file.exists(opt$`gene-list-file`)) {
    cat(sprintf("Error: Gene list file not found: %s\n", opt$`gene-list-file`))
    quit(status=1)
  }
  
  if (!dir.exists(opt$`model-dir`)) {
    cat(sprintf("Error: Model directory not found: %s\n", opt$`model-dir`))
    quit(status=1)
  }
  
  if (!file.exists(opt$`feature-file`)) {
    cat(sprintf("Error: Feature file not found: %s\n", opt$`feature-file`))
    quit(status=1)
  }
  
  # Create output directories
  if (!dir.exists(opt$`output-dir`)) {
    dir.create(opt$`output-dir`, recursive=TRUE, showWarnings=FALSE)
    cat(sprintf("Created output directory: %s\n", opt$`output-dir`))
  }
  
  # Set plot dir to output dir if not specified
  if (is.null(opt$`plot-dir`)) {
    opt$`plot-dir` <- opt$`output-dir`
  }
  
  if (!dir.exists(opt$`plot-dir`)) {
    dir.create(opt$`plot-dir`, recursive=TRUE, showWarnings=FALSE)
    cat(sprintf("Created plot directory: %s\n", opt$`plot-dir`))
  }
  
  return(opt)
}

opt <- validate_arguments(opt)

# ============================================================================
# SETUP
# ============================================================================

tissue_name <- opt$tissue
gene_list_file <- opt$`gene-list-file`
model_dir <- opt$`model-dir`
feature_file <- opt$`feature-file`
output_dir <- opt$`output-dir`
plot_dir <- opt$`plot-dir`
rsq_threshold <- opt$`rsq-threshold`
pval_threshold <- opt$`pval-threshold`

cat("\n")
cat("================================================================================\n")
cat(" Combining Model Results and Creating Summary\n")
cat("================================================================================\n")
cat(sprintf("  Tissue:          %s\n", tissue_name))
cat(sprintf("  Gene list:       %s\n", gene_list_file))
cat(sprintf("  Model dir:       %s\n", model_dir))
cat(sprintf("  Feature file:    %s\n", feature_file))
cat(sprintf("  Output dir:      %s\n", output_dir))
cat(sprintf("  Plot dir:        %s\n", plot_dir))
cat(sprintf("  R² threshold:    %.2f\n", rsq_threshold))
cat(sprintf("  P-val threshold: %.2f\n", pval_threshold))
cat("================================================================================\n\n")

# ============================================================================
# LOAD GENE FEATURES
# ============================================================================

cat("Loading gene features...\n")
featuredt <- as.data.frame(fread(feature_file, header=T, sep="\t", stringsAsFactors=F))
cat(sprintf("  Loaded: %d genes\n", nrow(featuredt)))

# ============================================================================
# LOAD GENE LIST
# ============================================================================

cat("\nLoading gene list...\n")
genes_here <- read_lines(gene_list_file)
cat(sprintf("  Genes to process: %d\n", length(genes_here)))

# ============================================================================
# COMBINE MODEL SUMMARIES
# ============================================================================

cat("\nCombining model summaries...\n")

model_all <- NULL
n_found <- 0
n_missing <- 0
progress_interval <- max(1, floor(length(genes_here) / 10))

for (i in seq_along(genes_here)) {
  gene_here <- genes_here[i]
  
  # Show progress
  if (i %% progress_interval == 0) {
    cat(sprintf("  Progress: %d/%d genes (%.1f%%)\n", i, length(genes_here), 100*i/length(genes_here)))
  }
  
  # Construct model summary filename
  model_file <- sprintf("%s/%s_%s_model_summary.txt.gz", model_dir, tissue_name, gene_here)
  
  if (file.exists(model_file)) {
    model_here <- tryCatch({
      as.data.frame(fread(model_file, header=T, sep="\t", stringsAsFactors=F))
    }, error = function(e) {
      cat(sprintf("  Warning: Error reading %s: %s\n", model_file, e$message))
      return(NULL)
    })
    
    if (!is.null(model_here)) {
      model_all <- rbind(model_all, model_here)
      n_found <- n_found + 1
    }
  } else {
    n_missing <- n_missing + 1
  }
}

cat(sprintf("  Models found: %d\n", n_found))
cat(sprintf("  Models missing: %d\n", n_missing))

if (is.null(model_all) || nrow(model_all) == 0) {
  cat("\nError: No model summaries found!\n")
  cat("  Check that model files exist in the model directory.\n")
  quit(status=1)
}

# ============================================================================
# FILTER FOR SIGNIFICANT MODELS
# ============================================================================

cat("\nFiltering for significant models...\n")
cat(sprintf("  Criteria: R² > %.2f, P < %.2f (training and validation)\n", 
            rsq_threshold, pval_threshold))

model_sig <- model_all %>%
  filter(!is.na(n_predictors)) %>%
  filter(as.numeric(Rsq) > rsq_threshold) %>%
  filter(as.numeric(p_Rsq) < pval_threshold) %>%
  filter(as.numeric(R2_rep) > rsq_threshold) %>%
  filter(as.numeric(p_rep) < pval_threshold)

n_total <- nrow(model_all)
n_sig <- nrow(model_sig)

cat(sprintf("  Total genes tested:   %d\n", n_total))
cat(sprintf("  Significant models:   %d (%.1f%%)\n", n_sig, 100*n_sig/n_total))

if (n_sig == 0) {
  cat("\nWarning: No significant models found!\n")
  cat("  You may need to adjust thresholds or check model quality.\n")
  quit(status=0)
}

# ============================================================================
# ANNOTATE WITH GENE FEATURES
# ============================================================================

cat("\nAnnotating with gene features...\n")
sig_annot1 <- left_join(model_sig, featuredt, by=c('Gene'='Name'))

# ============================================================================
# LOAD AND COMBINE COEFFICIENTS
# ============================================================================

cat("\nLoading coefficients for significant models...\n")

sig_genes <- unique(sig_annot1$Gene)
cat(sprintf("  Genes with significant models: %d\n", length(sig_genes)))

coef_all <- NULL
n_coef_found <- 0
progress_interval_coef <- max(1, floor(length(sig_genes) / 10))

for (i in seq_along(sig_genes)) {
  gene_here <- sig_genes[i]
  
  # Show progress
  if (i %% progress_interval_coef == 0) {
    cat(sprintf("  Progress: %d/%d genes (%.1f%%)\n", i, length(sig_genes), 100*i/length(sig_genes)))
  }
  
  # Construct coefficient filename
  coef_file <- sprintf("%s/%s_%s_betas.txt.gz", model_dir, tissue_name, gene_here)
  
  if (file.exists(coef_file)) {
    coef_here <- tryCatch({
      as.data.frame(fread(coef_file, header=T, sep="\t", stringsAsFactors=F))
    }, error = function(e) {
      cat(sprintf("  Warning: Error reading %s: %s\n", coef_file, e$message))
      return(NULL)
    })
    
    if (!is.null(coef_here)) {
      coef_here$Gene <- gene_here
      coef_all <- rbind(coef_all, coef_here)
      n_coef_found <- n_coef_found + 1
    }
  }
}

cat(sprintf("  Coefficient files loaded: %d\n", n_coef_found))

# ============================================================================
# ANNOTATE COEFFICIENTS
# ============================================================================

cat("\nAnnotating coefficients...\n")

# Add gene features
coef_annot <- left_join(coef_all, featuredt, by=c('Gene'='Name'))
coef_annot$tissue <- tissue_name

# Count predictors per gene
n_predict_df <- coef_annot %>%
  group_by(tissue, Gene) %>%
  summarise(n_predictors = n(), .groups='drop')

# Join predictor counts
beta_int_annot <- left_join(coef_annot, n_predict_df, by=c('tissue', 'Gene'))

cat(sprintf("  Total predictors: %d\n", nrow(beta_int_annot)))
cat(sprintf("  Predictors per gene (mean): %.1f\n", mean(n_predict_df$n_predictors)))
cat(sprintf("  Predictors per gene (median): %.0f\n", median(n_predict_df$n_predictors)))

# Save annotated coefficients
out_betas <- sprintf("%s/%s_r2_%.2f_p_%.2f_betas.txt.gz", 
                     output_dir, tissue_name, rsq_threshold, pval_threshold)
write.table(beta_int_annot, file = gzfile(out_betas), sep = "\t", 
            col.names=T, row.names = F, quote=F)
cat(sprintf("\nSaved: %s\n", out_betas))

# ============================================================================
# CATEGORIZE PREDICTOR TYPES
# ============================================================================

cat("\nCategorizing predictor types...\n")

# Identify genes with interactions
interactions <- coef_all %>%
  filter(predictor_type == 'Interaction Effect') %>%
  .$Gene %>%
  unique()

cat(sprintf("  Genes with interactions: %d\n", length(interactions)))

# Identify genes with only main effects
mains <- setdiff(unique(sig_annot1$Gene), interactions)

# Check which main effect genes include BMI
main_df <- coef_all %>% filter(Gene %in% mains)
main_BMI <- main_df %>%
  filter(name == 'BMI') %>%
  .$Gene %>%
  unique()

cat(sprintf("  Genes with BMI main effect (no interaction): %d\n", length(main_BMI)))

# Genes with SNP main effects only
main_noBMI <- setdiff(mains, main_BMI)
cat(sprintf("  Genes with SNP main effects only: %d\n", length(main_noBMI)))

# Add predictor type category
sig_annot1$Predictor_types <- ifelse(
  sig_annot1$Gene %in% interactions,
  'Interaction',
  ifelse(sig_annot1$Gene %in% main_BMI,
         'Includes BMI as a Main Effect',
         'Has SNP Main Effects')
)

cat("\n  Summary:\n")
print(summary(factor(sig_annot1$Predictor_types)))

# Save annotated models
out_models <- sprintf("%s/%s_r2_%.2f_p_%.2f_models.txt.gz",
                      output_dir, tissue_name, rsq_threshold, pval_threshold)
write.table(sig_annot1, file = gzfile(out_models), sep = "\t",
            col.names=T, row.names = F, quote=F)
cat(sprintf("\nSaved: %s\n", out_models))

# ============================================================================
# CREATE VALIDATION PLOT
# ============================================================================

cat("\nCreating validation R² plot...\n")

# Reorder so interactions are plotted last (on top)
sig_annot1$order <- ifelse(
  sig_annot1$Gene %in% interactions, '3',
  ifelse(sig_annot1$Gene %in% main_BMI, '2', '1')
)
sig_annot1 <- sig_annot1 %>% arrange(order)

# Fit linear model
model_out <- summary(lm(R2_rep ~ Rsq, data=sig_annot1))
beta <- round(model_out$coefficients[2, 1], 3)
intercept <- round(model_out$coefficients[1, 1], 3)
rsq <- round(model_out$r.squared, 3)
pval <- model_out$coefficients[2, 4]

# Format annotation text
formula_syntax <- paste0('y = ', intercept, ' + ', beta, 'x\n')
rsq_syntax <- paste0('R² = ', rsq, '\n')
pval_syntax <- paste0('P = ', format(pval, scientific=TRUE, digits=3))
text_add <- paste0(formula_syntax, rsq_syntax, pval_syntax)

cat(sprintf("  Linear model: y = %.3f + %.3fx\n", intercept, beta))
cat(sprintf("  R² = %.3f, P = %.2e\n", rsq, pval))

# Create plot
pdf(sprintf("%s/%s_R2_holdout_%.2f_%.2f.pdf", plot_dir, tissue_name, rsq_threshold, pval_threshold))

p <- ggplot(sig_annot1, aes(x=Rsq, y=R2_rep, color=Predictor_types)) +
  geom_point(alpha=0.7, size=2) +
  geom_smooth(method = "lm", formula = y ~ x, color='black', linetype='dashed') +
  geom_abline(intercept=0, slope=1, linetype='dotted', color='gray50') +
  theme_classic(base_size=15) +
  annotate("text", x=-Inf, y=Inf, hjust=0, vjust=1, label=text_add, size=5) +
  xlim(0, 1) + ylim(0, 1) +
  scale_color_manual(
    values=c('#004D40', '#FFC107', '#D81B60'),
    name='Predictor Type'
  ) +
  labs(
    x='Training R² (3-fold CV)',
    y='Validation R² (Held-out Fold)',
    title=paste0(tissue_name, ': Model Validation')
  ) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill='white', color='black')
  )

print(p)
dev.off()

cat(sprintf("\nSaved: %s/%s_R2_holdout_%.2f_%.2f.pdf\n",
            plot_dir, tissue_name, rsq_threshold, pval_threshold))

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\n")
cat("================================================================================\n")
cat(" Summary Statistics\n")
cat("================================================================================\n")
cat(sprintf("  Tissue:                        %s\n", tissue_name))
cat(sprintf("  Total genes tested:            %d\n", n_total))
cat(sprintf("  Significant models:            %d (%.1f%%)\n", n_sig, 100*n_sig/n_total))
cat("\n")
cat("  By predictor type:\n")
cat(sprintf("    Interaction models:          %d (%.1f%%)\n",
            length(interactions), 100*length(interactions)/n_sig))
cat(sprintf("    BMI main effect models:      %d (%.1f%%)\n",
            length(main_BMI), 100*length(main_BMI)/n_sig))
cat(sprintf("    SNP-only models:             %d (%.1f%%)\n",
            length(main_noBMI), 100*length(main_noBMI)/n_sig))
cat("\n")
cat("  Model performance (training):\n")
cat(sprintf("    Mean R²:                     %.3f\n", mean(sig_annot1$Rsq)))
cat(sprintf("    Median R²:                   %.3f\n", median(sig_annot1$Rsq)))
cat(sprintf("    Range:                       %.3f - %.3f\n",
            min(sig_annot1$Rsq), max(sig_annot1$Rsq)))
cat("\n")
cat("  Model performance (validation):\n")
cat(sprintf("    Mean R²:                     %.3f\n", mean(sig_annot1$R2_rep)))
cat(sprintf("    Median R²:                   %.3f\n", median(sig_annot1$R2_rep)))
cat(sprintf("    Range:                       %.3f - %.3f\n",
            min(sig_annot1$R2_rep), max(sig_annot1$R2_rep)))
cat("\n")
cat("  Predictors per model:\n")
cat(sprintf("    Mean:                        %.1f\n", mean(sig_annot1$n_predictors)))
cat(sprintf("    Median:                      %.0f\n", median(sig_annot1$n_predictors)))
cat(sprintf("    Range:                       %d - %d\n",
            min(sig_annot1$n_predictors), max(sig_annot1$n_predictors)))
cat("================================================================================\n")

# ============================================================================
# OUTPUT FILES SUMMARY
# ============================================================================

cat("\n")
cat("Output files:\n")
cat(sprintf("  1. Model summaries:  %s\n", out_models))
cat(sprintf("  2. All coefficients: %s\n", out_betas))
cat(sprintf("  3. Validation plot:  %s/%s_R2_holdout_%.2f_%.2f.pdf\n",
            plot_dir, tissue_name, rsq_threshold, pval_threshold))
cat("\n")
cat("Done!\n\n")
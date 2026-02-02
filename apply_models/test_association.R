#!/usr/bin/Rscript

#' Test Association Between Predicted Expression and Phenotype
#' 
#' This script tests for associations between predicted gene expression (GREX)
#' and a phenotype of interest. Supports both case/control (logistic regression)
#' and quantitative traits (linear regression).
#' 
#' Uses parts system to match prediction output and enable parallel processing.
#'
#' @author Rebecca Signer
#' @date 01/31/26

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
              help="Tissue name (e.g., Adipose-Subcutaneous)",
              metavar="CHARACTER"),
  
  make_option(c("--chr"), type="integer", default=NULL,
              help="Chromosome number (1-22)",
              metavar="INTEGER"),
  
  make_option(c("--part"), type="integer", default=NULL,
              help="Part number (1-10)",
              metavar="INTEGER"),
  
  make_option(c("--expression-dir"), type="character", default=NULL,
              help="Directory containing predicted expression files",
              metavar="PATH"),
  
  make_option(c("--phenotype-file"), type="character", default=NULL,
              help="Path to phenotype file (columns: SUBJID, phenotype)",
              metavar="PATH"),
  
  make_option(c("--covariate-file"), type="character", default=NULL,
              help="Path to covariate file (columns: SUBJID, cov1, cov2, ...)",
              metavar="PATH"),
  
  make_option(c("--output-dir"), type="character", default=NULL,
              help="Output directory for results",
              metavar="PATH"),
  
  make_option(c("--phenotype-type"), type="character", default="binary",
              help="Phenotype type: 'binary' for case/control or 'quantitative' [default: %default]",
              metavar="TYPE"),
  
  make_option(c("--phenotype-name"), type="character", default=NULL,
              help="Name of phenotype column in phenotype file",
              metavar="NAME"),
  
  make_option(c("--covariates"), type="character", default=NULL,
              help="Comma-separated list of covariate names (e.g., 'Age,Sex,PC1,PC2,PC3')",
              metavar="LIST"),
  
  make_option(c("--test-with-bmi"), type="logical", default=FALSE,
              help="If TRUE, includes BMI as a covariate in the model [default: %default]",
              metavar="LOGICAL")
)

opt_parser <- OptionParser(
  option_list=option_list,
  description="Test association between predicted expression and phenotype",
  epilogue="Example:\n  Rscript test_association.R --tissue Adipose-Subcutaneous --chr 1 --part 1 --expression-dir results/ --phenotype-file phenotypes.txt --covariate-file covariates.txt --output-dir associations/ --phenotype-type binary --phenotype-name disease_status --covariates 'Age,Sex,PC1,PC2,PC3,PC4,PC5'"
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
  if (is.null(opt$`expression-dir`)) errors <- c(errors, "  --expression-dir is required")
  if (is.null(opt$`phenotype-file`)) errors <- c(errors, "  --phenotype-file is required")
  if (is.null(opt$`output-dir`)) errors <- c(errors, "  --output-dir is required")
  if (is.null(opt$`phenotype-name`)) errors <- c(errors, "  --phenotype-name is required")
  
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
  if (opt$part < 1 || opt$part > 10) {
    cat("Error: Part must be between 1 and 10\n")
    quit(status=1)
  }
  
  # Check files/directories exist
  if (!dir.exists(opt$`expression-dir`)) {
    cat(sprintf("Error: Expression directory not found: %s\n", opt$`expression-dir`))
    quit(status=1)
  }
  
  if (!file.exists(opt$`phenotype-file`)) {
    cat(sprintf("Error: Phenotype file not found: %s\n", opt$`phenotype-file`))
    quit(status=1)
  }
  
  if (!is.null(opt$`covariate-file`) && !file.exists(opt$`covariate-file`)) {
    cat(sprintf("Error: Covariate file not found: %s\n", opt$`covariate-file`))
    quit(status=1)
  }
  
  # Validate phenotype type
  if (!(opt$`phenotype-type` %in% c("binary", "quantitative"))) {
    cat("Error: --phenotype-type must be 'binary' or 'quantitative'\n")
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
expression_dir <- opt$`expression-dir`
phenotype_file <- opt$`phenotype-file`
covariate_file <- opt$`covariate-file`
output_dir <- opt$`output-dir`
phenotype_type <- opt$`phenotype-type`
phenotype_name <- opt$`phenotype-name`
covariate_names <- if (!is.null(opt$covariates)) strsplit(opt$covariates, ",")[[1]] else NULL
test_with_bmi <- opt$`test-with-bmi`

cat("\n")
cat("================================================================================\n")
cat(" Testing Association: Predicted Expression vs Phenotype\n")
cat("================================================================================\n")
cat(sprintf("  Tissue:           %s\n", tissue_name))
cat(sprintf("  Chromosome:       %d\n", chr_here))
cat(sprintf("  Part:             %d\n", part_here))
cat(sprintf("  Expression dir:   %s\n", expression_dir))
cat(sprintf("  Phenotype file:   %s\n", phenotype_file))
cat(sprintf("  Covariate file:   %s\n", ifelse(is.null(covariate_file), "None", covariate_file)))
cat(sprintf("  Output dir:       %s\n", output_dir))
cat(sprintf("  Phenotype type:   %s\n", phenotype_type))
cat(sprintf("  Phenotype name:   %s\n", phenotype_name))
cat(sprintf("  Covariates:       %s\n", ifelse(is.null(covariate_names), "None", paste(covariate_names, collapse=", "))))
cat(sprintf("  Include BMI:      %s\n", test_with_bmi))
cat("================================================================================\n\n")

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading data...\n")

# Load predicted expression for this chromosome and part
cat("  Loading predicted expression...\n")
expression_file <- sprintf("%s/%s_%d_%d_BGREX.txt.gz", expression_dir, tissue_name, chr_here, part_here)

if (!file.exists(expression_file)) {
  cat(sprintf("Error: Expression file not found: %s\n", expression_file))
  cat("  Expected format: {tissue}_{chr}_{part}_BGREX.txt.gz\n")
  quit(status=1)
}

bgrex_df <- tryCatch({
  as.data.frame(fread(expression_file, header=TRUE, sep="\t", stringsAsFactors=FALSE))
}, error = function(e) {
  cat(sprintf("Error reading expression file: %s\n", e$message))
  quit(status=1)
})

cat(sprintf("    Loaded: %d individuals, %d genes\n", nrow(bgrex_df), ncol(bgrex_df)-1))

# Load phenotype
cat("  Loading phenotype data...\n")
pheno_df <- tryCatch({
  as.data.frame(fread(phenotype_file, header=TRUE, sep="\t", stringsAsFactors=FALSE))
}, error = function(e) {
  cat(sprintf("Error reading phenotype file: %s\n", e$message))
  quit(status=1)
})

if (!("SUBJID" %in% colnames(pheno_df))) {
  cat("Error: Phenotype file must contain 'SUBJID' column\n")
  cat(sprintf("  Found columns: %s\n", paste(colnames(pheno_df), collapse=", ")))
  quit(status=1)
}

if (!(phenotype_name %in% colnames(pheno_df))) {
  cat(sprintf("Error: Phenotype '%s' not found in phenotype file\n", phenotype_name))
  cat(sprintf("  Available columns: %s\n", paste(colnames(pheno_df), collapse=", ")))
  quit(status=1)
}

cat(sprintf("    Loaded: %d individuals\n", nrow(pheno_df)))

# Load covariates (if provided)
if (!is.null(covariate_file)) {
  cat("  Loading covariates...\n")
  cov_df <- tryCatch({
    as.data.frame(fread(covariate_file, header=TRUE, sep="\t", stringsAsFactors=FALSE))
  }, error = function(e) {
    cat(sprintf("Error reading covariate file: %s\n", e$message))
    quit(status=1)
  })
  
  if (!("SUBJID" %in% colnames(cov_df))) {
    cat("Error: Covariate file must contain 'SUBJID' column\n")
    cat(sprintf("  Found columns: %s\n", paste(colnames(cov_df), collapse=", ")))
    quit(status=1)
  }
  
  # Check all requested covariates exist
  missing_covs <- setdiff(covariate_names, colnames(cov_df))
  if (length(missing_covs) > 0) {
    cat("Error: Covariates not found in covariate file:\n")
    cat(sprintf("  %s\n", paste(missing_covs, collapse=", ")))
    quit(status=1)
  }
  
  cat(sprintf("    Loaded: %d individuals, %d covariates\n", nrow(cov_df), length(covariate_names)))
}

# ============================================================================
# MERGE DATA
# ============================================================================

cat("Merging datasets...\n")

# Merge phenotype with expression
dat <- bgrex_df %>%
  inner_join(pheno_df, by = "SUBJID")

cat(sprintf("  After merging expression + phenotype: %d individuals\n", nrow(dat)))

# Merge with covariates if provided
if (!is.null(covariate_file)) {
  dat <- dat %>%
    inner_join(cov_df %>% select(SUBJID, all_of(covariate_names)), by = "SUBJID")
  
  cat(sprintf("  After merging covariates: %d individuals\n", nrow(dat)))
}

# ============================================================================
# PREPARE FOR TESTING
# ============================================================================

cat("Preparing for association testing...\n")

# Get list of genes
genes <- setdiff(colnames(bgrex_df), "SUBJID")
cat(sprintf("  Testing %d genes (chr %d, part %d)\n", length(genes), chr_here, part_here))

# Build model formula
if (!is.null(covariate_names)) {
  cov_syntax <- paste(covariate_names, collapse=" + ")
  if (test_with_bmi && "BMI" %in% colnames(dat)) {
    model_formula <- paste("phenotype ~ BGREX + BMI +", cov_syntax)
  } else {
    model_formula <- paste("phenotype ~ BGREX +", cov_syntax)
  }
} else {
  if (test_with_bmi && "BMI" %in% colnames(dat)) {
    model_formula <- "phenotype ~ BGREX + BMI"
  } else {
    model_formula <- "phenotype ~ BGREX"
  }
}

cat(sprintf("  Model formula: %s\n", model_formula))

# Prepare phenotype
if (phenotype_type == "binary") {
  # For binary phenotype, ensure it's a factor
  unique_vals <- unique(dat[[phenotype_name]])
  if (length(unique_vals) != 2) {
    cat(sprintf("Warning: Binary phenotype has %d unique values: %s\n", 
                length(unique_vals), paste(unique_vals, collapse=", ")))
  }
  
  # Convert to factor (0/1 or Case/Control)
  dat$phenotype <- factor(dat[[phenotype_name]])
  cat(sprintf("  Phenotype levels: %s (reference: %s)\n", 
              paste(levels(dat$phenotype), collapse=", "),
              levels(dat$phenotype)[1]))
  
  n_case <- sum(dat$phenotype == levels(dat$phenotype)[2], na.rm=TRUE)
  n_control <- sum(dat$phenotype == levels(dat$phenotype)[1], na.rm=TRUE)
  cat(sprintf("  Sample size: %d cases, %d controls\n", n_case, n_control))
  
} else {
  # For quantitative phenotype
  dat$phenotype <- as.numeric(dat[[phenotype_name]])
  cat(sprintf("  Phenotype range: %.2f to %.2f\n", 
              min(dat$phenotype, na.rm=TRUE),
              max(dat$phenotype, na.rm=TRUE)))
}

# ============================================================================
# RUN ASSOCIATION TESTS
# ============================================================================

cat("\nRunning association tests...\n")

all_results <- NULL
progress_interval <- max(1, floor(length(genes) / 10))

for (i in seq_along(genes)) {
  gene <- genes[i]
  
  # Show progress
  if (i %% progress_interval == 0) {
    cat(sprintf("  Progress: %d/%d genes (%.1f%%)\n", i, length(genes), 100*i/length(genes)))
  }
  
  # Prepare data for this gene
  dat_gene <- dat %>%
    mutate(BGREX = .data[[gene]]) %>%
    filter(!is.na(BGREX))
  
  # Skip if all GREX values are NA
  if (nrow(dat_gene) == 0 || all(is.na(dat_gene$BGREX))) {
    result_row <- data.frame(
      Gene = gene,
      Beta = NA,
      SE = NA,
      P = NA,
      N = 0
    )
  } else {
    # Run regression
    result_row <- tryCatch({
      if (phenotype_type == "binary") {
        # Logistic regression
        model <- glm(as.formula(model_formula), family = "binomial", data = dat_gene)
      } else {
        # Linear regression
        model <- lm(as.formula(model_formula), data = dat_gene)
      }
      
      # Extract results
      summary_model <- summary(model)
      
      if ("BGREX" %in% rownames(summary_model$coefficients)) {
        bgrex_coef <- summary_model$coefficients["BGREX", ]
        
        if (phenotype_type == "binary") {
          data.frame(
            Gene = gene,
            Beta = bgrex_coef["Estimate"],
            SE = bgrex_coef["Std. Error"],
            P = bgrex_coef["Pr(>|z|)"],
            N = nrow(dat_gene)
          )
        } else {
          data.frame(
            Gene = gene,
            Beta = bgrex_coef["Estimate"],
            SE = bgrex_coef["Std. Error"],
            P = bgrex_coef["Pr(>|t|)"],
            N = nrow(dat_gene)
          )
        }
      } else {
        data.frame(
          Gene = gene,
          Beta = NA,
          SE = NA,
          P = NA,
          N = nrow(dat_gene)
        )
      }
    }, error = function(e) {
      data.frame(
        Gene = gene,
        Beta = NA,
        SE = NA,
        P = NA,
        N = nrow(dat_gene)
      )
    })
  }
  
  all_results <- rbind(all_results, result_row)
}

cat(sprintf("  Completed: %d/%d genes\n", length(genes), length(genes)))

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("\nSaving results...\n")

# Output filename includes chromosome and part
output_file <- sprintf("%s/%s_%d_%d_associations.txt.gz", 
                       output_dir, tissue_name, chr_here, part_here)

write.table(all_results,
            file = gzfile(output_file),
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

cat(sprintf("  Saved: %s\n", output_file))

# ============================================================================
# SUMMARY
# ============================================================================

cat("\nSummary:\n")
cat("================================================================================\n")
cat(sprintf("  Chromosome:               %d\n", chr_here))
cat(sprintf("  Part:                     %d\n", part_here))
cat(sprintf("  Total genes tested:       %d\n", nrow(all_results)))
cat(sprintf("  Genes with results:       %d\n", sum(!is.na(all_results$P))))
cat(sprintf("  Significant (P < 0.05):   %d\n", sum(all_results$P < 0.05, na.rm=TRUE)))

if (sum(!is.na(all_results$P)) > 0) {
  cat(sprintf("  Minimum P-value:          %.2e\n", min(all_results$P, na.rm=TRUE)))
}
cat("================================================================================\n")
cat("\nNote: Use combine_results.R to join all parts and calculate genome-wide FDR.\n")
cat("Done!\n\n")

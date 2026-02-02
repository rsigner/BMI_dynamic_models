#!/usr/bin/Rscript

#' Combine Association Results and Calculate Genome-Wide FDR
#' 
#' This script combines association results across all chromosomes and parts,
#' then calculates genome-wide False Discovery Rate (FDR).
#' 
#' CRITICAL: FDR must be calculated on ALL tests genome-wide to be valid.
#' Per-chromosome or per-part FDR is incorrect and overly liberal.
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
              help="Tissue name (e.g., Adipose-Subcutaneous)",
              metavar="CHARACTER"),
  
  make_option(c("--input-dir"), type="character", default=NULL,
              help="Directory containing association result files",
              metavar="PATH"),
  
  make_option(c("--output-file"), type="character", default=NULL,
              help="Output file for combined results with genome-wide FDR",
              metavar="PATH"),
  
  make_option(c("--n-chr"), type="integer", default=22,
              help="Number of chromosomes to combine [default: %default]",
              metavar="INTEGER"),
  
  make_option(c("--n-parts"), type="integer", default=10,
              help="Number of parts per chromosome [default: %default]",
              metavar="INTEGER")
)

opt_parser <- OptionParser(
  option_list=option_list,
  description="Combine association results and calculate genome-wide FDR",
  epilogue="Example:\n  Rscript combine_results.R --tissue Adipose-Subcutaneous --input-dir associations/ --output-file associations/Adipose-Subcutaneous_all_associations_FDR.txt.gz"
)

opt <- parse_args(opt_parser)

# ============================================================================
# VALIDATE ARGUMENTS
# ============================================================================

validate_arguments <- function(opt) {
  errors <- c()
  
  # Check required arguments
  if (is.null(opt$tissue)) errors <- c(errors, "  --tissue is required")
  if (is.null(opt$`input-dir`)) errors <- c(errors, "  --input-dir is required")
  if (is.null(opt$`output-file`)) errors <- c(errors, "  --output-file is required")
  
  if (length(errors) > 0) {
    cat("Error: Missing required arguments:\n")
    cat(paste(errors, collapse="\n"), "\n\n")
    print_help(opt_parser)
    quit(status=1)
  }
  
  # Check input directory exists
  if (!dir.exists(opt$`input-dir`)) {
    cat(sprintf("Error: Input directory not found: %s\n", opt$`input-dir`))
    quit(status=1)
  }
}

validate_arguments(opt)

# ============================================================================
# SETUP
# ============================================================================

tissue_name <- opt$tissue
input_dir <- opt$`input-dir`
output_file <- opt$`output-file`
n_chr <- opt$`n-chr`
n_parts <- opt$`n-parts`

cat("\n")
cat("================================================================================\n")
cat(" Combining Association Results and Calculating Genome-Wide FDR\n")
cat("================================================================================\n")
cat(sprintf("  Tissue:        %s\n", tissue_name))
cat(sprintf("  Input dir:     %s\n", input_dir))
cat(sprintf("  Output file:   %s\n", output_file))
cat(sprintf("  Chromosomes:   %d\n", n_chr))
cat(sprintf("  Parts/chr:     %d\n", n_parts))
cat("================================================================================\n\n")

# ============================================================================
# COMBINE RESULTS ACROSS ALL CHROMOSOMES AND PARTS
# ============================================================================

cat("Combining association results...\n")

all_results <- NULL
n_files_found <- 0
n_files_missing <- 0
missing_files <- c()

for (chr in 1:n_chr) {
  for (part in 1:n_parts) {
    # Construct expected filename
    input_file <- sprintf("%s/%s_%d_%d_associations.txt.gz", 
                         input_dir, tissue_name, chr, part)
    
    if (file.exists(input_file)) {
      # Read file
      results_here <- tryCatch({
        as.data.frame(fread(input_file, header=TRUE, sep="\t", stringsAsFactors=FALSE))
      }, error = function(e) {
        cat(sprintf("Warning: Error reading %s: %s\n", input_file, e$message))
        return(NULL)
      })
      
      if (!is.null(results_here)) {
        # Add chromosome and part information
        results_here$Chr <- chr
        results_here$Part <- part
        
        # Combine with all results
        all_results <- rbind(all_results, results_here)
        n_files_found <- n_files_found + 1
        
        if (n_files_found %% 10 == 0) {
          cat(sprintf("  Processed %d files...\n", n_files_found))
        }
      }
    } else {
      n_files_missing <- n_files_missing + 1
      missing_files <- c(missing_files, sprintf("Chr%d_Part%d", chr, part))
    }
  }
}

cat(sprintf("  Files found: %d\n", n_files_found))
cat(sprintf("  Files missing: %d\n", n_files_missing))

if (n_files_missing > 0) {
  cat("\nMissing files:\n")
  cat(sprintf("  %s\n", paste(head(missing_files, 10), collapse=", ")))
  if (n_files_missing > 10) {
    cat(sprintf("  ... and %d more\n", n_files_missing - 10))
  }
}

if (is.null(all_results) || nrow(all_results) == 0) {
  cat("\nError: No results found to combine!\n")
  cat("  Check that:\n")
  cat(sprintf("    1. Input directory is correct: %s\n", input_dir))
  cat(sprintf("    2. Files follow naming: %s_CHR_PART_associations.txt.gz\n", tissue_name))
  cat("    3. Association testing completed successfully\n")
  quit(status=1)
}

cat(sprintf("\nTotal genes combined: %d\n", nrow(all_results)))

# ============================================================================
# CALCULATE GENOME-WIDE FDR
# ============================================================================

cat("\nCalculating genome-wide FDR...\n")

# Remove duplicate genes (shouldn't happen, but just in case)
n_before <- nrow(all_results)
all_results <- all_results %>% distinct(Gene, .keep_all = TRUE)
n_after <- nrow(all_results)

if (n_before != n_after) {
  cat(sprintf("  Warning: Removed %d duplicate gene entries\n", n_before - n_after))
}

# Calculate genome-wide FDR
all_results$FDR <- p.adjust(all_results$P, method = "fdr")

cat(sprintf("  Calculated FDR for %d genes\n", nrow(all_results)))

# ============================================================================
# SAVE COMBINED RESULTS
# ============================================================================

cat(sprintf("\nSaving combined results to: %s\n", output_file))

# Create output directory if needed
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
}

write.table(all_results,
            file = gzfile(output_file),
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\n")
cat("================================================================================\n")
cat(" Within-tissue Transcriptome-Wide Association Summary\n")
cat("================================================================================\n")
cat(sprintf("  Total genes tested:              %d\n", nrow(all_results)))
cat(sprintf("  Genes with valid P-values:       %d (%.1f%%)\n", 
            sum(!is.na(all_results$P)),
            100 * sum(!is.na(all_results$P)) / nrow(all_results)))
cat(sprintf("  Genes with missing P-values:     %d (%.1f%%)\n",
            sum(is.na(all_results$P)),
            100 * sum(is.na(all_results$P)) / nrow(all_results)))
cat("\n")
cat(sprintf("  Significant at P < 0.05:         %d (%.1f%%)\n",
            sum(all_results$P < 0.05, na.rm=TRUE),
            100 * sum(all_results$P < 0.05, na.rm=TRUE) / sum(!is.na(all_results$P))))
cat(sprintf("  Significant at P < 0.01:         %d (%.1f%%)\n",
            sum(all_results$P < 0.01, na.rm=TRUE),
            100 * sum(all_results$P < 0.01, na.rm=TRUE) / sum(!is.na(all_results$P))))
cat(sprintf("  Significant at P < 0.001:        %d (%.1f%%)\n",
            sum(all_results$P < 0.001, na.rm=TRUE),
            100 * sum(all_results$P < 0.001, na.rm=TRUE) / sum(!is.na(all_results$P))))
cat("\n")
cat(sprintf("  Significant at FDR < 0.05:       %d (%.1f%%)\n",
            sum(all_results$FDR < 0.05, na.rm=TRUE),
            100 * sum(all_results$FDR < 0.05, na.rm=TRUE) / sum(!is.na(all_results$P))))
cat(sprintf("  Significant at FDR < 0.10:       %d (%.1f%%)\n",
            sum(all_results$FDR < 0.10, na.rm=TRUE),
            100 * sum(all_results$FDR < 0.10, na.rm=TRUE) / sum(!is.na(all_results$P))))
cat(sprintf("  Significant at FDR < 0.20:       %d (%.1f%%)\n",
            sum(all_results$FDR < 0.20, na.rm=TRUE),
            100 * sum(all_results$FDR < 0.20, na.rm=TRUE) / sum(!is.na(all_results$P))))
cat("\n")

if (sum(!is.na(all_results$P)) > 0) {
  cat(sprintf("  Minimum P-value:                 %.2e\n", min(all_results$P, na.rm=TRUE)))
  cat(sprintf("  Minimum FDR:                     %.4f\n", min(all_results$FDR, na.rm=TRUE)))
}

cat("================================================================================\n")

# ============================================================================
# OPTIONAL: SAVE TOP HITS
# ============================================================================

# Save top hits by P-value
top_hits_p <- all_results %>%
  filter(!is.na(P)) %>%
  arrange(P) %>%
  head(50)

top_hits_file <- gsub("\\.txt\\.gz$", "_top50_byP.txt", output_file)
write.table(top_hits_p,
            file = top_hits_file,
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

cat(sprintf("\nSaved top 50 genes (by P-value): %s\n", top_hits_file))

# Save significant hits by FDR
sig_hits_fdr <- all_results %>%
  filter(FDR < 0.05) %>%
  arrange(FDR)

if (nrow(sig_hits_fdr) > 0) {
  sig_hits_file <- gsub("\\.txt\\.gz$", "_significant_FDR0.05.txt", output_file)
  write.table(sig_hits_fdr,
              file = sig_hits_file,
              sep = "\t",
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE)
  
  cat(sprintf("Saved significant genes (FDR < 0.05): %s\n", sig_hits_file))
} else {
  cat("\nNo genes significant at FDR < 0.05\n")
}

cat("\nDone!\n\n")
# Frequently Asked Questions (FAQ)

## General Questions

### What are BMI-dynamic models?

BMI-dynamic models predict gene expression using genotype data while accounting for Body Mass Index (BMI) as a modulating factor. Unlike traditional expression prediction models that assume genetic effects are constant, these models allow SNP effects to vary with BMI through interaction terms.

### What data do I need?

**To apply pre-built models:**

- Genotype data (imputed or WGS)
- BMI measurements
- Phenotype of interest (for association testing)

**To build new models:**

- Paired genotype and gene expression data
- BMI measurements

### Which GTEx tissues are available?

Pre-trained models are available for 8 GTEx v8 tissues:

- Adipose-Subcutaneous
- Brain-Cerebellum
- Brain-Cortex
- Brain-NucleusAccumbens
- Colon-Transverse
- Esophagus-Mucosa
- Esophagus-Muscularis
- Nerve-Tibial

------

## Data Format Questions

### What genotype format is required?

**For applying models:**

- PLINK `.raw` format (dosage or hard calls)
- Column names must be GTEx variant IDs: `chr{N}_{position}_{ref}_{alt}_b38`
- First column must be named `IID`
- Tab-separated, can be gzipped

**For building models:**

- Same format as above
- Must include all SNPs in ±1 Mb windows around genes
- Quality filters: MAF > 0.05, no missing values, HWE

### Why do column names need to be in GTEx format?

The models were trained using GTEx variant IDs. Your SNPs must be mapped to the same naming convention so the models can identify which SNPs to use. The mapping file `models/mapping/model_snps_mapping.txt.gz` helps convert rsIDs or hg19 coordinates to GTEx format.

### My genotypes are in hg19. What should I do?

1. Use the mapping file to identify positions in hg19 see `models/Readme.md` and `models/mapping/${tissue}_model_variants.txt.gz`
2. See `apply_models/README.md` Step 1 for detailed instructions

### Can I use hard-called genotypes instead of dosages?

Yes! Hard-called genotypes (0, 1, 2) work fine. The models were trained on dosages but accept hard calls.

------

## Technical Questions

### Why are results split into parts?

**Memory management**: Large biobank cohorts with thousands of samples require splitting chromosomes into parts to avoid memory issues.

**Parallelization**: Parts can be run simultaneously on a cluster, dramatically reducing wall-clock time.

**For predictions**: 10 parts per chromosome (220 jobs total) **For model building**: 30 parts per chromosome (660 jobs total)

### How many samples do I need?

**To apply models**: No minimum, but statistical power for association testing increases with sample size. Typical biobanks have 10,000+ samples.

### What is BGREX?

**BGREX** = **B**MI-aware **G**enetically **R**egulated **Ex**pression

This distinguishes our BMI-aware predictions from standard GREX (Genetically Regulated Expression) which doesn't account for BMI.

------

## Results Interpretation

### Should I test with or without BMI as a covariate?

**Both!** They answer different questions:

**Without BMI covariate**: Does BMI-aware genetically-predicted expression (which includes BMI effects) associate with disease?

**With BMI covariate**: Does BMI-aware genetically-predicted expression associate with disease while adjusting for the direct effects of BMI on the disease

### What's the difference between training and validation R²?

**Training R²**: How well the model predicts expression in the 3 folds used for training (with cross-validation)

**Validation R²**: How well the model predicts expression in the held-out 4th fold (never seen during training)

Validation R² is the true measure of model performance. We require both R² > 0.01 and P < 0.05 for a model to be considered significant.

### How do I interpret interaction effects?

A SNP × BMI interaction means the SNP's effect on gene expression changes with BMI. For example:

- SNP may increase expression at high BMI but decrease it at low BMI
- SNP effect may be significant at only high BMI and vice versa
- Genetic regulatory capacity of the SNP is BMI-dependent
- See the supplement of Signer 2026 for the full directionality decision tree

### Why don't all genes have predictions?

Genes may lack predictions because:

- Model R² < 0.01 (poor predictive performance)
- Model didn't validate on held-out fold
- Gene not expressed in that tissue

------

## Troubleshooting

### "Error: File not found"

**Check:**

1. File path is correct (use absolute paths)
2. File exists and is readable
3. File is gzipped if script expects `.gz` extension

### "Error: Column not found"

**Check:**

1. Required columns are present (SUBJID, BMI, IID, etc.)
2. Column names match exactly (case-sensitive)
3. File has header row
4. No extra spaces in column names

### "Memory allocation failed"

**Solutions:**

1. Request more memory in cluster job
2. Reduce chromosome into more parts
3. Process fewer genes per job
4. Use a machine with more RAM

### "No SNPs in cis-window"

This is expected for some genes. Possible reasons:

- Gene on chromosome not in your genotype file
- SNPs filtered out due to QC (MAF, missingness, etc.)
- Genotype file doesn't cover that genomic region

### "Model convergence issues"

For model building, convergence issues can occur when:

- Too few samples
- Too many predictors (SNPs)
- Multicollinearity among SNPs
- Very low variance in expression

**Solution**: These genes are automatically marked as failed and skipped.

### Results are all NA

**Check:**

1. Sample IDs match across files (SUBJID/IID)
2. Gene IDs match expression file column names
3. Genotype SNPs are in correct format
4. No systematic data issues (all zero, all NA, etc.)

------

### Should I keep intermediate files?

**Keep**:

- CV fold assignments
- Final combined results
- Significant models

**Can delete after combining**:

- Individual chromosome-part predictions
- Individual chromosome-part associations
- Individual gene model files (keep combined betas)

------

## Getting Help

### Script isn't working. What should I do?

1. Check this FAQ
2. Review the relevant README (`apply_models/` or `build_models/`)
3. Check error messages carefully
4. Verify file formats match exactly
5. Try with a small test dataset
6. Open a GitHub issue with:
   - Error message
   - Command you ran
   - Description of your data
   - Steps you've tried

### How do I cite this work?

```
[Citation will be added after publication]
```

### Can I modify the code for my needs?

Yes! This is MIT licensed - you can modify and redistribute. If you make improvements, please consider contributing back via pull request.

------

## Future Development

### How can I contribute?

See `CONTRIBUTING.md` for guidelines on:

- Reporting bugs
- Suggesting features
- Submitting pull requests
- Adding documentation
# Applying Pre-Built BMI-Dynamic Models to Biobank Data

## Overview

This workflow allows you to **predict gene expression** in biobank cohorts using pre-built BMI-dynamic models. The models were built using paired DNA and RNA data from GTEx and can be **applied** to biobanks that have only DNA data.

**What the models accept**:

- Genotype data (DNA)
- BMI values

**What the models predict**:

- Tissue-specific BMI-dynamic gene expression across thousands of genes

**What you can do with predictions**:

- Test association of  BMI-dynamic predicted gene expression with disease status in large biobanks
- Test association with quantitative traits

**Original validation**: We applied these models to imputed genotype files from the UK Biobank and Mount Sinai BioMe biobank.

------

## Pipeline Overview

```
Your Biobank Genotype Data
         ↓
Step 1: Prepare Genotype Files
         ↓
Step 2: Prepare BMI File  
         ↓
Step 3: Predict Expression (by chromosome and part)
         ↓
Step 4: Test Associations (by chromosome and part)
         ↓
Step 5: Assign Genome-Wide Significance
         ↓
     Results!
```

**Why parts?** Large biobank cohorts require splitting each chromosome into 10 parts for memory management and parallel processing.

------

## Step 1: Prepare the Genotype Files

### Required Final Format

The final formatted file needs to be a **data frame** with:

- **Rows**: Individuals from your biobank
- **Columns**: SNPs from the BMI-dynamic models
- **Column names**: GTEx variant IDs (format: `chr{N}_{position}_{ref}_{alt}_b38`)
- **Cells**: Values from 0-2 representing the number of alternate alleles
  - Can be hard calls from WGS (0, 1, or 2)
  - Or dosages from genotype imputation (0.0 to 2.0)
- **Format**: Tab-separated .raw file (can be gzipped)
- **Required column**: `IID` containing individual identifiers

### Detailed Preparation Steps

The next steps depend on the genetic data available in your biobank. In our original study, we used imputed genotype files from the UK Biobank and Mount Sinai BioMe biobank.

#### 1.1 Map Your Biobank SNP IDs to BMI-Dynamic Model GTEx IDs

**If your biobank is already mapped to build hg38**:

- You can construct GTEx IDs directly using the format: `chr{CHR}_{POS}_{REF}_{ALT}_b38`
- Example: `chr1_752566_A_G_b38`

**If your biobank uses hg19/GRCh37 or rsIDs**:

- Use the provided mapping file: `models/mapping/{tissue}_model_variants.txt.gz`
- This file maps model variants to rsIDs (file size: ~6 MB)

**Mapping file format**:

```
CHR  pos_37   ref  alt  variant_id_hg38              rs_id        dbSNP151_GRCh38p7  tissue
1    14677    G    A    chr1_14677_G_A_b38          rs2123123    ...                Adipose-Subcutaneous
1    256498   A    G    chr1_256498_A_G_b38         rs369556846  ...                Adipose-Subcutaneous
```

**Columns**:

- `CHR`: Chromosome number
- `pos_37`: Position in hg19/GRCh37 (if you need to lift over)
- `ref`: Reference allele
- `alt`: Alternate allele (this is what should be counted!)
- `variant_id_hg38`: GTEx ID in hg38 format (this is what you need as column names!)
- `rs_id`: rsID from dbSNP
- `dbSNP151_GRCh38p7`: dbSNP build information
- `tissue`: Which tissue(s) this SNP is used in

**Quality control recommendation**:

- We applied an **INFO > 0.8** filter to ensure high-quality imputed variants
- You should apply appropriate QC filters based on your biobank's data quality
- Save the mapping file for referemce and filtered SNP list in a text file (no header) to use for plink2 extraction in the next step

#### 1.2 Ensure the Alternate Alleles Are Counted

**Critical**: The models expect counts of the **alternate allele**, not the reference allele.

To ensure the alternate alleles are counted correctly with plink2:

1. Create a file for the `--export-allele` flag
2. Format: Two columns with **no column names**
   - Column 1: SNP ID (as used in your biobank files)
   - Column 2: Alternate allele (from the mapping file)

**Example allele file**:

```
rs2123123      A
rs369556846    G
rs111497954    A
```

This ensures plink2 counts the correct allele orientation.

#### 1.3 Subset and Reformat the Dosage Files Using plink2

Extract the model SNPs from your biobank genotype files and convert to the required format:

```bash
plink2 --bfile your_biobank \
       --extract snp_list.txt \
       --recode A \
       --export-allele allele_file.txt \
       --rm-dup force-first \
       --out biobank_model_snps
```

**Required plink2 flags**:

- `--extract`: Extract only the SNPs present in the models
- `--recode A`: Convert to additive (0/1/2) dosage format
- `--export-allele`: Specify which allele to count (ensures alternate allele)
- `--rm-dup force-first`: Handle any duplicate SNP IDs

**Output**: `biobank_model_snps.raw`

#### 1.4 Change Column Names to GTEx Format

In the final step, you will need to:

1. Change the .raw file column names from your biobank IDs to GTEx `chr{N}_{position}_{ref}_{alt}_b38` format using the mapping from step 1.1 to rename columns
2. **Make sure there is still a column named `IID`** - this will be used by the prediction script

**Final format example**:

```
IID           chr1_14677_G_A_b38    chr1_256498_A_G_b38    chr1_286478_A_C_b38
SAMPLE001     0.0                   1.2                    2.0
SAMPLE002     1.0                   0.0                    1.8
SAMPLE003     2.0                   1.0                    0.2
```

**Notes**:

- Tab-separated format
- Can be gzipped (`.raw.gz`)
- First column **must** be named `IID`
- All other columns are SNPs with GTEx variant IDs
- Values are dosages (0.0 to 2.0)

------

## Step 2: Prepare the BMI File

The last file we need to predict expression is a file with BMI values for each individual.

### Required Format

**Two columns with header**:

- `SUBJID`: Individual identifiers (must match `IID` from genotype file)
- `BMI`: Body Mass Index values (numeric)

**Example**:

```
SUBJID        BMI
SAMPLE001     25.3
SAMPLE002     31.2
SAMPLE003     24.8
SAMPLE004     27.6
```

**Notes**:

- Tab-separated text file
- Header row is required
- `SUBJID` must exactly match `IID` in your genotype file
- BMI should be numeric (typical range 15-50)
- File can be gzipped

------

## Step 3: Predict Gene Expression

Now use the prediction script to generate predicted expression values.

### Basic Usage

```bash
Rscript apply_models/predict_expression.R \
    --tissue Adipose-Subcutaneous \
    --chr 1 \
    --part 1 \
    --genotype-file biobank_chr1.raw.gz \
    --bmi-file biobank_bmi.txt \
    --output-dir results/
```

### Required Arguments

| Argument          | Description                                 | Example                |
| ----------------- | ------------------------------------------- | ---------------------- |
| `--tissue`        | Tissue name                                 | `Adipose-Subcutaneous` |
| `--chr`           | Chromosome number (1-22)                    | `1`                    |
| `--part`          | Part number (1-10)                          | `1`                    |
| `--genotype-file` | Path to formatted genotype file from Step 1 | `biobank_chr1.raw.gz`  |
| `--bmi-file`      | Path to BMI file from Step 2                | `biobank_bmi.txt`      |
| `--output-dir`    | Directory where output files will be saved  | `results/`             |

### Optional Arguments

| Argument      | Description                                 | Default   | Example      |
| ------------- | ------------------------------------------- | --------- | ------------ |
| `--model-dir` | Path to directory containing models         | `models/` | `../models/` |
| `--n-parts`   | Total number of parts (for splitting genes) | `10`      | `10`         |

### Available Tissues

Pre-built models are available for the following GTEx tissues:

- Adipose-Subcutaneous
- Brain-Cerebellum
- Brain-Cortex
- Brain-NucleusAccumbens
- Colon-Transverse
- Esophagus-Mucosa
- Esophagus-Muscularis
- Nerve-Tibial

### Parts System

**Why parts?** For memory management in large cohorts, genes on each chromosome are split into 10 equal parts.

- Each part processes ~10% of genes on that chromosome
- Parts can be run in parallel on a cluster
- Total jobs per tissue: 22 chromosomes × 10 parts = 220 jobs

### Processing All Chromosomes and Parts

To predict expression across all chromosomes for a tissue:

```bash
# Set your parameters
TISSUE="Adipose-Subcutaneous"
BMI_FILE="/path/to/biobank_bmi.txt"
GENO_DIR="/path/to/genotype_files"
OUTPUT_DIR="/path/to/results"

# Loop through all chromosomes and parts
for CHR in {1..22}; do
    for PART in {1..10}; do
        echo "Processing chromosome ${CHR}, part ${PART}..."
        
        Rscript apply_models/predict_expression.R \
            --tissue ${TISSUE} \
            --chr ${CHR} \
            --part ${PART} \
            --genotype-file ${GENO_DIR}/biobank_chr${CHR}.raw.gz \
            --bmi-file ${BMI_FILE} \
            --output-dir ${OUTPUT_DIR}
    done
done

echo "All chromosomes and parts completed for ${TISSUE}"
```

### Output Files

For each chromosome and part, the script generates **two output files**:

#### 1. Predicted Expression: `{tissue}_{chr}_{part}_BGREX.txt.gz`

**BGREX** = **B**MI-aware **G**enetically **R**egulated **Ex**pression*

*A reminder that while all the genes can select BMI as a predictor, not all of them do. You can see the predictors selected in the beta files in the model directory 

**Format**:

```
SUBJID         ENSG00000000003.15    ENSG00000000005.6    ENSG00000000419.12
SAMPLE001      -0.234                0.891                1.234
SAMPLE002      0.567                 -1.234               -0.456
SAMPLE003      1.234                 0.234                0.678
```

- **Rows**: Individuals (SUBJID)
- **Columns**: Genes (Ensembl gene IDs with version numbers)
- **Values**: Predicted expression levels
  - Can be positive or negative (these are residualized predictions)
  - `NA` indicates the gene could not be predicted (no available SNPs)

#### 2. Imputation Statistics: `{tissue}_{chr}_{part}_imputestats.txt.gz`

**Format**:

```
gene                 N_pred_all    N_pred_used    pred_used
ENSG00000000003.15   5             5              chr1_752566_A_G_b38,chr1_753405_G_A_b38,...
ENSG00000000005.6    3             2              chr1_11234_C_T_b38,chr1_11567_A_G_b38
ENSG00000000419.12   8             0              NA
```

- `gene`: Ensembl gene ID
- `N_pred_all`: Total number of predictors in the original model (SNPs + BMI + interactions)
- `N_pred_used`: Number of predictors available in your biobank data
- `pred_used`: Comma-separated list of predictors actually used
  - Includes SNP IDs, "BMI", and interaction terms (e.g., "chr1_12345_A_G_b38*BMI")
  - `NA` if no predictors available

**Use this file to**:

- Check prediction quality (higher N_pred_used = better)
- Identify genes with poor prediction (N_pred_used = 0)
- Filter genes for downstream analysis based on number of predictors

------

## Step 4: Test Associations with Phenotypes

After predicting expression, test for associations between BMI-aware predicted expression (BGREX) and your phenotype of interest.

### Prepare Input Files

You need three files for association testing:

#### 4.1 Predicted Expression Files

These are the outputs from Step 3:

- Files: `{tissue}_{chr}_{part}_BGREX.txt.gz`
- Already prepared!

#### 4.2 Phenotype File

**Required format**:

- **Columns**: `SUBJID`, `[your_phenotype_name]`
- **SUBJID**: Must match `SUBJID` from predicted expression files
- **Phenotype column**: Name it whatever you want (you'll specify with `--phenotype-name`)

**For binary phenotypes (case/control)**:

```
SUBJID          disease_status
SAMPLE001       Case
SAMPLE002       Control
SAMPLE003       Case
SAMPLE004       Control
```

**For quantitative phenotypes**:

```
SUBJID          LDL_cholesterol
SAMPLE001       120.5
SAMPLE002       145.2
SAMPLE003       98.3
SAMPLE004       132.7
```

**Notes**:

- Tab-separated text file
- Header row required
- Binary values can be: `Case`/`Control`, `1`/`0`, or any two distinct values
- File can be gzipped

#### 4.3 Covariate File (Recommended)

**Required format**:

- **Columns**: `SUBJID`, `Age`, `Sex`, `PC1`, `PC2`, ..., `PC10`, (optionally `BMI`)
- **SUBJID**: Must match `SUBJID` from predicted expression files

**Example**:

```
SUBJID          Age     Sex     PC1         PC2         PC3         PC4         PC5         BMI
SAMPLE001       45      Male    -0.001234   0.005678    -0.002345   0.001234    0.000567    25.3
SAMPLE002       52      Female  0.002345    -0.003456   0.001234    -0.000789   0.002345    31.2
SAMPLE003       38      Male    -0.003456   0.001234    0.000567    0.002345    -0.001234   24.8
```

**Standard covariates to include**:

- **Age**: Age at assessment
- **Sex**: Male/Female or 0/1
- **PC1-PC10**: First 10 genetic principal components (for ancestry adjustment)
- **BMI**: Body Mass Index (optional - use `--test-with-bmi` flag to include)

**Notes**:

- Tab-separated text file
- Header row required
- You choose which covariates to include in the model
- File can be gzipped
- **All files use `SUBJID` for consistency**

### Basic Usage

```bash
Rscript apply_models/test_association.R \
    --tissue Adipose-Subcutaneous \
    --chr 1 \
    --part 1 \
    --expression-dir results/ \
    --phenotype-file phenotypes.txt \
    --covariate-file covariates.txt \
    --output-dir associations/ \
    --phenotype-type binary \
    --phenotype-name disease_status \
    --covariates 'Age,Sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10'
```

### Required Arguments

| Argument           | Description                                     | Example                |
| ------------------ | ----------------------------------------------- | ---------------------- |
| `--tissue`         | Tissue name                                     | `Adipose-Subcutaneous` |
| `--chr`            | Chromosome number (1-22)                        | `1`                    |
| `--part`           | Part number (1-10)                              | `1`                    |
| `--expression-dir` | Directory containing predicted expression files | `results/`             |
| `--phenotype-file` | Phenotype file (with `SUBJID` column)           | `phenotypes.txt`       |
| `--output-dir`     | Output directory for results                    | `associations/`        |
| `--phenotype-type` | `binary` or `quantitative`                      | `binary`               |
| `--phenotype-name` | Name of phenotype column                        | `disease_status`       |

### Optional Arguments

| Argument           | Description                           | Default | Example                 |
| ------------------ | ------------------------------------- | ------- | ----------------------- |
| `--covariate-file` | Covariate file (with `SUBJID` column) | None    | `covariates.txt`        |
| `--covariates`     | Comma-separated covariate names       | None    | `'Age,Sex,PC1,PC2,PC3'` |
| `--test-with-bmi`  | Include BMI as covariate              | FALSE   | `TRUE`                  |

### Testing With and Without BMI as a Covariate

You can test associations with or without BMI as a covariate. This is important because:

- **Without BMI**: Tests if BGREX (which includes BMI effects) associates with phenotype
- **With BMI**: Tests if BGREX associates with phenotype *independent* of BMI

Run the script twice to get both results:

```bash
# Test WITHOUT BMI as covariate
Rscript apply_models/test_association.R \
    --tissue Adipose-Subcutaneous \
    --chr 1 \
    --part 1 \
    --expression-dir results/ \
    --phenotype-file phenotypes.txt \
    --covariate-file covariates.txt \
    --output-dir associations/noBMI/ \
    --phenotype-type binary \
    --phenotype-name disease_status \
    --covariates 'Age,Sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10'

# Test WITH BMI as covariate
Rscript apply_models/test_association.R \
    --tissue Adipose-Subcutaneous \
    --chr 1 \
    --part 1 \
    --expression-dir results/ \
    --phenotype-file phenotypes.txt \
    --covariate-file covariates.txt \
    --output-dir associations/withBMI/ \
    --phenotype-type binary \
    --phenotype-name disease_status \
    --covariates 'Age,Sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10' \
    --test-with-bmi TRUE
```

### Processing All Chromosomes and Parts

To test associations for all chromosomes and parts:

```bash
TISSUE="Adipose-Subcutaneous"
EXPR_DIR="results/"
PHENO_FILE="phenotypes.txt"
COV_FILE="covariates.txt"
PHENO_NAME="disease_status"
COVARIATES="Age,Sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
OUTPUT_DIR="associations/"

for CHR in {1..22}; do
    for PART in {1..10}; do
        echo "Testing associations: chr ${CHR}, part ${PART}..."
        
        Rscript apply_models/test_association.R \
            --tissue ${TISSUE} \
            --chr ${CHR} \
            --part ${PART} \
            --expression-dir ${EXPR_DIR} \
            --phenotype-file ${PHENO_FILE} \
            --covariate-file ${COV_FILE} \
            --output-dir ${OUTPUT_DIR} \
            --phenotype-type binary \
            --phenotype-name ${PHENO_NAME} \
            --covariates "${COVARIATES}"
    done
done
```

### Output File Format

The association testing script produces per-part results:

**Filename**: `{tissue}_{chr}_{part}_associations.txt.gz`

```
Gene         Beta       SE         P           N
ENSG...15    0.123      0.045      0.0063      5000
ENSG...06   -0.234      0.067      0.0005      5000
ENSG...12    0.045      0.089      0.6134      5000
```

**Columns**:

- `Gene`: Ensembl gene ID
- `Beta`: Effect size (log odds ratio for binary, beta for quantitative)
- `SE`: Standard error
- `P`: P-value
- `N`: Sample size (individuals with non-NA BGREX)

**Note**: FDR is NOT included in per-part results. Genome-wide FDR is calculated in Step 5.

### Ancestry Stratification

For multi-ancestry cohorts, we recommend stratifying by ancestry:

**Why stratify?**

- Controls for population stratification
- Allows ancestry-specific effect sizes
- More powerful than simply adjusting for PCs in mixed-ancestry samples

**Workflow**:

1. **Define ancestry groups** using genetic PCs, self-reported ancestry, or reference panels
2. **Predict expression separately** for each ancestry (if needed) OR use same predictions
3. **Test associations separately** for each ancestry group
4. **Meta-analyze** results across ancestries using METAL or similar tools

**Note**: In our original study, we used genetically-defined ancestry groups based on PC clustering and meta-analyzed results with METAL.

------

## Step 5: Assign Genome-Wide Significance

**CRITICAL**: You must combine results across all chromosomes and parts to calculate genome-wide FDR. Per-chromosome or per-part p-values do not account for multiple testing across the entire genome.

### Basic Usage

```bash
Rscript apply_models/assign_significance.R \
    --tissue Adipose-Subcutaneous \
    --input-dir associations/ \
    --output-file associations/Adipose-Subcutaneous_genome_wide_FDR.txt.gz
```

### Required Arguments

| Argument        | Description                               | Example                |
| --------------- | ----------------------------------------- | ---------------------- |
| `--tissue`      | Tissue name                               | `Adipose-Subcutaneous` |
| `--input-dir`   | Directory containing association results  | `associations/`        |
| `--output-file` | Output file for combined results with FDR | `final_results.txt.gz` |

### Optional Arguments

| Argument    | Description                      | Default |
| ----------- | -------------------------------- | ------- |
| `--n-chr`   | Number of chromosomes to combine | `22`    |
| `--n-parts` | Number of parts per chromosome   | `10`    |

### What This Script Does

1. **Combines** all association results across chromosomes and parts
2. **Calculates** genome-wide FDR using Benjamini-Hochberg method
3. **Saves** combined results with FDR column
4. **Exports** top hits and significant genes

### Output Files

**Main output**: `{tissue}_genome_wide_FDR.txt.gz`

```
Gene         Beta    SE      P        N     Chr  Part  FDR
ENSG...15    0.123   0.045   0.0001   5000  1    1     0.023
ENSG...06   -0.234   0.067   0.0005   5000  2    3     0.045
```

**Additional outputs**:

- `{tissue}_top50_byP.txt` - Top 50 genes ranked by p-value
- `{tissue}_significant_FDR0.05.txt` - All genes with FDR < 0.05

### Why Genome-Wide FDR Matters

**Example**: If you test 20,000 genes across 22 chromosomes × 10 parts:

- Total tests: ~220 separate jobs, but ~20,000 unique genes
- FDR calculated on just one part (~100 genes) would be **too liberal**
- FDR calculated on just one chromosome (~1,000 genes) would be **too liberal**
- **Always calculate FDR on the complete set of genome-wide tests**

### Summary Statistics

The script outputs comprehensive summary statistics:

```
================================================================================
 Genome-Wide Association Summary
================================================================================
  Total genes tested:              18,523
  Genes with valid P-values:       17,891 (96.6%)
  Genes with missing P-values:     632 (3.4%)

  Significant at P < 0.05:         1,245 (7.0%)
  Significant at P < 0.01:         423 (2.4%)
  Significant at P < 0.001:        87 (0.5%)

  Significant at FDR < 0.05:       256 (1.4%)
  Significant at FDR < 0.10:       489 (2.7%)
  Significant at FDR < 0.20:       892 (5.0%)

  Minimum P-value:                 3.24e-12
  Minimum FDR:                     0.0001
================================================================================
```

------

## Complete Workflow Example

Here's a complete end-to-end example for testing association with a disease:

```bash
#!/bin/bash
# Complete workflow: Predict expression and test associations

# ============================================================================
# SETUP
# ============================================================================
TISSUE="Adipose-Subcutaneous"
BIOBANK_DIR="/data/my_biobank"
WORK_DIR="/analysis/bmi_predictions"
PHENO_FILE="${WORK_DIR}/disease_phenotypes.txt"
COV_FILE="${WORK_DIR}/covariates.txt"
PHENO_NAME="disease_status"
PHENO_TYPE="binary"
COVARIATES="Age,Sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"

# ============================================================================
# STEP 3: PREDICT EXPRESSION (220 jobs - can run in parallel)
# ============================================================================
echo "Step 3: Predicting expression..."
mkdir -p ${WORK_DIR}/results

for CHR in {1..22}; do
    for PART in {1..10}; do
        Rscript apply_models/predict_expression.R \
            --tissue ${TISSUE} \
            --chr ${CHR} \
            --part ${PART} \
            --genotype-file ${BIOBANK_DIR}/biobank_chr${CHR}.raw.gz \
            --bmi-file ${WORK_DIR}/biobank_bmi.txt \
            --output-dir ${WORK_DIR}/results/
    done
done

echo "Step 3 complete: Expression predicted"

# ============================================================================
# STEP 4: TEST ASSOCIATIONS (220 jobs - can run in parallel)
# ============================================================================
echo "Step 4: Testing associations..."
mkdir -p ${WORK_DIR}/associations

for CHR in {1..22}; do
    for PART in {1..10}; do
        Rscript apply_models/test_association.R \
            --tissue ${TISSUE} \
            --chr ${CHR} \
            --part ${PART} \
            --expression-dir ${WORK_DIR}/results/ \
            --phenotype-file ${PHENO_FILE} \
            --covariate-file ${COV_FILE} \
            --output-dir ${WORK_DIR}/associations/ \
            --phenotype-type ${PHENO_TYPE} \
            --phenotype-name ${PHENO_NAME} \
            --covariates "${COVARIATES}"
    done
done

echo "Step 4 complete: Associations tested"

# ============================================================================
# STEP 5: ASSIGN GENOME-WIDE SIGNIFICANCE (1 job)
# ============================================================================
echo "Step 5: Calculating genome-wide FDR..."

Rscript apply_models/assign_significance.R \
    --tissue ${TISSUE} \
    --input-dir ${WORK_DIR}/associations/ \
    --output-file ${WORK_DIR}/associations/${TISSUE}_genome_wide_FDR.txt.gz

echo "Step 5 complete: FDR assigned"

# ============================================================================
# RESULTS
# ============================================================================
echo "========================================="
echo "Workflow complete!"
echo "========================================="
echo "Results:"
echo "  Full results:       ${WORK_DIR}/associations/${TISSUE}_genome_wide_FDR.txt.gz"
echo "  Top 50 hits:        ${WORK_DIR}/associations/${TISSUE}_top50_byP.txt"
echo "  Significant genes:  ${WORK_DIR}/associations/${TISSUE}_significant_FDR0.05.txt"
echo "========================================="
```

------

## Citation

If you use these BMI-dynamic models in your research, please cite:

```
[Citation information will be added upon publication]
```

------

## Need Help?

- Review the [FAQ](https://claude.ai/docs/FAQ.md)
- Check existing [GitHub Issues](https://github.com/rsigner/BMI_dynamic_models/issues)
- Open a new issue with:
  - Your error message
  - Description of your biobank data format
  - Steps you've already tried
# Building BMI-Dynamic Gene Expression Predictor Models with GTEx

## Overview

This workflow allows you to **train new BMI-dynamic gene expression prediction models** from paired DNA and RNA data. The models identify genes whose expression is predicted by:
- SNP main effects
- BMI main effect
- SNP × BMI interaction effects

**Original training**: These methods were used to build models from GTEx v8 data (~700 individuals, 8 tissues).

---

## Pipeline Overview

```
Your Paired DNA + RNA Data
         ↓
Step 1: Create Cross-Validation Folds
         ↓
Step 2: Build Models (by chromosome and part)
         ↓
Step 3: Calculate Holdout Statistics
         ↓
     New Models!
```

**Why parts?** Training glinternet models is memory-intensive. Each chromosome is split into 30 parts for parallel processing.

---

## Step 1: Create Cross-Validation Folds

Create stratified k-fold splits based on BMI categories to ensure balanced representation across folds.

### Prepare Input File

You need a phenotype file with:
- **Columns**: `SUBJID`, `BMI`
- **SUBJID**: Subject identifiers for individuals in your study
- **BMI**: Body Mass Index values

**Example**:
```
SUBJID          BMI
GTEX-1117F      24.5
GTEX-111CU      31.2
GTEX-111FC      27.8
```

**Notes**:
- Tab-separated text file
- Header row required
- File can be gzipped
- Only include samples that passed QC

### Basic Usage

```bash
python build_models/create_folds.py \
    --tissue Brain-Cortex \
    --phenotype-file data/Brain-Cortex_phenoData.txt.gz \
    --output-dir output/folds/
```

### Arguments

| Argument | Description | Default | Example |
|----------|-------------|---------|---------|
| `--tissue` | Tissue name | Required | `Brain-Cortex` |
| `--phenotype-file` | Phenotype file path | Required | `phenotypes.txt.gz` |
| `--output-dir` | Output directory | Required | `output/folds/` |
| `--n-folds` | Number of CV folds | `4` | `4` |
| `--random-state` | Random seed | `42` | `42` |

### BMI Categories

Samples are stratified into four categories based on CDC definitions:
- **Underweight**: BMI < 18.5
- **Normal Weight**: 18.5 ≤ BMI < 25
- **Overweight**: 25 ≤ BMI < 30
- **Obese**: BMI ≥ 30

### Output File

**Filename**: `{tissue}_CV_test_ID_map_weightcategory.txt.gz`

```
SUBJID          BMI     Weight_Category    Fold_in_test
GTEX-1117F      24.5    Normal Weight      1
GTEX-111CU      31.2    Obese              2
GTEX-111FC      27.8    Overweight         3
```

**Columns**:
- `SUBJID`: Subject identifier
- `BMI`: BMI value
- `Weight_Category`: BMI category
- `Fold_in_test`: Fold assignment (1 to 4)

**Note**: The script uses 4 folds. Folds 1-3 are used for training (with internal 3-fold CV), and fold 4 is held out for validation.

---

## Step 2: Build Prediction Models

Train glinternet models for each gene using genotype data, gene expression, and BMI.

### Prepare Input Files

You need six types of input files:

#### 2.1 Gene List File

A text file listing genes to process (one Ensembl ID per line):

```
ENSG00000000003.15
ENSG00000000005.6
ENSG00000000419.12
```

**Notes**:
- One gene per line
- Use full Ensembl IDs with version numbers
- Typically create separate lists per tissue and chromosome
- Only include genes expressed in the tissue of interest

#### 2.2 Fold Assignment File

This is the output from Step 1:
- File: `{tissue}_CV_test_ID_map_weightcategory.txt.gz`
- Already prepared!

#### 2.3 Phenotype File

Same file used in Step 1:
- **Columns**: `SUBJID`, `BMI`
- Tab-separated, can be gzipped

#### 2.4 Expression Residuals File

**Required format**:
- **Rows**: Samples (individuals)
- **Columns**: `ID` (sample IDs), then one column per gene (Ensembl IDs)
- **Values**: Expression residuals (corrected for covariates and/or latent factor variables)

**Important preprocessing**:

- Remove unwanted variation from expression data BEFORE model building
- Keep autosomes only (remove sex chromosomes)
- We corrected GTEx expression for:
  - Surrogate variables (SVA) with "be" method
  - Genetic principal components
- Expression values should be normalized and residualized
- Tab-separated, can be gzipped

**Example**:
```
ID              ENSG00000000003.15    ENSG00000000005.6    ENSG00000000419.12
GTEX-1117F      -0.234                0.456                1.123
GTEX-111CU      0.567                 -0.234               -0.456
GTEX-111FC      1.234                 0.123                0.789
```

#### 2.5 Genotype Dosage File

**Required format**:
- **Rows**: Samples (individuals)
- **Columns**: `ID` (sample IDs), then one column per SNP
- **Column names**: GTEx variant IDs (format: `chr{N}_{position}_{ref}_{alt}_b38`)
- **Values**: Dosages (0.0 to 2.0) or hardcall 0/1/2 (counting alternate alleles)

**Quality control filters (apply BEFORE model building)**:

- **MAF > 0.05**: Remove rare variants
- **No ambiguous SNPs**: Remove A/T, T/A, G/C, C/G SNPs
- **Missingness**: Remove SNPs with any missing values
- **HWE**: Standard Hardy-Weinberg filters
- Create separate files per chromosome

**Notes**:

- Tab-separated, can be gzipped
- One file per chromosome

**Example**:

```
ID              chr1_752566_A_G_b38    chr1_753405_G_A_b38    chr1_754182_A_G_b38
GTEX-1117F      0.0                    1.2                    2.0
GTEX-111CU      1.0                    0.0                    1.8
GTEX-111FC      2.0                    1.0                    0.2
```

#### 2.6 Gene Feature File

**Required format**:
- Gene annotations with coordinates

**Example**:

```
Name                    Chr    start        end          gene_type        Description
ENSG00000000419.12      20     50934867     50958555     protein_coding   DPM1
ENSG00000000457.13      1      169849631    169894267    protein_coding   SCYL3
ENSG00000000460.16      1      169662007    169854080    protein_coding   C1orf112
```

Required columns**:
- `Name`: Ensembl gene ID (must match expression file columns)
- `Chr`: Chromosome (1-22)
- `start`: Gene start position (for defining cis-window)
- `end`: Gene end position
- `gene_type`: Gene biotype (protein_coding, lncRNA, etc.)
- `Description`: Gene symbol or description (optional)

### Basic Usage

```bash
Rscript build_models/build_glinternet_model.R \
    --tissue Brain-Cortex \
    --chr 1 \
    --part 1 \
    --gene-list-file genes/Brain-Cortex_chr1.txt \
    --fold-file folds/Brain-Cortex_CV_test_ID_map_weightcategory.txt.gz \
    --phenotype-file data/Brain-Cortex_phenoData.txt.gz \
    --expression-file data/Brain-Cortex_Residuals.txt.gz \
    --genotype-file data/chr1_WGS.txt.gz \
    --feature-file data/Brain-Cortex_featureData.txt.gz \
    --output-dir output/models/
```

### Required Arguments

| Argument | Description | Example |
|----------|-------------|---------|
| `--tissue` | Tissue name | `Brain-Cortex` |
| `--chr` | Chromosome (1-22) | `1` |
| `--part` | Part number (1-30) | `1` |
| `--gene-list-file` | Gene list file | `genes_chr1.txt` |
| `--fold-file` | CV fold assignments | `folds.txt.gz` |
| `--phenotype-file` | Phenotype file | `phenotypes.txt.gz` |
| `--expression-file` | Expression residuals | `expression.txt.gz` |
| `--genotype-file` | Genotype dosages for this chr | `chr1_WGS.txt.gz` |
| `--feature-file` | Gene features/annotations | `features.txt.gz` |
| `--output-dir` | Output directory | `output/models/` |

### Optional Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--n-parts` | Total number of parts per chromosome | `30` |
| `--cis-window` | Cis-window size in bp | `1000000` (1 Mb) |
| `--rsq-threshold` | Minimum R² to save model | `0.01` |
| `--n-folds-cv` | Internal CV folds for glinternet | `3` |
| `--random-seed` | Random seed | `1702` |

### Parts System

**Why 30 parts?** Training glinternet models on genes with thousands of SNPs is very memory-intensive.

- Each chromosome is split into 30 equal parts
- Parts can be run in parallel on a cluster
- Total jobs per tissue: 22 chromosomes × 30 parts = 660 jobs

### How the Model Building Works

For each gene:

1. **Define cis-window**: ±1 Mb from gene start position
2. **Extract SNPs**: All SNPs within the cis-window
3. **Train model**: Use glinternet with 3-fold internal CV
   - Predictors: SNPs + BMI
   - Interactions: SNP × BMI
   - Regularization: Automatic lambda selection
4. **Validate**: Predict on held-out 4th fold
5. **Save**: If training R² > 0.01, save model

### Processing All Chromosomes and Parts

To build models across all chromosomes:

```bash
TISSUE="Brain-Cortex"
FOLD_FILE="folds/${TISSUE}_CV_test_ID_map_weightcategory.txt.gz"
PHENO_FILE="data/${TISSUE}_phenoData.txt.gz"
EXPR_FILE="data/${TISSUE}_Residuals.txt.gz"
FEATURE_FILE="data/${TISSUE}_featureData.txt.gz"
OUTPUT_DIR="output/models/${TISSUE}/"

for CHR in {1..22}; do
    GENE_LIST="genes/${TISSUE}_chr${CHR}.txt"
    GENO_FILE="data/chr${CHR}_WGS.txt.gz"
    
    for PART in {1..30}; do
        echo "Processing chromosome ${CHR}, part ${PART}..."
        
        Rscript build_models/build_glinternet_model.R \
            --tissue ${TISSUE} \
            --chr ${CHR} \
            --part ${PART} \
            --gene-list-file ${GENE_LIST} \
            --fold-file ${FOLD_FILE} \
            --phenotype-file ${PHENO_FILE} \
            --expression-file ${EXPR_FILE} \
            --genotype-file ${GENO_FILE} \
            --feature-file ${FEATURE_FILE} \
            --output-dir ${OUTPUT_DIR}
    done
done
```

### Output Files

For each gene, three output files are created:

#### 1. Model Summary: `{tissue}_{gene}_model_summary.txt.gz`

```
Gene                Tissue        Rsq     p_Rsq    CV_Lambda  CV_MStEr  CV_MSE   R2_rep   p_rep    n_predictors
ENSG00000000003.15  Brain-Cortex  0.123   0.0001   0.0234     0.045     0.123    0.098    0.002    12
```

**Columns**:
- `Gene`: Ensembl gene ID
- `Tissue`: Tissue name
- `Rsq`: Training R² (3-fold CV)
- `p_Rsq`: P-value for training R²
- `CV_Lambda`: Best lambda from cross-validation
- `CV_MStEr`: CV mean standard error
- `CV_MSE`: CV mean squared error
- `R2_rep`: Validation R² (held-out fold 4)
- `p_rep`: P-value for validation R²
- `n_predictors`: Number of predictors in final model

#### 2. Model Coefficients: `{tissue}_{gene}_betas.txt.gz`

```
name                   Beta       predictor_type       Predictor
BMI                    0.123      Main Effect          BMI
chr1_752566_A_G_b38    0.456      Main Effect          chr1_752566_A_G_b38
chr1_753405_G_A_b38    0.234      Interaction Effect   chr1_753405_G_A_b38*BMI
```

**Columns**:
- `name`: Predictor name (SNP ID or "BMI")
- `Beta`: Coefficient value
- `predictor_type`: "Main Effect" or "Interaction Effect"
- `Predictor`: Full predictor name (includes "*BMI" for interactions)

#### 3. Model Fit Object: `{tissue}_{gene}_fit.rds`

R object containing the full glinternet model fit. Can be loaded for predictions:

```r
fit <- readRDS("Brain-Cortex_ENSG00000000003.15_fit.rds")
predictions <- predict(fit, new_data)
```

---

## Step 3: Calculate Holdout Statistics and Combine Results

After building models for all genes, combine results, calculate validation statistics, and create plots.

### Basic Usage

```bash
Rscript build_models/holdout_statistics.R \
    --tissue Brain-Cortex \
    --gene-list-file genes/Brain-Cortex_all_genes.txt \
    --model-dir output/models/Brain-Cortex/ \
    --feature-file data/Brain-Cortex_featureData.txt.gz \
    --output-dir output/results/ \
    --plot-dir output/plots/
```

### Required Arguments

| Argument | Description | Example |
|----------|-------------|---------|
| `--tissue` | Tissue name | `Brain-Cortex` |
| `--gene-list-file` | List of all genes tested | `all_genes.txt` |
| `--model-dir` | Directory with model files | `models/` |
| `--feature-file` | Gene feature file | `features.txt.gz` |
| `--output-dir` | Output directory for results | `results/` |

### Optional Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--plot-dir` | Plot output directory | Same as `output-dir` |
| `--rsq-threshold` | R² threshold for significance | `0.01` |
| `--pval-threshold` | P-value threshold | `0.05` |

### What This Script Does

1. **Combines** all model summaries across genes
2. **Filters** for significant models (R² > 0.01, P < 0.05 in training AND validation)
3. **Categorizes** models by predictor types:
   - Models with SNP × BMI interactions
   - Models with BMI main effect only
   - Models with SNP main effects only
4. **Creates** validation plot (training R² vs validation R²)
5. **Exports** combined results

### Output Files

#### 1. Combined Model Summaries: `{tissue}_r2_0.01_p_0.05_models.txt.gz`

All significant models with annotations:
- Gene information
- Model statistics
- Predictor type category

#### 2. All Coefficients: `{tissue}_r2_0.01_p_0.05_betas.txt.gz`

All coefficients for all significant models:
- SNP IDs
- Beta values
- Predictor types
- Gene annotations

#### 3. Validation Plot: `{tissue}_R2_holdout_0.01_0.05.pdf`

Scatter plot showing:
- X-axis: Training R² (3-fold CV)
- Y-axis: Validation R² (held-out fold)
- Colors: Predictor types
- Diagonal line: Perfect concordance
- Regression line: Actual relationship

---

## Complete Workflow Example

Here's a complete end-to-end example for building models for one tissue:

```bash
#!/bin/bash
# Complete model building workflow

TISSUE="Brain-Cortex"
DATA_DIR="data"
OUTPUT_BASE="output"

# ============================================================================
# STEP 1: CREATE CV FOLDS
# ============================================================================
echo "Step 1: Creating cross-validation folds..."

python build_models/create_folds.py \
    --tissue ${TISSUE} \
    --phenotype-file ${DATA_DIR}/${TISSUE}_phenoData.txt.gz \
    --output-dir ${OUTPUT_BASE}/folds/ \
    --n-folds 4

echo "Step 1 complete"

# ============================================================================
# STEP 2: BUILD MODELS (660 jobs - run in parallel on cluster)
# ============================================================================
echo "Step 2: Building models..."

FOLD_FILE="${OUTPUT_BASE}/folds/${TISSUE}_CV_test_ID_map_weightcategory.txt.gz"
PHENO_FILE="${DATA_DIR}/${TISSUE}_phenoData.txt.gz"
EXPR_FILE="${DATA_DIR}/${TISSUE}_Residuals.txt.gz"
FEATURE_FILE="${DATA_DIR}/${TISSUE}_featureData.txt.gz"
MODEL_DIR="${OUTPUT_BASE}/models/${TISSUE}/"

mkdir -p ${MODEL_DIR}

for CHR in {1..22}; do
    GENE_LIST="${DATA_DIR}/gene_lists/${TISSUE}_chr${CHR}.txt"
    GENO_FILE="${DATA_DIR}/genotypes/chr${CHR}_WGS.txt.gz"
    
    for PART in {1..30}; do
        # Submit to cluster or run locally
        Rscript build_models/build_glinternet_model.R \
            --tissue ${TISSUE} \
            --chr ${CHR} \
            --part ${PART} \
            --gene-list-file ${GENE_LIST} \
            --fold-file ${FOLD_FILE} \
            --phenotype-file ${PHENO_FILE} \
            --expression-file ${EXPR_FILE} \
            --genotype-file ${GENO_FILE} \
            --feature-file ${FEATURE_FILE} \
            --output-dir ${MODEL_DIR} &
    done
    
    # Wait for all parts of this chromosome to finish
    wait
done

echo "Step 2 complete"

# ============================================================================
# STEP 3: CALCULATE HOLDOUT STATISTICS
# ============================================================================
echo "Step 3: Calculating holdout statistics and combining results..."

# Create gene list with all genes tested
cat ${DATA_DIR}/gene_lists/${TISSUE}_chr*.txt > ${OUTPUT_BASE}/${TISSUE}_all_genes.txt

Rscript build_models/holdout_statistics.R \
    --tissue ${TISSUE} \
    --gene-list-file ${OUTPUT_BASE}/${TISSUE}_all_genes.txt \
    --model-dir ${MODEL_DIR} \
    --feature-file ${FEATURE_FILE} \
    --output-dir ${OUTPUT_BASE}/results/ \
    --plot-dir ${OUTPUT_BASE}/plots/

echo "Step 3 complete"

# ============================================================================
# RESULTS
# ============================================================================
echo "========================================="
echo "Workflow complete!"
echo "========================================="
echo "Results:"
echo "  Models:        ${OUTPUT_BASE}/results/${TISSUE}_r2_0.01_p_0.05_models.txt.gz"
echo "  Coefficients:  ${OUTPUT_BASE}/results/${TISSUE}_r2_0.01_p_0.05_betas.txt.gz"
echo "  Plot:          ${OUTPUT_BASE}/plots/${TISSUE}_R2_holdout_0.01_0.05.pdf"
echo "========================================="
```

---

## Citation

If you use these methods to build new models, please cite:

```
[Citation information will be added upon publication]
```

---

## Need Help?

- Review the [FAQ](../docs/FAQ.md)
- Check existing [GitHub Issues](https://github.com/rsigner/BMI_dynamic_models/issues)
- Open a new issue with:
  - Description of your data
  - Error messages
  - Steps you've already tried

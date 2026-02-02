# Pre-Trained BMI-Dynamic Models

This directory contains pre-trained gene expression prediction models for 8 GTEx v8 tissues.

**Total size**: ~16 MB

- Model coefficients: ~9.5 MB
- SNP mapping files: ~6 MB
- Feature data: ~0.5 MB

------

## Model Files

### Betas File Format

Each tissue model file (`{tissue}_r2_0.01_p_0.05_betas.txt.gz`) contains coefficients for all genes:

**Columns**:

- `name`: Short predictor name (SNP ID or "BMI")
- `Beta`: Model coefficient (effect size)
- `predictor_type`: Type of predictor ("main" or "interaction")
- `Predictor`: Full predictor name (includes interaction notation like "SNP*BMI")
- `Gene`: Ensembl gene ID with version
- `Description`: Gene symbol
- `Chr`: Chromosome
- `start`: Gene start position (hg38)
- `end`: Gene end position (hg38)
- `gene_type`: Gene biotype (protein_coding, lincRNA, etc.)
- `tissue`: Tissue name
- `n_predictors`: Total number of predictors for this gene

Predictor types**:

- **SNP main effect**: `predictor_type = "main"`, `name` is GTEx variant ID
- **BMI main effect**: `predictor_type = "main"`, `name = "BMI"`
- **SNP × BMI interaction**: `predictor_type = "interaction"`, `Predictor` contains "*BMI"

### Model Summary File Format

Model summary files (`{tissue}_r2_0.01_p_0.05_models.txt.gz`) contain performance metrics for each gene:

**Columns**:

- `Gene`: Ensembl gene ID with version
- `Tissue`: Tissue name
- `Rsq`: R-squared on crossfold validation set
- `p_Rsq`: P-value for crossfold R-squared
- `CV_Lambda`: Optimal lambda from cross-validation
- `CV_MStEr`: Mean standard error from CV
- `CV_MSE`: Mean squared error from CV
- `R2_rep`: 4th fold hold out R-squared
- `p_rep`: 4th fold hold out p-value
- `n_predictors`: Total number of predictors in model
- `Description`: Gene symbol
- `Chr`: Chromosome
- `start`: Gene start position (hg38)
- `end`: Gene end position (hg38)
- `gene_type`: Gene biotype

### Selection Criteria

Models included only genes meeting these criteria:

- **R² > 0.01, p-value < 0.05** on cross fold and hold out analyses
- At least one significant predictor



------

## Feature Data Files

Gene annotation files for each tissue containing:

**Columns**:

- `Name`: Ensembl gene ID with version (matches model files)
- `Chr`: Chromosome
- `start`: Gene start position (hg38)
- `end`: Gene end position (hg38)
- `Description`: Gene symbol
- `gene_type`: Gene biotype (protein_coding, lincRNA, etc.)

**Example**:

```
Name                   Chr  start      end        Description  gene_type
ENSG00000000003.15     X    100627108  100639991  TSPAN6      protein_coding
```

**Usage**:

- Define cis-windows when building models
- Annotate results with gene names and locations
- Filter by gene type for analyses

------

## Mapping Files

### SNP Mapping File

**File**: ``models/mapping/{tissue}_model_variants.txt.gz` (~6 MB)

Maps all SNPs used across all models to multiple ID formats:

```
CHR  pos_37   ref  alt  variant_id_hg38              rs_id        dbSNP151_GRCh38p7  tissue
1    14677    G    A    chr1_14677_G_A_b38          rs2123123    ...                Adipose-Subcutaneous
```

**Columns**:

- `CHR`: Chromosome number
- `pos_37`: Position in hg19/GRCh37 (for liftover)
- `ref`: Reference allele
- `alt`: Alternate allele (what models count!)
- `variant_id_hg38`: GTEx format ID (hg38 coordinates)
- `rs_id`: rsID from dbSNP
- `dbSNP151_GRCh38p7`: dbSNP build
- `tissue`: Which tissue(s) use this SNP

**Usage**:

- Convert biobank SNP IDs to GTEx format
- Liftover between genome builds
- Extract SNPs for specific tissues

See `mapping/README.md` for detailed usage instructions.

------

## Usage

### Using with Prediction Script

The prediction script automatically loads the correct model:

```bash
Rscript apply_models/predict_expression.R \
    --tissue Adipose-Subcutaneous \
    --chr 1 \
    --genotype-file genotypes.raw.gz \
    --bmi-file bmi.txt \
    --output-dir results/
```

**Note**: The prediction script uses the **betas file** (`*_betas.txt.gz`), not the models summary file.

------

## Model Training Details

**Algorithm**: GLInternet (group-lasso interaction network)

**Software versions**:

- R 4.0+
- glinternet package
- GTEx v8 data

See [`docs/MODEL_DETAILS.md`](https://claude.ai/docs/MODEL_DETAILS.md) for technical details.

------

## Citation

If you use these pre-trained models, please cite:

```
[Citation information will be added upon publication]
```

------

## Questions?

- See [`apply_models/README.md`](https://claude.ai/apply_models/README.md) for usage instructions
- Open an issue on GitHub for support
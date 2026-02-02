# BMI-Dynamic Models

Gene expression prediction models with BMI-gene interactions using glinternet.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

------

## Overview

This repository provides pre-trained models for **predicting BMI-dynamic tissue-specific gene expression** from genotype data, incorporating BMI (Body Mass Index) as a dynamic variable that modulates genetic effects from the Signer et al. 2026 publication.

**Key features**:

- Pre-trained models for 8 GTEx v8 tissues:  Adipose-Subcutaneous, Brain-Cerebellum, Brain-Cortex, Brain-Nucleus Accumbens, Colon-Transverse, Esophagus-Mucosa, Esophagus-Muscularis, Nerve-Tibial

------

## Workflows

### Apply Pre-Built BMI-Dynamic Models

**You have**: DNA data from a biobank (UK Biobank, etc.) **You want** **to** : Predict BMI-dynamic gene expression and test associations with disease in a biobank

👉 **Go to**: [`apply_models/README.md`]

### Build New Models

**You have**: Paired DNA + RNA data (like GTEx) **You want**: Train new BMI-dynamic models on your data using glinternet

👉 **Go to**: [`build_models/README.md`]

### View pre-built models

👉 **Go to**: [`models/README.md`]

------

## Installation

### Requirements

**R (required)**:

- R ≥ 4.0
- Packages: `tidyverse`, `data.table`, `glinternet`, `optparse`

**Python (for building models only)**:

- Python ≥ 3.7
- Packages: `pandas`, `numpy`, `scikit-learn`

**Other tools (for genotype preprocessing)**:

- PLINK2 (recommended for applying models)

### Install

```bash
# Clone repository (includes pre-trained models)
git clone https://github.com/rsigner/BMI_dynamic_models.git
cd BMI_dynamic_models

# Install R dependencies
R -e "install.packages(c('tidyverse', 'data.table', 'glinternet', 'optparse'))"

# Install Python dependencies (if building models)
pip install pandas numpy scikit-learn
```

------

## Use Cases

### 1. Common Disease Association Studies

Predict BMI-aware gene expression expression from SNPs in large scale biobanke for logistic and linear association testing 

```bash
# Predict expression for cases and controls
Rscript apply_models/predict_expression.R \
    --tissue Adipose-Subcutaneous \
    --chr 1 \
    --part 1 \
    --genotype-file genotypes.raw.gz \
    --bmi-file bmi.txt \
    --output-dir results/

# Test association between predicted expression and disease status
Rscript apply_models/test_association.R \
    --tissue Adipose-Subcutaneous \
    --chr 1 \
    --part 1 \
    --expression-dir results/ \
    --phenotype-file disease_status.txt \
    --output-dir associations/

# Assign genome-wide significance
Rscript apply_models/assign_significance.R \
    --tissue Adipose-Subcutaneous \
    --input-dir associations/ \
    --output-file final_results.txt.gz
```

### 2. Cross-Tissue Analysis

Compare BMI effects across multiple tissues:

```bash
TISSUES=("Adipose-Subcutaneous" "Brain-Cortex" "Colon-Transverse")
for TISSUE in "${TISSUES[@]}"; do
    # Predict expression for each tissue
    for CHR in {1..22}; do
        for PART in {1..10}; do
            Rscript apply_models/predict_expression.R \
                --tissue ${TISSUE} \
                --chr ${CHR} \
                --part ${PART} \
                --genotype-file genotypes_chr${CHR}.raw.gz \
                --bmi-file bmi.txt \
                --output-dir results/${TISSUE}/
        done
    done
done
```

### Model Features

- **Algorithm**: GLInternet (group-lasso interaction network)
- **Available Predictors**:
  - SNP main effects
  - BMI main effect
  - SNP × BMI interaction effects
- **Selection**: R² > 0.01, p < 0.05 on held-out validation set

### Performance

- **Training**: 3-fold cross-validation
- **Validation**: Independent 4th fold
- **Metrics**: R², prediction p-value

------

## Citation

If you use these models in your research, please cite:

```
[Citation information will be added upon publication]
```

------

### Additional Resources

- [`INSTALLATION.md`](https://claude.ai/chat/INSTALLATION.md) - Installation guide
- [`docs/FAQ.md`](https://claude.ai/chat/docs/FAQ.md) - Frequently asked questions
- [`CONTRIBUTING.md`](https://claude.ai/chat/CONTRIBUTING.md) - How to contribute

------

## License

This project is licensed under the MIT License - see the [LICENSE](https://claude.ai/chat/LICENSE) file for details.

------

## Contact

- **Issues**: [GitHub Issues](https://github.com/rsigner/BMI_dynamic_models/issues)

------

## Acknowledgments

- GTEx Consortium for providing data
- LMH acknowledges funding from NIMH (R01MH124839, R01MH118278, R01MH125938, RM1MH132648, R01MH136149), NIEHS (R01ES033630), and the Department of Defense (TP220451). This work was supported in part through the computational and data resources and staff expertise provided by Scientific Computing and Data at the Icahn School of Medicine at Mount Sinai and supported by the Clinical and Translational Science Awards (CTSA) grant UL1TR004419 from the National Center for Advancing Translational Sciences. Research reported in this publication was also supported by the Office of Research Infrastructure of the National Institutes of Health under award number S10OD026880 and S10OD030463. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.

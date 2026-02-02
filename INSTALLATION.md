# Installation Guide

## Requirements

### R (Required)

**Version**: R ≥ 4.0

**Required packages**:

```r
install.packages(c(
  "tidyverse",
  "data.table",
  "glinternet",
  "optparse"
))
```

**Check installation**:

```r
# Check R version
R.version.string

# Check packages
library(tidyverse)
library(data.table)
library(glinternet)
library(optparse)
```

### Python (Required for building models only)

**Version**: Python ≥ 3.7

**Required packages**:

```bash
pip install pandas numpy scikit-learn
```

**Check installation**:

```bash
python --version
python -c "import pandas, numpy, sklearn; print('All packages installed')"
```

### PLINK2 (Recommended for data preprocessing)

**For genotype preprocessing** (applying models):

**Install**:

```bash
# Download from https://www.cog-genomics.org/plink/2.0/
wget https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_x86_64_latest.zip
unzip plink2_linux_x86_64_latest.zip
chmod +x plink2

# Add to PATH
export PATH=$PATH:/path/to/plink2/directory
```

**Check installation**:

```bash
plink2 --version
```

------

## Installation Methods

### Method 1: Clone from GitHub (Recommended)

```bash
# Clone repository
git clone https://github.com/rsigner/BMI_dynamic_models.git
cd BMI_dynamic_models

# Install R dependencies
Rscript -e "install.packages(c('tidyverse', 'data.table', 'glinternet', 'optparse'), repos='https://cloud.r-project.org/')"

# Install Python dependencies (if building models)
pip install pandas numpy scikit-learn

# Test installation
Rscript apply_models/predict_expression.R --help
python build_models/create_folds.py --help
```

### Method 2: Download ZIP

```bash
# Download ZIP from GitHub
# Extract to desired location
unzip BMI_dynamic_models-main.zip
cd BMI_dynamic_models-main

# Install dependencies (same as Method 1)
```

------

## Verify Installation

### Test R Scripts

```bash
# Test prediction script
Rscript apply_models/predict_expression.R --help

# Should show usage information
```

### Test Python Scripts

```bash
# Test fold creation script
python build_models/create_folds.py --help

# Should show usage information
```

### Check Pre-trained Models

```bash
# Check if model files are present
ls models/tissues/

# Should see tissue-specific model files (if models are included)
```

------

## Platform-Specific Instructions

### Linux

Most cluster computing environments use Linux. The installation above should work directly.

**Module system** (if available on your cluster):

```bash
module load R/4.2.0
module load python/3.9
```

### macOS

**Install R**:

- Download from https://cran.r-project.org/bin/macosx/
- Or use Homebrew: `brew install r`

**Install Python**:

- Download from https://www.python.org/downloads/
- Or use Homebrew: `brew install python`

**Install R packages** (same as above):

```r
install.packages(c("tidyverse", "data.table", "glinternet", "optparse"))
```

### Windows

**Install R**:

- Download from https://cran.r-project.org/bin/windows/base/

**Install Python**:

- Download from https://www.python.org/downloads/

**Install R packages** (open R and run):

```r
install.packages(c("tidyverse", "data.table", "glinternet", "optparse"))
```

**Note**: Some scripts use Unix-style paths. You may need to modify file paths for Windows.

------

## HPC Cluster Setup

### Loading Modules

Most HPC systems use environment modules:

```bash
# Check available modules
module avail

# Load R and Python
module load R/4.2.0
module load python/3.9

# Load additional dependencies if needed
module load gcc/11.2.0  # For compiling R packages
```

### Installing R Packages on Cluster

```bash
# Start R
R

# Set personal library location
.libPaths("/home/your_username/R/library")

# Install packages
install.packages(c("tidyverse", "data.table", "glinternet", "optparse"),
                 repos='https://cloud.r-project.org/',
                 lib="/home/your_username/R/library")
```

### Installing Python Packages on Cluster

```bash
# Create virtual environment
python -m venv ~/bmi_models_env

# Activate environment
source ~/bmi_models_env/bin/activate

# Install packages
pip install pandas numpy scikit-learn

# Add to your job scripts:
# source ~/bmi_models_env/bin/activate
```

------

## Troubleshooting

### R Package Installation Fails

**Problem**: Package compilation errors

**Solution**:

```bash
# Install development tools (Linux)
sudo apt-get install build-essential libcurl4-openssl-dev libssl-dev libxml2-dev

# Or request admin to install on cluster
```

### glinternet Package Issues

**Problem**: glinternet fails to install

**Solution**:

```r
# Install from source
install.packages("glinternet", type="source")

# Or install dependencies first
install.packages("Rcpp")
install.packages("RcppEigen")
install.packages("glinternet")
```

### Python sklearn Not Found

**Problem**: `ModuleNotFoundError: No module named 'sklearn'`

**Solution**:

```bash
pip install --upgrade scikit-learn
```

### Permission Denied

**Problem**: Cannot write to R library directory

**Solution**:

```r
# Create personal library
dir.create("~/R/library", recursive=TRUE)
.libPaths("~/R/library")

# Then install packages
install.packages("tidyverse", lib="~/R/library")
```

### PLINK2 Not in PATH

**Problem**: `command not found: plink2`

**Solution**:

```bash
# Add to ~/.bashrc or ~/.bash_profile
export PATH=$PATH:/path/to/plink2

# Reload
source ~/.bashrc
```

------

## System Requirements

### Minimum Requirements

**For applying models:**

- CPU: 1 core
- RAM: 8 GB (16 GB recommended)
- Disk: 50 GB for input/output files

**For building models:**

- CPU: 1 core per job
- RAM: 32 GB (64 GB for chromosome 1, 6)
- Disk: 100 GB for model files

### Recommended Setup

**HPC cluster with**:

- Job scheduler (SLURM, PBS, SGE)
- 64+ GB RAM per node
- Ability to submit array jobs
- Shared filesystem for data

**Software versions tested**:

- R 4.2.0
- Python 3.9
- PLINK 2.0

------

## Optional Tools

### For Analysis and Visualization

```r
# Additional R packages for analysis
install.packages(c(
  "ggplot2",     # Plotting (included in tidyverse)
  "gridExtra",   # Multi-panel plots
  "viridis",     # Color scales
  "cowplot"      # Publication-ready plots
))
```

### For Pipeline Management

```bash
# Snakemake (for workflow management)
pip install snakemake

# Nextflow (alternative)
curl -s https://get.nextflow.io | bash
```

------

## Updating

### Update from GitHub

```bash
cd BMI_dynamic_models
git pull origin main
```

### Update R Packages

```r
update.packages(ask=FALSE)
```

### Update Python Packages

```bash
pip install --upgrade pandas numpy scikit-learn
```

------

## Getting Help

**Installation issues**:

1. Check this guide
2. Review error messages
3. Search existing GitHub issues
4. Open new issue with:
   - Your operating system
   - R/Python versions
   - Complete error message
   - Steps you've tried

**Contact**: [GitHub Issues](https://github.com/rsigner/BMI_dynamic_models/issues)
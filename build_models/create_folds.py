#!/usr/bin/env python3
"""
Create Cross-Validation Folds Stratified by BMI

This script creates stratified k-fold cross-validation splits based on BMI 
categories to ensure balanced representation across folds for model training.

Author: Rebecca Signer
Date: 01/31/26
"""

import argparse
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold
from pathlib import Path
import sys


def categorize_bmi(bmi_value):
    """
    Categorize BMI into weight categories based on standard CDC definitions.
    
    Parameters
    ----------
    bmi_value : float
        BMI value
    
    Returns
    -------
    str
        Weight category: 'Underweight', 'Normal Weight', 'Overweight', or 'Obese'
    """
    if pd.isna(bmi_value):
        return None
    elif bmi_value < 18.5:
        return 'Underweight'
    elif bmi_value < 25:
        return 'Normal Weight'
    elif bmi_value < 30:
        return 'Overweight'
    else:
        return 'Obese'


def create_folds(pheno_path, output_dir, tissue_name, n_folds=4, random_state=42):
    """
    Create stratified cross-validation folds based on BMI categories.
    
    Parameters
    ----------
    pheno_path : str
        Path to phenotype data file (tab-separated with SUBJID and BMI columns)
    output_dir : str
        Directory where output file will be saved
    tissue_name : str
        Tissue name for organizing outputs
    n_folds : int, default=4
        Number of cross-validation folds
    random_state : int, default=42
        Random seed for reproducibility
    
    Returns
    -------
    None
        Saves output file: {tissue}_CV_test_ID_map_weightcategory.txt.gz
    """
    
    print("\n" + "="*80)
    print(" Creating Cross-Validation Folds")
    print("="*80)
    print(f"  Tissue:         {tissue_name}")
    print(f"  Phenotype file: {pheno_path}")
    print(f"  Output dir:     {output_dir}")
    print(f"  Number of folds: {n_folds}")
    print("="*80 + "\n")
    
    # Read phenotype file
    print("Loading phenotype data...")
    try:
        phenodt_raw = pd.read_csv(pheno_path, header=0, sep="\t")
    except FileNotFoundError:
        print(f"Error: Phenotype file not found: {pheno_path}")
        sys.exit(1)
    
    # Check required columns
    if 'SUBJID' not in phenodt_raw.columns or 'BMI' not in phenodt_raw.columns:
        print("Error: Phenotype file must contain 'SUBJID' and 'BMI' columns")
        print(f"  Found columns: {', '.join(phenodt_raw.columns)}")
        sys.exit(1)
    
    # Subset to required columns
    phen_here = phenodt_raw[['SUBJID', 'BMI']].copy()
    
    print(f"  Loaded: {len(phen_here)} individuals")
    
    # Create weight category column
    print("\nCategorizing BMI...")
    phen_here['Weight_Category'] = phen_here['BMI'].apply(categorize_bmi)
    
    # Remove any rows with missing BMI
    n_before = len(phen_here)
    phen_here = phen_here.dropna(subset=['Weight_Category'])
    n_after = len(phen_here)
    
    if n_before != n_after:
        print(f"  Warning: Removed {n_before - n_after} individuals with missing BMI")
    
    # Print category distribution
    print("\n  BMI category distribution:")
    for category, count in phen_here['Weight_Category'].value_counts().sort_index().items():
        pct = 100 * count / len(phen_here)
        print(f"    {category:15s}: {count:5d} ({pct:5.1f}%)")
    
    # Reset index
    phen_here_reset = phen_here.reset_index(drop=True)
    
    # Create stratified k-fold splits
    print(f"\nCreating {n_folds}-fold stratified splits...")
    stratifiedfold = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=random_state)
    
    # Define target column for stratification
    target = phen_here_reset['Weight_Category']
    
    # Loop through folds
    fold_here = 1
    test_all = []
    for train_index, test_index in stratifiedfold.split(phen_here_reset, target):
        print(f"  Creating fold {fold_here}...")
        
        # Get test set for this fold
        test = phen_here_reset.loc[test_index].copy()
        test['Fold_in_test'] = fold_here
        test_all.append(test)
        
        fold_here += 1
    
    # Combine all folds
    test_final = pd.concat(test_all)
    
    # Verify fold balance
    print("\n  Fold sizes:")
    for fold in range(1, n_folds + 1):
        n_fold = len(test_final[test_final['Fold_in_test'] == fold])
        pct = 100 * n_fold / len(test_final)
        print(f"    Fold {fold}: {n_fold:5d} ({pct:5.1f}%)")
    
    # Create output directory if needed
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Save output
    output_file = f"{output_dir}/{tissue_name}_CV_test_ID_map_weightcategory.txt.gz"
    test_final.to_csv(output_file, sep='\t', header=True, index=False, compression="gzip")
    
    print(f"\nSaved: {output_file}")
    print("\nDone!\n")


def main():
    """Parse arguments and run fold creation."""
    parser = argparse.ArgumentParser(
        description="Create stratified cross-validation folds based on BMI categories",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python create_folds.py \\
      --tissue Brain-Cortex \\
      --phenotype-file data/Brain-Cortex_phenoData.txt.gz \\
      --output-dir output/folds/

Output:
  {tissue}_CV_test_ID_map_weightcategory.txt.gz
  
  Columns:
    - SUBJID: Subject identifiers
    - BMI: Body Mass Index values
    - Weight_Category: BMI category (Underweight/Normal Weight/Overweight/Obese)
    - Fold_in_test: Fold assignment (1 to n_folds)
        """
    )
    
    parser.add_argument(
        '--tissue',
        type=str,
        required=True,
        help='Tissue name (e.g., Brain-Cortex, Adipose-Subcutaneous)'
    )
    
    parser.add_argument(
        '--phenotype-file',
        type=str,
        required=True,
        help='Path to phenotype file (tab-separated with SUBJID and BMI columns)'
    )
    
    parser.add_argument(
        '--output-dir',
        type=str,
        required=True,
        help='Output directory for fold assignments'
    )
    
    parser.add_argument(
        '--n-folds',
        type=int,
        default=4,
        help='Number of cross-validation folds (default: 4)'
    )
    
    parser.add_argument(
        '--random-state',
        type=int,
        default=42,
        help='Random seed for reproducibility (default: 42)'
    )
    
    args = parser.parse_args()
    
    # Run fold creation
    create_folds(
        pheno_path=args.phenotype_file,
        output_dir=args.output_dir,
        tissue_name=args.tissue,
        n_folds=args.n_folds,
        random_state=args.random_state
    )


if __name__ == '__main__':
    main()
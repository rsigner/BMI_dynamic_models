# Contributing to BMI-Dynamic Models

Thank you for your interest in contributing! This document provides guidelines for contributing to the BMI-dynamic models repository.

## Ways to Contribute

### 1. Report Bugs

Found a bug? Please open a [GitHub Issue](https://github.com/rsigner/BMI_dynamic_models/issues) with:

- **Clear title** describing the problem
- **Steps to reproduce** the issue
- **Expected behavior** vs actual behavior
- **Error messages** (complete text)
- **System information**:
  - Operating system
  - R version (`R.version.string`)
  - Python version (`python --version`)
  - Package versions
- **Example data** (if possible, small reproducible example)

**Example**:

```
Title: predict_expression.R fails with "column not found" error

Description:
When running predict_expression.R on chromosome 22, I get:
Error in select: column 'chr22_12345_A_G_b38' not found

Steps to reproduce:
1. Run: Rscript predict_expression.R --chr 22 --part 1 ...
2. See error above

System:
- Ubuntu 20.04
- R version 4.2.0
- tidyverse 2.0.0

Expected: Script should run without error
Actual: Column not found error
```

### 2. Suggest Features

Have an idea for improvement? Open an issue with:

- **Feature description**: What do you want to add?
- **Use case**: Why is this useful?
- **Example**: How would it work?
- **Alternatives**: Other ways to achieve this?

### 3. Improve Documentation

Documentation contributions are very welcome! You can:

- Fix typos or unclear explanations
- Add examples
- Improve README organization
- Add FAQ entries
- Create tutorials

### 4. Submit Code

Ready to contribute code? Great! Please follow these guidelines:

#### Before You Start

1. **Check existing issues** - Someone might already be working on this
2. **Open an issue** - Discuss your idea first
3. **Fork the repository** - Create your own copy
4. **Create a branch** - Don't work on `main`

#### Coding Standards

**R Code**:

- Use `tidyverse` style guide: https://style.tidyverse.org/
- Include comments for complex logic
- Use meaningful variable names
- Add error checking and validation

**Python Code**:

- Follow PEP 8: https://pep8.org/
- Use docstrings for functions
- Type hints appreciated
- Include error handling

**General**:

- Keep functions focused and modular
- Add command-line arguments instead of hardcoding
- Include progress messages for long-running operations
- Test your code thoroughly

#### Pull Request Process

1. **Fork** the repository

2. **Clone** your fork

3. **Create** a new branch:

   ```bash
   git checkout -b feature/your-feature-name
   ```

4. **Make your changes**:

   - Write clean, documented code
   - Add tests if applicable
   - Update documentation

5. **Test** your changes:

   - Run existing tests
   - Test with example data
   - Check for errors

6. **Commit** with clear messages:

   ```bash
   git add .
   git commit -m "Add feature: brief description
   
   - Detailed point 1
   - Detailed point 2
   - Fixes #123"
   ```

7. **Push** to your fork:

   ```bash
   git push origin feature/your-feature-name
   ```

8. **Open** a Pull Request:

   - Describe what you changed and why
   - Reference related issues
   - Explain how you tested it

#### Pull Request Template

```markdown
## Description
Brief description of changes

## Motivation
Why is this change needed?

## Changes Made
- Change 1
- Change 2

## Testing
How did you test this?

## Checklist
- [ ] Code follows style guidelines
- [ ] Documentation updated
- [ ] Tests added/updated
- [ ] No breaking changes (or documented)
```

------

## Development Setup

### 1. Fork and Clone

```bash
# Fork on GitHub, then clone your fork
git clone https://github.com/YOUR_USERNAME/BMI_dynamic_models.git
cd BMI_dynamic_models

# Add upstream remote
git remote add upstream https://github.com/rsigner/BMI_dynamic_models.git
```

### 2. Create Development Branch

```bash
# Update main
git checkout main
git pull upstream main

# Create feature branch
git checkout -b feature/my-new-feature
```

### 3. Make Changes

Edit files, add features, fix bugs...

### 4. Test Changes

```bash
# Test R scripts
Rscript apply_models/predict_expression.R --help

# Test Python scripts
python build_models/create_folds.py --help

# Run with small test data
```

### 5. Commit and Push

```bash
# Check what changed
git status
git diff

# Stage changes
git add file1.R file2.py

# Commit with message
git commit -m "Add feature X

- Detailed description
- Why this is needed"

# Push to your fork
git push origin feature/my-new-feature
```

### 6. Open Pull Request

Go to GitHub and open a pull request from your branch to the main repository.

------

## Code Review Process

1. **Maintainers review** your pull request
2. **Feedback** may be requested
3. **Iterate** based on feedback
4. **Approval** when ready
5. **Merge** into main branch

**Typical timeline**: 1-2 weeks for review

------

## Reporting Security Issues

**Do not** open public issues for security vulnerabilities.

Instead:

- Email: [maintainer email]
- Describe the vulnerability
- Include steps to reproduce
- Allow time for fix before public disclosure

------

## Style Guidelines

### R Scripts

```r
# Good
calculate_prediction <- function(genotype_matrix, bmi_values) {
  # Validate inputs
  if (nrow(genotype_matrix) != length(bmi_values)) {
    stop("Genotype matrix and BMI values must have same length")
  }
  
  # Calculate predictions
  predictions <- genotype_matrix %*% beta_coefficients
  
  return(predictions)
}

# Bad
calc <- function(g,b) {
  g%*%beta
}
```

### Python Scripts

```python
# Good
def create_folds(phenotype_data: pd.DataFrame, 
                 n_folds: int = 4) -> pd.DataFrame:
    """
    Create stratified cross-validation folds.
    
    Args:
        phenotype_data: DataFrame with SUBJID and BMI columns
        n_folds: Number of folds to create
        
    Returns:
        DataFrame with fold assignments
    """
    # Implementation
    pass

# Bad
def cf(p,n):
    # no docstring
    pass
```

### Commit Messages

```bash
# Good
git commit -m "Fix: Handle missing BMI values in fold creation

- Add NA checking before stratification
- Return informative error message
- Add test case for missing values
- Fixes #123"

# Bad
git commit -m "fix bug"
```

------

## Adding New Features

### Before Starting

1. Open an issue to discuss
2. Get feedback from maintainers
3. Agree on approach

### During Development

1. Keep changes focused
2. Update documentation
3. Add tests
4. Follow coding standards

### When Submitting

1. Describe use case clearly
2. Provide examples
3. Document any breaking changes
4. Update README if needed

------

## Questions?

- **General questions**: Open a GitHub issue
- **Development help**: Comment on your pull request
- **Private questions**: Email [maintainer email]

------

## Recognition

Contributors will be:

- Listed in CONTRIBUTORS.md
- Acknowledged in release notes
- Credited in relevant publications (for major contributions)

------

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

------

Thank you for contributing to BMI-dynamic models! 🎉
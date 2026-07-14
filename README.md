# Compact serum miRNA panels for colorectal cancer classification

This repository contains the R and Python code supplied for the study of serum miRNA biomarkers for colorectal cancer detection and stage-based risk stratification using XGBoost.

The code is organized as simple, sequential scripts and Jupyter notebooks. The analytical methods and model definitions were not redesigned. Variable names, comments, file paths, notebook structure, and unused exploratory code were cleaned for readability.

## Repository structure

```text
R/
  01_data_preparation_DE_LASSO_enrichment.R
notebooks/
  00_train_test_split.ipynb
  01_crc_vs_non_cancer_model.ipynb
  02_risk_stratification_model.ipynb
  03_selected_mirna_expression.ipynb
  04_single_mirna_roc.ipynb
data/
  raw/
  external/
  processed/
results/
  figures/
  tables/
  intermediate/
environment/
  python_requirements.txt
  R_packages.txt
  report_python_versions.py
  report_R_versions.R
```

## Analysis workflow

1. Retrieve sample annotations for GSE211692.
2. Combine the CRC and non-cancer expression matrices.
3. Create stratified 75:25 training and hold-out test sets using `random_state=1`.
4. Quantile-normalize the training expression matrix for differential-expression and LASSO analyses.
5. Perform limma differential-expression analysis.
6. Perform LASSO feature selection with `lambda.1se` and 10-fold cross-validation.
7. Prepare the selected-miRNA matrices for Python.
8. Train the five-miRNA CRC-versus-non-cancer XGBoost model.
9. Rank the 24 candidate risk miRNAs and select the top 12 using five-fold macro-F1 cross-validation.
10. Train and evaluate the weighted 12-miRNA risk-stratification model.
11. Generate ROC, precision–recall, individual-miRNA ROC, expression, and SHAP outputs.

## Software

The supplied notebooks were created with Python 3.11.9. The manuscript reports R 4.3.2. Package names are listed in the `environment` directory. Exact package versions should be captured from the original analysis environment before the public release.

## Required data

See `data/README.md` for filenames and locations. Large expression matrices should not be committed to GitHub.

## Running the code

Run the R script from the repository root:

```bash
Rscript R/01_data_preparation_DE_LASSO_enrichment.R
```

Launch Jupyter from the `notebooks` directory so the relative paths resolve as written:

```bash
cd notebooks
jupyter notebook
```

Run the notebooks in this order:

1. `00_train_test_split.ipynb`
2. `01_crc_vs_non_cancer_model.ipynb`
3. `02_risk_stratification_model.ipynb`
4. `03_selected_mirna_expression.ipynb`
5. `04_single_mirna_roc.ipynb`

## Model definitions

The exact model settings present in the supplied code are documented in `MODEL_HYPERPARAMETERS.md`.

## Important status note

The attached files did not contain every step required for an end-to-end reconstruction. The missing items are listed in `CODE_AUDIT.md`. These gaps should be resolved before the repository is made public or cited in the reviewer response.

## Data source

The expression data and sample annotations are from GEO accession GSE211692.

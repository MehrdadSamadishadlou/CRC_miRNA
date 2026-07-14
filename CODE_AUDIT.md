# Code audit and outstanding items

This repository is a cleaned version of the files supplied by the author. No new analytical functions or model architectures were introduced.

## Removed material

- Hard-coded paths from the author's computer.
- Unused Python and R imports.
- Duplicate feature-subset scoring code.
- The Gaussian-process hyperparameter-search section because its returned parameters were not used by the final risk model.
- Abandoned or undefined notebook cells.
- Classification-report loops in the single-miRNA notebook that did not contribute to a reported figure or final result.
- PCA plotting code.
- Before/after-normalization boxplots.
- Enrichment plotting code; enrichment result-table generation was retained.
- Stored notebook outputs.

## Execution fixes

- Relative repository paths replaced laptop-specific paths.
- Notebook code was separated into small sequential cells.
- Plot axes were explicitly created in the single-miRNA ROC notebook because the supplied cells referred to an undefined `ax` object.
- Output directories and filenames were made explicit.

## Resolved item

- The original train/test splitter has now been supplied and cleaned as
  `notebooks/00_train_test_split.ipynb`.
- It uses `train_test_split` with `test_size=0.25`, `random_state=1`, and
  stratification by the class labels.

## Missing code that must be supplied before public release

1. **Complete stage-specific differential-expression and LASSO workflow**
   - The supplied R script contains an explicit stage 0 differential-expression block.
   - It then reads a pre-existing stage 4 differential-expression file for LASSO.
   - Equivalent executable blocks for the overall CRC comparison and stages 1–4 were not supplied.
   - The all-CRC LASSO file required by enrichment is also read as an existing file rather than generated in the supplied script.

2. **TargetScanHuman–HCMDB analysis**
   - No code for target filtering, gene-set intersection, or table generation was included in the attached files.

3. **Exact package versions**
   - Python 3.11.9 is recorded in the notebook metadata and R 4.3.2 is reported in the manuscript.
   - Exact package versions were not present in the attached files.

## Method/code points requiring author confirmation

1. The diagnostic precision–recall plot in the supplied code uses probabilities obtained after fitting on the complete training set, whereas a separate later cell calculates out-of-fold probabilities. The code was preserved rather than silently changing which probabilities are plotted.
2. The risk-model feature-subset cross-validation does not pass the class-dependent sample weights. The weights are applied in the final model and in the later ROC/out-of-fold calculations.
3. The final risk XGBoost model does not specify `random_state`.
4. The single-miRNA models use XGBoost defaults without a fixed random seed.
5. The R script contains fixed sample counts such as 4,232 healthy training samples and 90 stage 4 samples. These should be checked against the final split.

These points do not necessarily mean that the results are incorrect, but they must be resolved or documented for a defensible reproducibility statement.
